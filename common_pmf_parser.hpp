#ifndef COMMON_PMF_PARSER_HPP
#define COMMON_PMF_PARSER_HPP

#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>

#include "common_tools.hpp"
#include "third_party/ndarray.hpp"
#include "third_party/ndarray_io.hpp"

// parsing pmf files
// usage:
//   // read NAMD formatted PMF file
//   auto a = Pmf<double>("file.pmf")
//   // read plain PMF file
//   // the last par (pbc) can be omitted
//   auto a = Pmf<double>("file.pmf",{-20,0},{0.2,0.1},{20,3},{0,0})
//   // write NAMD formmatted PMF file
//   a.writePmfFile("file2.pmf")
//   // get data
//   a[{-20,0}]
//   a.getPmfData()
//   a.getLowerboundary()
//   a.getUpperboundary()
//   a.getWidth()
//   a.getShape()
//   a.getDimension()
//   a.RCToInternal()
//   a.internalToRC()
//   // calculate RMSD between two pmfs
//   auto b = Pmf<double>("file2.pmf")
//   a.calcRMSD(b)
//   // return a region by a point in this region and a cutoff energy
//   auto r1 = a.determineRegion({0.2, 1.3}, 5)
//   auto r2 = a.determineRegion({-10, 3.4}, 5)
//   // calculate free-energy difference between two regions
//   a.calcFreeEnergyDiff(r1, r2, 298.15)
//

namespace pmf_parser {

    // pmf (T=double) or count (T=int) data
    template <typename T>
    class Pmf {
       public:
        // copy constructor
        Pmf(const Pmf& another_Pmf) {
            this->lowerboundary_ = another_Pmf.lowerboundary_;
            this->upperboundary_ = another_Pmf.upperboundary_;
            this->width_ = another_Pmf.width_;
            this->pbc_ = another_Pmf.pbc_;
            this->shape_ = another_Pmf.shape_;
            this->dimension_ = another_Pmf.dimension_;
            this->data_ = new ndarray::NdArray<T>(*(another_Pmf.data_));
        }

        // initialize the pmf using NAMD formatted file
        Pmf(const std::string& pmf_file) {
            static_assert(
                std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be a kind of number");

            std::ifstream read_file;
            read_file.open(pmf_file, std::ios::in);
            if (!read_file.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            // used for parsing
            std::string line;
            std::vector<std::string> splited_line;

            getline(read_file, line);
            if (!pystring::startswith(line, "#")) {
                std::cerr << "This is not an NAMD PMF file!" << std::endl;
                exit(1);
            }

            // dimension
            pystring::split(line, splited_line);
            this->dimension_ = std::stoi(splited_line[1]);

            // lb, width and ub
            this->lowerboundary_ = std::vector<double>(this->dimension_);
            this->upperboundary_ = std::vector<double>(this->dimension_);
            this->width_ = std::vector<double>(this->dimension_);
            this->pbc_ = std::vector<bool>(this->dimension_);
            this->shape_ = std::vector<int>(this->dimension_);
            for (int i = 0; i < this->dimension_; i++) {
                getline(read_file, line);
                pystring::split(line, splited_line);
                this->shape_[i] = std::stoi(splited_line[3]);
                this->width_[i] = std::stod(splited_line[2]);
                this->lowerboundary_[i] =
                    std::stod(splited_line[1]) + 0.5 * this->width_[i];
                this->upperboundary_[i] =
                    this->lowerboundary_[i] +
                    this->width_[i] * (this->shape_[i] - 1);
                this->pbc_[i] = bool(std::stoi(splited_line[4]));
            }

            // reading data
            this->data_ = new ndarray::NdArray<T>(this->shape_);
            std::vector<double> RCPosition(this->dimension_);
            while (getline(read_file, line)) {
                pystring::split(line, splited_line);
                if (splited_line.size() == 0) {
                    continue;
                }
                for (int i = 0; i < this->dimension_; i++) {
                    RCPosition[i] = std::stod(splited_line[i]);
                }
                (*(this->data_))[this->RCToInternal(RCPosition)] =
                    std::stod(splited_line[this->dimension_]);
            }

            read_file.close();
        }

        // initialize the pmf given lb, ub, width and pbc
        Pmf(const std::string& pmf_file,
            const std::vector<double>& lowerboundary,
            const std::vector<double>& width,
            const std::vector<double>& upperboundary,
            std::vector<bool> pbc = {}) {
            static_assert(
                std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be a kind of number");

            assert(lowerboundary.size() == width.size());
            assert(lowerboundary.size() == upperboundary.size());

            this->lowerboundary_ = lowerboundary;
            this->upperboundary_ = upperboundary;
            this->width_ = width;
            this->dimension_ = lowerboundary.size();

            if (pbc == std::vector<bool>{}) {
                this->pbc_ = std::vector<bool>(lowerboundary.size(), false);
            } else {
                this->pbc_ = pbc;
            }

            // the shape of internal data
            this->shape_ = std::vector<int>(this->dimension_, 0);
            for (int i = 0; i < this->dimension_; i++) {
                // +1 means the boundaries are included
                this->shape_[i] = int((upperboundary[i] - lowerboundary[i] +
                                       common_tools::kAccuracy) /
                                      width[i]) +
                                  1;
            }
            this->data_ = new ndarray::NdArray<T>(this->shape_);

            // read pmfFile into memory
            auto dat_format_data = ndarray::ReadDat(pmf_file, T(0));
            auto dat_shape = dat_format_data.GetShape();
            // convert it into internal data
            std::vector<double> RCPosition;
            for (int row = 0; row < dat_shape[0]; row++) {
                RCPosition = std::vector<double>(this->dimension_, 0);
                for (int col = 0; col < dat_shape[1]; col++) {
                    if (col < this->dimension_) {
                        RCPosition[col] = dat_format_data[{row, col}];
                    } else {
                        (*(this->data_))[this->RCToInternal(RCPosition)] =
                            T(dat_format_data[{row, col}]);
                        break;
                    }
                }
            }
        }

        // initialize an empty pmf (data = 0) given lb, ub, width and pbc
        Pmf(const std::vector<double>& lowerboundary,
            const std::vector<double>& width,
            const std::vector<double>& upperboundary,
            std::vector<bool> pbc = {}) {
            static_assert(
                std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be a kind of number");

            assert(lowerboundary.size() == width.size());
            assert(lowerboundary.size() == upperboundary.size());

            this->lowerboundary_ = lowerboundary;
            this->upperboundary_ = upperboundary;
            this->width_ = width;
            this->dimension_ = lowerboundary.size();

            if (pbc == std::vector<bool>{}) {
                this->pbc_ = std::vector<bool>(lowerboundary.size(), false);
            } else {
                this->pbc_ = pbc;
            }

            // the shape of internal data
            this->shape_ = std::vector<int>(this->dimension_, 0);
            for (int i = 0; i < this->dimension_; i++) {
                // +1 means the boundaries are included
                this->shape_[i] = int((upperboundary[i] - lowerboundary[i] +
                                       common_tools::kAccuracy) /
                                      width[i]) +
                                  1;
            }
            this->data_ = new ndarray::NdArray<T>(this->shape_);
        }

        // initialize the pmf using NAMD ifstream
        // used in parsing .hist files
        Pmf(std::ifstream& pmf_stream) {
            static_assert(
                std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be a kind of number");

            // used for parsing
            std::string line;
            std::vector<std::string> splited_line;

            getline(pmf_stream, line);
            if (!pystring::startswith(line, "#")) {
                std::cerr << "This is not an NAMD PMF stream!" << std::endl;
                exit(1);
            }

            // dimension
            pystring::split(line, splited_line);
            this->dimension_ = std::stoi(splited_line[1]);

            // number of bins
            int num_of_bins = 1;
            // lb, width and ub
            this->lowerboundary_ = std::vector<double>(this->dimension_);
            this->upperboundary_ = std::vector<double>(this->dimension_);
            this->width_ = std::vector<double>(this->dimension_);
            this->pbc_ = std::vector<bool>(this->dimension_);
            this->shape_ = std::vector<int>(this->dimension_);
            for (int i = 0; i < this->dimension_; i++) {
                getline(pmf_stream, line);
                pystring::split(line, splited_line);
                this->shape_[i] = std::stoi(splited_line[3]);
                num_of_bins *= shape_[i];
                this->width_[i] = std::stod(splited_line[2]);
                this->lowerboundary_[i] =
                    std::stod(splited_line[1]) + 0.5 * this->width_[i];
                this->upperboundary_[i] =
                    this->lowerboundary_[i] +
                    this->width_[i] * (this->shape_[i] - 1);
            }

            // reading data
            this->data_ = new ndarray::NdArray<T>(this->shape_);
            std::vector<double> RCPosition(this->dimension_);
            for (int i = 0; i < num_of_bins;) {
                getline(pmf_stream, line);
                pystring::split(line, splited_line);
                if (splited_line.size() == 0) {
                    continue;
                }
                for (int j = 0; j < this->dimension_; j++) {
                    RCPosition[j] = std::stod(splited_line[j]);
                }
                (*(this->data_))[this->RCToInternal(RCPosition)] =
                    std::stod(splited_line[this->dimension_]);

                i++;
            }
            // an empty line
            getline(pmf_stream, line);
        }

        // write internal data to a pmf file
        // in NAMD pmf format!
        // note: PBCs are not recorded! So they are zeroes!
        void WritePmfFile(const std::string& file) const {
            std::ofstream write_file;
            write_file.open(file, std::ios::out);
            if (!write_file.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            // the head of NAMD pmf file
            write_file << std::setw(2) << "# " << this->dimension_ << "\n";
            for (int i = 0; i < this->dimension_; i++) {
                write_file << "# " << std::setw(10)
                          << this->lowerboundary_[i] - 0.5 * this->width_[i]
                          << std::setw(10) << this->width_[i] << std::setw(10)
                          << this->shape_[i] << " " << int(this->pbc_[i])
                          << "\n";
            }
            write_file << "\n";

            // iterate over any dimension
            int n = 0;
            std::vector<int> loop_flag(this->dimension_, 0);
            while (n >= 0) {
                auto RC = this->InternalToRC(loop_flag);
                for (auto coor : RC) {
                    write_file << common_tools::Round(
                                     coor, common_tools::kDecimalAccuracy)
                              << " ";
                }
                write_file << (*(this->data_))[loop_flag] << "\n";

                // mimic an nD for loop
                n = this->dimension_ - 1;
                while (n >= 0) {
                    loop_flag[n] += 1;
                    if (loop_flag[n] > this->shape_[n] - 1) {
                        loop_flag[n] = 0;
                        n--;
                        write_file << "\n";
                    } else {
                        break;
                    }
                }
            }

            write_file.close();
        }

        // get the data (ndarray) of the pmf
        const ndarray::NdArray<T>& GetPmfData() const { return *(this->data_); }

        // set the internal data using an NdArray
        void SetPMFData(const ndarray::NdArray<T>& data) {
            assert(this->data_->GetShape() == data.GetShape());
            delete this->data_;
            this->data_ = new ndarray::NdArray<T>(data);
        }

        // get lowerboundary, upperboundary, width, shape and dimension
        const std::vector<double>& GetLowerboundary() const {
            return this->lowerboundary_;
        }

        const std::vector<double>& GetUpperboundary() const {
            return this->upperboundary_;
        }

        const std::vector<double>& GetWidth() const { return this->width_; }

        const std::vector<bool>& GetPbc() const { return this->pbc_; }

        const std::vector<int>& GetShape() const { return this->shape_; }

        int GetDimension() const { return this->dimension_; }

        // operator[] to get the desired item at a given RCPosition
        // one may note that the var type of operator[] is extremely important
        T& operator[](const std::vector<double>& rc_position) const {
            return (*(this->data_))[this->RCToInternal(rc_position)];
        }

        // get the desired item using internal RC
        // this will make parsing PMF easier
        T& operator[](const std::vector<int>& internal_position) const {
            return (*(this->data_))[internal_position];
        }

        // set the value corresponding a given RCPosition
        // void setData(const std::vector<double>& RCPosition, T value) {
        //    (*(this->data))[this->RCToInternal(RCPosition)] = value;
        //}

        // void setData(const std::vector<int>& internalPosition, T value) {
        //     (*(this->data))[internalPosition] = value;
        // }

        // convert external/real reaction coordinate into the internal
        // coordinate
        std::vector<int> RCToInternal(
            const std::vector<double>& RCPosition) const {
            assert(RCPosition.size() == this->dimension_);

            std::vector<int> internal_position(this->dimension_);
            for (int i = 0; i < this->dimension_; i++) {
                internal_position[i] =
                    int((RCPosition[i] - this->lowerboundary_[i] +
                         common_tools::kAccuracy) /
                        this->width_[i]);
            }
            return internal_position;
        }

        // convert internal coordinate into external/real reaction coordinate
        std::vector<double> InternalToRC(
            const std::vector<int>& internal_position) const {
            assert(internal_position.size() == this->dimension_);

            std::vector<double> RCPosition(this->dimension_);
            for (int i = 0; i < this->dimension_; i++) {
                RCPosition[i] = double((internal_position[i] * this->width_[i]) +
                                       this->lowerboundary_[i]);
            }
            return RCPosition;
        }

        // calculate RMSD between two pmfs
        double CalcRMSD(const Pmf<T>& another_pmf_data) const {
            assert(this->shape_ == another_pmf_data.GetShape());

            // number of bins
            int num_of_bins = 0;
            double sum = 0;
            // iterate over any dimension
            int n = 0;
            std::vector<int> loop_flag(this->dimension_, 0);
            while (n >= 0) {
                sum += pow(
                    (*(this->data_))[loop_flag] - another_pmf_data[loop_flag], 2);
                num_of_bins++;

                // mimic an nD for loop
                n = this->dimension_ - 1;
                while (n >= 0) {
                    loop_flag[n] += 1;
                    if (loop_flag[n] > this->shape_[n] - 1) {
                        loop_flag[n] = 0;
                        n--;
                    } else {
                        break;
                    }
                }
            }

            sum /= num_of_bins;
            return sqrt(sum);
        }

        // calculate free-energy difference between two regions
        // the region is provided as a set of points
        double CalcFreeEnergyDiff(
            const std::vector<std::vector<double> >& region1,
            const std::vector<std::vector<double> >& region2,
            double temperature) const {
            assert(region1.size() != 0);
            assert(region2.size() != 0);

            for (const auto& point : region1) {
                assert(point.size() == this->dimension_);
            }
            for (const auto& point : region2) {
                assert(point.size() == this->dimension_);
            }

            // sum of partition function
            double sum1 = 0;
            double sum2 = 0;

            // delta RC in integration
            double dRC = 1;
            for (const auto& item : this->width_) {
                dRC *= item;
            }

            for (const auto& point : region1) {
                sum1 += exp(-(*this)[point] /
                            (common_tools::kBoltzmann * temperature)) *
                        dRC;
            }
            for (const auto& point : region2) {
                sum2 += exp(-(*this)[point] /
                            (common_tools::kBoltzmann * temperature)) *
                        dRC;
            }

            return -common_tools::kBoltzmann * temperature * log(sum1 / sum2);
        }

        // determine a region based on a point in this region and a cutoff
        // energy
        std::vector<std::vector<double> > DetermineRegion(
            const std::vector<double>& point, double cutoff_energy) const {
            std::vector<std::vector<double> > region;
            std::vector<std::vector<int> > close_list;
            std::vector<std::vector<int> > open_list;

            if ((*this)[point] < cutoff_energy) {
                open_list.push_back(this->RCToInternal(point));
            } else {
                return {};
            }

            while (open_list.size() > 0) {
                auto temp_open_list = open_list;
                open_list.clear();
                for (auto item : temp_open_list) {
                    close_list.push_back(item);
                    for (int i = 0; i < this->dimension_; i++) {
                        // left side
                        auto p = item;
                        if (item[i] - 1 >= 0) {
                            p[i] = item[i] - 1;
                        } else {
                            if (this->pbc_[i]) {
                                p[i] = this->shape_[i] - 1;
                            } else {
                                p = {};
                            }
                        }

                        // right side
                        auto q = item;
                        if (item[i] + 1 < this->shape_[i]) {
                            q[i] = item[i] + 1;
                        } else {
                            if (this->pbc_[i]) {
                                q[i] = 0;
                            } else {
                                q = {};
                            }
                        }

                        if (p.size() != 0) {
                            if ((*this)[p] < cutoff_energy &&
                                !common_tools::VectorInVectorOfVector(
                                    p, close_list) &&
                                !common_tools::VectorInVectorOfVector(
                                    p, open_list)) {
                                open_list.push_back(p);
                            }
                        }
                        if (q.size() != 0) {
                            if ((*this)[q] < cutoff_energy &&
                                !common_tools::VectorInVectorOfVector(
                                    q, close_list) &&
                                !common_tools::VectorInVectorOfVector(
                                    q, open_list)) {
                                open_list.push_back(q);
                            }
                        }
                    }
                }
            }

            for (const auto& item : close_list) {
                region.push_back(this->InternalToRC(item));
            }

            return region;
        }

        ~Pmf() { delete this->data_; }

       private:
        // the data (free energy) of the pmf
        ndarray::NdArray<T>* data_;
        // lowerboundary, upperboundary, width, and dimension
        std::vector<double> lowerboundary_;
        std::vector<double> upperboundary_;
        std::vector<double> width_;
        // pbc
        std::vector<bool> pbc_;
        // the shape of internal data
        std::vector<int> shape_;
        int dimension_;
    };
}  // namespace pmf_parser

#endif
