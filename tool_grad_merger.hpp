#ifndef TOOL_GRAD_MERGER_HPP
#define TOOL_GRAD_MERGER_HPP

#include <cassert>
#include <iomanip>
#include <string>
#include <vector>

#include "common_tools.hpp"
#include "third_party/ndarray.hpp"
#include "third_party/ndarray_io.hpp"

// An minimal implementation to merge grad files
// We are not going to parse grad file extensively, since it is ABF-only
// and even ABF calculations generate PMF files
// Hence, we recommend that the users analyse PMF directely.
// The most common task of dealing with grad files is merging windows.
// Althrough this can be done using Colvars by running a fake ABF calculation,
// using this program should make it simple.
//
// Usage:
//   auto a = GradMerge(lowerboundary,
//                      upperboundary,
//                      width,
//                      grad_files,
//                      count_files,
//                      output_prefix)

namespace grad_merger {

    class GradMerge {
       public:
        GradMerge(const std::vector<double>& lowerboundary,
                  const std::vector<double>& upperboundary,
                  const std::vector<double>& width,
                  const std::vector<std::string>& grad_files,
                  const std::vector<std::string>& count_files,
            const std::string& output_prefix) {

            this->dimension_ = lowerboundary.size();

            assert(upperboundary.size() == this->dimension_);
            assert(width.size() == this->dimension_);
            assert(grad_files.size() == count_files.size());

            this->lowerboundary_ = lowerboundary;
            this->upperboundary_ = upperboundary;
            this->width_ = width;
            this->output_prefix_ = output_prefix;

            // the shape of internal data
            this->shape_ = std::vector<int>(this->dimension_, 0);
            for (int i = 0; i < this->dimension_; i++) {
                this->shape_[i] = int((upperboundary[i] - lowerboundary[i] +
                                       common_tools::kAccuracy) /
                                      width[i]);
            }

            count_grid_ = new ndarray::NdArray<int>({ this->shape_ });

            for (int i = 0; i < this->dimension_; i++) {
                grad_grid_.push_back(
                    ndarray::NdArray<double>(this->shape_));
            }

            for (int i = 0; i < grad_files.size(); i++) {
                this->ReadFile(grad_files[i], count_files[i]);
            }

            this->WriteData();
        }

        ~GradMerge() {
        }

       private:
        // read a grad-count file pair
        void ReadFile(const std::string& grad_file,
                      const std::string& count_file) {
            std::ifstream grad_in;
            grad_in.open(grad_file, std::ios::in);
            if (!grad_in.is_open()) {
                std::cerr << "grad file cannot be opened!" << std::endl;
                exit(1);
            }

            std::ifstream count_in;
            count_in.open(count_file, std::ios::in);
            if (!count_in.is_open()) {
                std::cerr << "count file cannot be opened!" << std::endl;
                exit(1);
            }

            // used for parsing
            std::string grad_line;
            std::string count_line;
            std::vector<std::string> grad_splited_line;
            std::vector<std::string> count_splited_line;

            std::vector<double> rc_position(this->dimension_);

            // read grad file
            while (getline(grad_in, grad_line)) {
                if (pystring::startswith(grad_line, "#")) {
                    continue;
                }
                pystring::split(grad_line, grad_splited_line);
                if (grad_splited_line.size() == 0) {
                    continue;
                }

                for (int i = 0; i < this->dimension_; i++) {
                    rc_position[i] = std::stod(grad_splited_line[i]);
                }

                // read corresponding count file
                while (getline(count_in, count_line)) {
                    if (pystring::startswith(count_line, "#")) {
                        continue;
                    }
                    pystring::split(count_line, count_splited_line);
                    if (count_splited_line.size() == 0) {
                        continue;
                    }
                    for (int i = 0; i < this->dimension_; i++) {
                        if (abs(std::stod(count_splited_line[i]) -
                                rc_position[i]) > common_tools::kAccuracy) {
                            std::cout << count_splited_line[i] << std::endl;
                            std::cout << rc_position[i] << std::endl;
                            std::cerr << "grad and count files do not "
                                         "correspond to each other!"
                                      << std::endl;
                            exit(1);
                        }
                    }
                    break;
                }

                // check if the point is in [lb, ub]
                for (int i = 0; i < this->dimension_; i++) {
                    if (rc_position[i] < this->lowerboundary_[i] ||
                        rc_position[i] > this->upperboundary_[i]) {
                        continue;
                    }
                }
                for (int i = 0; i < this->dimension_; i++) {
                    (grad_grid_[i])[this->RCToInternal(rc_position)] =
                        (std::stod(grad_splited_line[this->dimension_ + i]) *
                             std::stod(count_splited_line[this->dimension_]) +
                         (*count_grid_)[this->RCToInternal(
                             rc_position)] *
                             (grad_grid_[i])[this->RCToInternal(
                                 rc_position)]) /
                        (std::stod(count_splited_line[this->dimension_]) +
                         (*count_grid_)[this->RCToInternal(
                             rc_position)]);
                }
                (*count_grid_)[this->RCToInternal(rc_position)] +=
                    std::stod(count_splited_line[this->dimension_]);
            }

            grad_in.close();
            count_in.close();
        }

        // write internal data to files
        // in NAMD grad and count format
        // note: PBCs are not recorded! So they are zeroes!
        void WriteData() const {
            std::ofstream write_count_file, write_grad_file;
            write_count_file.open(this->output_prefix_ + ".count", std::ios::out);
            if (!write_count_file.is_open()) {
                std::cerr << "count file cannot be opened!" << std::endl;
                exit(1);
            }

            write_grad_file.open(this->output_prefix_ + ".grad", std::ios::out);
            if (!write_grad_file.is_open()) {
                std::cerr << "grad file cannot be opened!" << std::endl;
                exit(1);
            }

            // the head of NAMD files
            write_count_file << std::setw(2) << "# " << this->dimension_ << "\n";
            write_grad_file << std::setw(2) << "# " << this->dimension_ << "\n";
            for (int i = 0; i < this->dimension_; i++) {
                write_count_file << "# " << std::setw(10)
                               << this->lowerboundary_[i] << std::setw(10)
                               << this->width_[i] << std::setw(10)
                               << this->shape_[i] << " " << 0 << "\n";
                write_grad_file << "# " << std::setw(10)
                              << this->lowerboundary_[i] << std::setw(10)
                              << this->width_[i] << std::setw(10)
                              << this->shape_[i] << " " << 0 << "\n";
            }
            write_count_file << "\n";
            write_grad_file << "\n";

            // iterate over any dimension
            int n = 0;
            std::vector<int> loop_flag(this->dimension_, 0);
            while (n >= 0) {
                auto RC = this->InternalToRC(loop_flag);
                for (auto coor : RC) {
                    write_count_file << common_tools::Round(
                                          coor, common_tools::kDecimalAccuracy)
                                   << " ";
                    write_grad_file << common_tools::Round(
                                         coor, common_tools::kDecimalAccuracy)
                                  << " ";
                }
                write_count_file << (*count_grid_)[loop_flag] << "\n";
                for (int i = 0; i < this->dimension_; i++) {
                    write_grad_file << (grad_grid_[i])[loop_flag]
                                  << " ";
                }
                write_grad_file << "\n";

                // mimic an nD for loop
                n = this->dimension_ - 1;
                while (n >= 0) {
                    loop_flag[n] += 1;
                    if (loop_flag[n] > this->shape_[n] - 1) {
                        loop_flag[n] = 0;
                        n--;
                        write_count_file << "\n";
                        write_grad_file << "\n";
                    } else {
                        break;
                    }
                }
            }

            write_count_file.close();
            write_grad_file.close();
        }

        // convert external/real reaction coordinate into the internal
        // coordinate
        std::vector<int> RCToInternal(
            const std::vector<double>& rc_position) const {
            assert(rc_position.size() == this->dimension_);

            std::vector<int> internal_position(this->dimension_);
            for (int i = 0; i < this->dimension_; i++) {
                internal_position[i] =
                    int((rc_position[i] - this->lowerboundary_[i] +
                         common_tools::kAccuracy) /
                        this->width_[i]);
            }
            return internal_position;
        }

        // convert internal coordinate into external/real reaction coordinate
        std::vector<double> InternalToRC(
            const std::vector<int>& internal_position) const {
            assert(internal_position.size() == this->dimension_);

            std::vector<double> rc_position(this->dimension_);
            for (int i = 0; i < this->dimension_; i++) {
                rc_position[i] =
                    double(((internal_position[i] + 0.5) * this->width_[i]) +
                           this->lowerboundary_[i]);
            }
            return rc_position;
        }

        // the grid of grad and count
        std::vector<ndarray::NdArray<double> > grad_grid_;
        ndarray::NdArray<int> * count_grid_;
        // lowerboundary, upperboundary, width, and dimension
        std::vector<double> lowerboundary_;
        std::vector<double> upperboundary_;
        std::vector<double> width_;
        // the shape data
        std::vector<int> shape_;
        int dimension_;
        // the prefix of output file
        std::string output_prefix_;
    };

}  // namespace grad_merger

#endif  // GRADMERGER_HPP
