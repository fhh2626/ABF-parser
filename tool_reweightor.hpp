#ifndef TOOL_REIWEIGHTOR_HPP
#define TOOL_REIWEIGHTOR_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>

#include "common_pmf_parser.hpp"
#include "common_tools.hpp"
#include "third_party/pystring.hpp"

namespace reweightor {

    // reweight ABF pmf from one to another RC
    class Reweightor {
       public:
        Reweightor(const std::string& old_pmf,
                   const std::vector<std::string>& cvtrj_files,
                   const std::vector<int>& old_columns,
                   const std::vector<double>& lowerboundary,
                   const std::vector<double>& width,
                   const std::vector<double>& upperboundary,
                   const std::vector<int>& new_columns, double temperature,
                   bool is_gawtm_eabf = false,
                   const std::vector<std::string>& namd_log_files = {},
                   int ignore_samples = 0) {
            this->old_columns_ = old_columns;
            this->lowerboundary_ = lowerboundary;
            this->width_ = width;
            this->upperboundary_ = upperboundary;
            this->new_columns_ = new_columns;
            this->temperature_ = temperature;
            this->is_gawtm_eabf_ = is_gawtm_eabf;
            this->ignore_samples_ = ignore_samples;

            this->old_pmf_ = new pmf_parser::Pmf<double>(old_pmf);
            this->new_pmf_ = new pmf_parser::Pmf<double>(lowerboundary, width,
                                                         upperboundary);
            this->histogram_ = new pmf_parser::Pmf<double>(lowerboundary, width,
                                                           upperboundary);

            ReadCvtrj(cvtrj_files, namd_log_files);
            CalculatePmf();
        }

        // write data into a file
        void WriteFiles(const std::string& outputFile) const {
            histogram_->WritePmfFile(outputFile + ".histogram");
            new_pmf_->WritePmfFile(outputFile);
        }

        ~Reweightor() {
            delete old_pmf_;
            delete new_pmf_;
            delete histogram_;
        }

       private:
        void ReadCvtrj(const std::vector<std::string>& cvtrj_files,
                       const std::vector<std::string>& namd_log_files = {}) {
            // read cvtrj file and do histogram
            std::ifstream cvtrj, namd_log;
            std::string line, line2;
            std::vector<std::string> splited_line, splited_line2;
            std::vector<double> old_rc_pos, new_rc_pos;
            bool read_log;
            bool available_rc = true;
            double gamd_energy = 0;

            for (size_t i = 0; i < cvtrj_files.size(); i++) {
                cvtrj.open(cvtrj_files[i], std::ios::in);
                if (!cvtrj.is_open()) {
                    std::cerr << "cvtrj file cannot be opened!" << std::endl;
                    exit(1);
                }

                if (is_gawtm_eabf_) {
                    namd_log.open(namd_log_files[i], std::ios::in);
                    if (!namd_log.is_open()) {
                        std::cerr << "NAMDlog file cannot be opened!" << std::endl;
                        exit(1);
                    }
                }

                while (getline(cvtrj, line)) {
                    old_rc_pos.clear();
                    new_rc_pos.clear();
                    if (pystring::startswith(line, "#")) {
                        continue;
                    }
                    pystring::split(line, splited_line);
                    if (splited_line.size() == 0) {
                        continue;
                    }
                    for (int i = 0; i < old_columns_.size(); i++) {
                        old_rc_pos.push_back(
                            std::stod(splited_line[old_columns_[i]]) +
                            0.5 * old_pmf_->GetWidth()[i]);
                    }
                    for (int i = 0; i < new_columns_.size(); i++) {
                        new_rc_pos.push_back(
                            std::stod(splited_line[new_columns_[i]]) +
                            0.5 * new_pmf_->GetWidth()[i]);
                    }

                    // if reweight GaWTM-eABF
                    // write NAMD log file
                    if (is_gawtm_eabf_) {
                        read_log = false;
                        while (getline(namd_log, line2)) {
                            if (!pystring::startswith(line2, "ACCELERATED")) {
                                continue;
                            } else {
                                pystring::split(line2, splited_line2);
                                if (available_rc) {
                                    gamd_energy = std::stod(splited_line2[5]);
                                }
                                read_log = true;
                                break;
                            }
                        }
                        if (!read_log) {
                            std::cout << "Error, colvars.traj and log files do "
                                         "not match!"
                                      << std::endl;
                            std::cout << "At step: " << splited_line[0]
                                      << std::endl;
                            exit(1);
                        }
                    }

                    // check the range of rc
                    available_rc = true;
                    for (int i = 0; i < old_columns_.size(); i++) {
                        if (old_rc_pos[i] >= old_pmf_->GetUpperboundary()[i] +
                                                 old_pmf_->GetWidth()[i] ||
                            old_rc_pos[i] < old_pmf_->GetLowerboundary()[i]) {
                            available_rc = false;
                            break;
                        }
                    }
                    for (int i = 0; i < new_columns_.size(); i++) {
                        if (new_rc_pos[i] >= new_pmf_->GetUpperboundary()[i] +
                                                 new_pmf_->GetWidth()[i] ||
                            new_rc_pos[i] < new_pmf_->GetLowerboundary()[i]) {
                            available_rc = false;
                            break;
                        }
                    }

                    if (std::stod(splited_line[0]) > ignore_samples_) {
                        if (available_rc) {
                            (*histogram_)[new_rc_pos] +=
                                exp((-(*old_pmf_)[old_rc_pos] + gamd_energy) /
                                    (common_tools::kBoltzmann *
                                     this->temperature_));
                        }
                    }
                }
                cvtrj.close();

                if (is_gawtm_eabf_) {
                    namd_log.close();
                }
            }
        }

        void CalculatePmf() {
            // iterate over any dimension
            int n = 0;
            auto dimension = this->lowerboundary_.size();
            std::vector<int> loop_flag(dimension, 0);
            while (n >= 0) {
                if ((*(this->histogram_))[loop_flag] == 0) {
                    (*(this->new_pmf_))[loop_flag] =
                        -common_tools::kBoltzmann * temperature_ *
                        log(common_tools::kAccuracy);
                } else {
                    (*(this->new_pmf_))[loop_flag] =
                        -common_tools::kBoltzmann * temperature_ *
                        log((*(this->histogram_))[loop_flag]);
                }

                // mimic an nD for loop
                n = dimension - 1;
                while (n >= 0) {
                    loop_flag[n] += 1;
                    if (loop_flag[n] > new_pmf_->GetShape()[n] - 1) {
                        loop_flag[n] = 0;
                        n--;
                    } else {
                        break;
                    }
                }
            }
        }

        pmf_parser::Pmf<double>* old_pmf_;
        pmf_parser::Pmf<double>* new_pmf_;

        // record effectove histogram
        pmf_parser::Pmf<double>* histogram_;

        std::vector<int> old_columns_;
        std::vector<double> lowerboundary_;
        std::vector<double> width_;
        std::vector<double> upperboundary_;
        std::vector<int> new_columns_;

        double temperature_;
        bool is_gawtm_eabf_;
        int ignore_samples_;
    };

}  // namespace reweightor

#endif