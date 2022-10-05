#ifndef TOOL_PAIR_INTERACTION_HPP
#define TOOL_PAIR_INTERACTION_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "common_pmf_parser.hpp"
#include "third_party/pystring.hpp"

// calculate pair interaction in CV space
//
// usage:
//    // initialize
//    // one specify lb, width, ub, namd log file, colvars.traj file,
//    // columns indication CVs in colvars.traj file, whether it is a self
//    interaction calculation
//    // and the number of strides of colvars.traj file
//    auto pair =
//    pair_interaction::PairInteraction({-180,-180},{2,2},{180,180},"namd.log","namd.colvars.traj",{1,2},false,0)
//    // write result
//    pair.WriteFiles("outputPrefix")

namespace pair_interaction {

    // pair interaction in CV space
    class PairInteraction {
       public:
        PairInteraction(const std::vector<double>& lowerboundary,
                        const std::vector<double>& width,
                        const std::vector<double>& upperboundary,
                        const std::string& namd_log_file,
                        const std::string& cvtrj_file,
                        const std::vector<int>& cols,
                        bool self_interaction = false, int stride = 0) {
            ele_energy_ = new pmf_parser::Pmf<double>(lowerboundary, width,
                                                    upperboundary);
            vdw_energy_ = new pmf_parser::Pmf<double>(lowerboundary, width,
                                                    upperboundary);
            total_energy_ = new pmf_parser::Pmf<double>(lowerboundary, width,
                                                      upperboundary);
            samples_ = new pmf_parser::Pmf<double>(lowerboundary, width,
                                                  upperboundary);

            ReadFiles(cvtrj_file, cols, namd_log_file, self_interaction,
                      stride);
        }

        // write data into files
        void WriteFiles(const std::string& output_prefix) const {
            ele_energy_->WritePmfFile(output_prefix + ".eleEnergy");
            vdw_energy_->WritePmfFile(output_prefix + ".vdwEnergy");
            total_energy_->WritePmfFile(output_prefix + ".totalEnergy");
        }

        ~PairInteraction() {
            delete ele_energy_;
            delete vdw_energy_;
            delete total_energy_;
            delete samples_;
        }

       private:
        // read cvtrj-log files
        void ReadFiles(const std::string& cvtrj_file,
                       const std::vector<int>& cols,
                       const std::string& namd_log_file,
                       bool self_interaction = false, int stride = 0) {
            std::ifstream cvtrj, namd_log;

            cvtrj.open(cvtrj_file, std::ios::in);
            if (!cvtrj.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            namd_log.open(namd_log_file, std::ios::in);
            if (!namd_log.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            // used for parsing
            std::string line;
            std::string line2;
            std::vector<std::string> splited_line;
            std::vector<std::string> splited_line2;
            std::vector<double> rc_pos;
            bool read_log;
            bool available_rc = true;

            // reading cvtrj file
            // cvtrj record the 0th frame, while log file not
            bool first_loop = true;
            while (getline(cvtrj, line)) {
                rc_pos.clear();
                if (pystring::startswith(line, "#")) {
                    continue;
                }
                pystring::split(line, splited_line);
                if (splited_line.size() == 0) {
                    continue;
                }
                for (int i = 0; i < cols.size(); i++) {
                    rc_pos.push_back(std::stod(splited_line[cols[i]]) +
                                    0.5 * ele_energy_->GetWidth()[i]);
                }

                if (!first_loop) {
                    read_log = false;
                    while (getline(namd_log, line2)) {
                        if (!pystring::startswith(line2, "ENERGY:")) {
                            continue;
                        } else {
                            pystring::split(line2, splited_line2);
                            for (int i = 0; i < rc_pos.size(); i++) {
                                if (rc_pos[i] >=
                                        ele_energy_->GetUpperboundary()[i] +
                                            ele_energy_->GetWidth()[i] ||
                                    rc_pos[i] <
                                        ele_energy_->GetLowerboundary()[i]) {
                                    available_rc = false;
                                }
                            }
                            if (available_rc) {
                                (*ele_energy_)[rc_pos] +=
                                    std::stod(splited_line2[6]);
                                (*vdw_energy_)[rc_pos] +=
                                    std::stod(splited_line2[7]);
                                (*samples_)[rc_pos] += 1;
                                if (self_interaction) {
                                    (*total_energy_)[rc_pos] +=
                                        std::stod(splited_line2[11]);
                                }
                            }
                            available_rc = true;
                            read_log = true;
                            break;
                        }
                    }
                    if (!read_log) {
                        std::cout
                            << "Error, colvars.traj and log files do not match!"
                            << std::endl;
                        std::cout << "At step: " << splited_line[0] << std::endl;
                        exit(1);
                    }
                } else {
                    first_loop = false;
                }

                for (int i = 0; i < stride; i++) {
                    getline(cvtrj, line);
                    if (pystring::startswith(line, "#")) {
                        i--;
                        continue;
                    }
                }
            }
            cvtrj.close();
            namd_log.close();
            Standarize(self_interaction);
        }

        // standrdaize internal data, used by readFile
        void Standarize(bool self_interaction = false) {
            auto dimension = ele_energy_->GetDimension();
            auto shape = ele_energy_->GetShape();

            // iterate over any dimension
            int n = 0;
            std::vector<int> loop_flag(dimension, 0);
            while (n >= 0) {
                (*ele_energy_)[loop_flag] /= (*samples_)[loop_flag];
                (*vdw_energy_)[loop_flag] /= (*samples_)[loop_flag];
                if (!self_interaction) {
                    (*total_energy_)[loop_flag] =
                        (*ele_energy_)[loop_flag] + (*vdw_energy_)[loop_flag];
                } else {
                    (*total_energy_)[loop_flag] /= (*samples_)[loop_flag];
                }

                // mimic an nD for loop
                n = dimension - 1;
                while (n >= 0) {
                    loop_flag[n] += 1;
                    if (loop_flag[n] > shape[n] - 1) {
                        loop_flag[n] = 0;
                        n--;
                    } else {
                        break;
                    }
                }
            }
        }

        // store energy along CV space
        pmf_parser::Pmf<double>* ele_energy_;
        pmf_parser::Pmf<double>* vdw_energy_;
        pmf_parser::Pmf<double>* total_energy_;
        // the total number of samples in each bin
        pmf_parser::Pmf<double>* samples_;
    };

}  // namespace pair_interaction

#endif  // PAIRINTERACTION_HPP
