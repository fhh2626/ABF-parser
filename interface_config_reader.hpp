#ifndef INTERFACE_CONFIG_READER_HPP
#define INTERFACE_CONFIG_READER_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "common_pmf_parser.hpp"
#include "interface.hpp"
#include "third_party/ini_reader.hpp"
#include "third_party/pystring.hpp"
#include "tool_hist_parser.hpp"

namespace interface {

    enum JobType {
        kCalcError,
        kFindPathway,
        kCalcRMSD,
        kCalcFreeEnergyDiff,
        kCalcEnergyInCVSpace,
        kMergeWindows,
        kReweightPmf,
        kProjectPmf
    };

    class ConfigReader {
       public:
        ConfigReader(const std::string& ini_file, JobType job_type) {
            reader = INIReader(ini_file);
            if (reader.ParseError() != 0) {
                std::cerr << "Can't load ini file\n";
                exit(1);
            }

            // interface
            auto face = interface::Interface();

            // temp variables that may be used in parsing the config file
            std::string hist_path, pmf_path, old_pmf, final_columns;
            std::vector<double> align_position;
            std::vector<double> lowerboundary, width, upperboundary;
            std::vector<double> initial_point, end_point;
            std::vector<bool> pbc;
            std::string output_prefix, output_path, namd_log_file, cvtrj_file;
            std::vector<std::vector<double> > targeted_points;
            std::vector<std::vector<double> > force_constant;
            std::vector<double> point1, point2;
            std::vector<int> cols, old_columns, new_columns;
            std::vector<std::string> grad_files, count_files, cvtrj_files,
                namd_log_files;
            double cutoff_energy1, cutoff_energy2, temperature;
            bool write_explored_points, self_interaction, is_gawtm_eabf;
            int stride, ignore_samples;

            switch (job_type) {
                case kCalcError:
                    ConfigCalcError(hist_path, align_position, output_path);
                    face.CalcError(hist_parser::Hist<double>(hist_path),
                                   align_position, output_path);
                    std::cout << "error calculation finished!" << std::endl;
                    break;
                case kFindPathway:
                    ConfigFindPathway(pmf_path, lowerboundary, width,
                                      upperboundary, initial_point, end_point,
                                      pbc, output_prefix, targeted_points,
                                      force_constant, write_explored_points);
                    face.FindPathway(pmf_parser::Pmf<double>(pmf_path),
                                     initial_point, end_point, pbc,
                                     output_prefix, targeted_points,
                                     force_constant, write_explored_points);
                    std::cout << "pathway calculation finished!" << std::endl;
                    break;
                case kCalcRMSD:
                    ConfigCalcRMSD(hist_path, pmf_path, output_path);
                    if (pmf_path == "") {
                        face.CalcRMSD(hist_parser::Hist<double>(hist_path),
                                      output_path);
                    } else {
                        face.CalcRMSD(hist_parser::Hist<double>(hist_path),
                                      pmf_parser::Pmf<double>(pmf_path),
                                      output_path);
                    }
                    std::cout << "RMSD calculation finished!" << std::endl;
                    break;
                case kCalcFreeEnergyDiff:
                    ConfigCalcFreeEnergyDiff(pmf_path, point1, cutoff_energy1,
                                             point2, cutoff_energy2,
                                             temperature);
                    free_energy_difference = face.CalcFreeEnergyDiff(
                        pmf_parser::Pmf<double>(pmf_path), point1,
                        cutoff_energy1, point2, cutoff_energy2, temperature);
                    std::cout << "free-energy difference is "
                              << free_energy_difference << std::endl;
                    std::cout << "free-energy difference calculation finished!"
                              << std::endl;
                    break;
                case kCalcEnergyInCVSpace:
                    ConfigCalcEnergyInCVSpace(lowerboundary, width,
                                              upperboundary, namd_log_file,
                                              cvtrj_file, cols, output_prefix,
                                              self_interaction, stride);
                    face.CalcEnergyInCVSpace(lowerboundary, width,
                                             upperboundary, namd_log_file,
                                             cvtrj_file, cols, output_prefix,
                                             self_interaction, stride);
                    std::cout
                        << "pair interaction in CV space calculation finished"
                        << std::endl;
                    break;
                case kMergeWindows:
                    ConfigMergeWindows(lowerboundary, upperboundary, width,
                                       grad_files, count_files, output_prefix);
                    face.MergeWindows(lowerboundary, upperboundary, width,
                                      grad_files, count_files, output_prefix);
                    std::cout << "merging windows finished" << std::endl;
                    break;
                case kReweightPmf:
                    ConfigReweightPmf(
                        old_pmf, cvtrj_files, old_columns, lowerboundary, width,
                        upperboundary, new_columns, temperature, is_gawtm_eabf,
                        namd_log_files, ignore_samples, output_path);
                    face.ReweightPmf(
                        old_pmf, cvtrj_files, old_columns, lowerboundary, width,
                        upperboundary, new_columns, temperature, is_gawtm_eabf,
                        namd_log_files, ignore_samples, output_path);
                    std::cout << "reweighting finished" << std::endl;
                    break;
                case kProjectPmf:
                    ConfigProjectPmf(pmf_path, cols, temperature, output_path);
                    face.ProjectPmf(pmf_parser::Pmf<double>(pmf_path), cols,
                                    temperature, output_path);
                    std::cout << "projecting finished" << std::endl;
                    break;

                default:
                    std::cerr << "Error, job type not supported!" << std::endl;
                    exit(1);
                    break;
            }
        }

        ~ConfigReader() {}

        inline double GetFreeEnergyDifference() const {
            return free_energy_difference;
        }

       private:
        // parse config file of CalcError
        void ConfigCalcError(std::string& hist_path,
                             std::vector<double>& align_position,
                             std::string& output_file) {
            hist_path = reader.Get("CalcError", "historyDirectory", "");

            auto align_position_string =
                reader.Get("CalcError", "alignPosition", "");
            align_position = StringToVector(align_position_string, 0.0);

            std::vector<std::string> temp_output_file;
            pystring::rpartition(hist_path, ".", temp_output_file);
            output_file = temp_output_file[0] + ".error";
        }

        // parse config file of FindPathway
        void ConfigFindPathway(
            std::string& pmf_path, std::vector<double>& lowerboundary,
            std::vector<double>& width, std::vector<double>& upperboundary,
            std::vector<double>& initial_point, std::vector<double>& end_point,
            std::vector<bool>& pbc, std::string& output_prefix,
            std::vector<std::vector<double> >& targeted_points,
            std::vector<std::vector<double> >& force_constant,
            bool& write_explored_points) {
            pmf_path = reader.Get("FindPathway", "directory", "");
            // these vars will be converted to std::vector
            auto lowerboundary_string =
                reader.Get("FindPathway", "lowerboundary", "");
            auto upperboundary_string =
                reader.Get("FindPathway", "upperboundary", "");
            auto width_string = reader.Get("FindPathway", "width", "");
            auto initial_string = reader.Get("FindPathway", "initial", "");
            auto end_string = reader.Get("FindPathway", "end", "");
            auto pbc_string = reader.Get("FindPathway", "pbc", "");
            auto targeted_points_and_fc_string =
                reader.Get("FindPathway", "target", "");

            write_explored_points =
                reader.GetBoolean("FindPathway", "writeExploredPoints", false);

            if (lowerboundary_string != "" && upperboundary_string != "" &&
                width_string != "") {
                lowerboundary = StringToVector(lowerboundary_string, 0.0);
                upperboundary = StringToVector(upperboundary_string, 0.0);
                width = StringToVector(width_string, 0.0);
            }

            initial_point = StringToVector(initial_string, 0.0);
            end_point = StringToVector(end_string, 0.0);
            pbc = StringToVector(end_string, false);

            // targeted points
            if (targeted_points_and_fc_string != "") {
                std::vector<double> point, fc;
                std::vector<std::string> splited_targeted_points_and_fc;
                pystring::split(targeted_points_and_fc_string,
                                splited_targeted_points_and_fc, ",");
                for (int i = 0; i < splited_targeted_points_and_fc.size() /
                                        initial_point.size() / 2;
                     i++) {
                    point = {};
                    fc = {};
                    for (int j = 0; j < initial_point.size(); j++) {
                        point.push_back(
                            std::stod(splited_targeted_points_and_fc
                                          [initial_point.size() * 2 * i + j]));
                        fc.push_back(
                            std::stod(splited_targeted_points_and_fc
                                          [initial_point.size() * 2 * i +
                                           initial_point.size() + j]));
                    }
                    targeted_points.push_back(point);
                    force_constant.push_back(fc);
                }
            }

            std::vector<std::string> temp_output_prefix;
            pystring::rpartition(pmf_path, ".", temp_output_prefix);
            output_prefix = temp_output_prefix[0];
        }

        // parse config file of CalcRMSD
        void ConfigCalcRMSD(std::string& hist_path, std::string& pmf_path,
                            std::string& output_path) {
            hist_path = reader.Get("CalcRMSD", "historyDirectory", "");
            pmf_path = reader.Get("CalcRMSD", "pmfDirectory", "");
            std::vector<std::string> temp_output_file;
            pystring::rpartition(hist_path, ".", temp_output_file);
            output_path = temp_output_file[0] + ".rmsd";
        }

        // parse config file of CalcFreeEnergyDiff
        void ConfigCalcFreeEnergyDiff(std::string& pmf_path,
                                      std::vector<double>& point1,
                                      double& cutoff_energy1,
                                      std::vector<double>& point2,
                                      double& cutoff_energy2,
                                      double& temperature) {
            pmf_path = reader.Get("CalcFreeEnergyDiff", "pmfDirectory", "");
            auto point1_string = reader.Get("CalcFreeEnergyDiff", "point1", "");
            auto point2_string = reader.Get("CalcFreeEnergyDiff", "point2", "");
            point1 = StringToVector(point1_string, 0.0);
            point2 = StringToVector(point2_string, 0.0);
            cutoff_energy1 =
                reader.GetReal("CalcFreeEnergyDiff", "cutoffEnergy1", 3.0);
            cutoff_energy2 =
                reader.GetReal("CalcFreeEnergyDiff", "cutoffEnergy2", 3.0);
            temperature =
                reader.GetReal("CalcFreeEnergyDiff", "temperature", 300);
        }

        // parse config file of CalcEnergyInCVSpace
        // or pair interaction in CV space
        void ConfigCalcEnergyInCVSpace(
            std::vector<double>& lowerboundary, std::vector<double>& width,
            std::vector<double>& upperboundary, std::string& namd_log_file,
            std::string& cvtrj_file, std::vector<int>& cols,
            std::string& output_prefix, bool& self_interaction, int& stride) {
            auto lowerboundary_string =
                reader.Get("CalcEnergyInCVSpace", "lowerboundary", "");
            auto upperboundary_string =
                reader.Get("CalcEnergyInCVSpace", "upperboundary", "");
            auto width_string = reader.Get("CalcEnergyInCVSpace", "width", "");
            auto cols_string = reader.Get("CalcEnergyInCVSpace", "columns", "");
            lowerboundary = StringToVector(lowerboundary_string, 0.0);
            upperboundary = StringToVector(upperboundary_string, 0.0);
            width = StringToVector(width_string, 0.0);
            cols = StringToVector(cols_string, 0);

            namd_log_file =
                reader.Get("CalcEnergyInCVSpace", "namdLogDirectory", "");
            cvtrj_file =
                reader.Get("CalcEnergyInCVSpace", "namdCvtrjDirectory", "");

            std::vector<std::string> temp_output_prefix;
            pystring::rpartition(namd_log_file, ".", temp_output_prefix);
            output_prefix = temp_output_prefix[0];

            self_interaction = reader.GetBoolean("CalcEnergyInCVSpace",
                                                 "selfInteraction", false);
            stride = reader.GetInteger("CalcEnergyInCVSpace", "stride", 0);
        }

        // parse config file of MergeWindows
        void ConfigMergeWindows(std::vector<double>& lowerboundary,
                                std::vector<double>& upperboundary,
                                std::vector<double>& width,
                                std::vector<std::string>& grad_files,
                                std::vector<std::string>& count_files,
                                std::string& output_prefix) {
            auto lowerboundary_string =
                reader.Get("MergeWindows", "lowerboundary", "");
            auto upperboundary_string =
                reader.Get("MergeWindows", "upperboundary", "");
            auto width_string = reader.Get("MergeWindows", "width", "");
            lowerboundary = StringToVector(lowerboundary_string, 0.0);
            upperboundary = StringToVector(upperboundary_string, 0.0);
            width = StringToVector(width_string, 0.0);

            auto grad_files_string =
                reader.Get("MergeWindows", "gradFiles", "");
            auto count_files_string =
                reader.Get("MergeWindows", "countFiles", "");
            grad_files = StringToVector(grad_files_string, std::string());
            count_files = StringToVector(count_files_string, std::string());

            std::vector<std::string> temp_output_prefix;
            pystring::rpartition(grad_files[0], ".", temp_output_prefix);
            output_prefix = temp_output_prefix[0] + ".merged";
        }

        // parse config file of ReweightPmf
        void ConfigReweightPmf(
            std::string& old_pmf, std::vector<std::string>& cvtrj_files,
            std::vector<int>& old_columns, std::vector<double>& lowerboundary,
            std::vector<double>& width, std::vector<double>& upperboundary,
            std::vector<int>& new_columns, double& temperature,
            bool& is_gawtm_eabf, std::vector<std::string>& namd_log_files,
            int& ignore_samples, std::string& output_path) {
            old_pmf = reader.Get("ReweightPmf", "oldPmf", "");
            auto cvtrj_files_string =
                reader.Get("ReweightPmf", "cvtrjFiles", "");
            cvtrj_files = StringToVector(cvtrj_files_string, std::string());
            auto lowerboundary_string =
                reader.Get("ReweightPmf", "lowerboundary", "");
            auto upperboundary_string =
                reader.Get("ReweightPmf", "upperboundary", "");
            auto width_string = reader.Get("ReweightPmf", "width", "");
            lowerboundary = StringToVector(lowerboundary_string, 0.0);
            upperboundary = StringToVector(upperboundary_string, 0.0);
            width = StringToVector(width_string, 0.0);

            auto old_columns_string =
                reader.Get("ReweightPmf", "oldColumns", "");
            auto new_columns_string =
                reader.Get("ReweightPmf", "newColumns", "");
            old_columns = StringToVector(old_columns_string, 0);
            new_columns = StringToVector(new_columns_string, 0);

            temperature = reader.GetReal("ReweightPmf", "temperature", 300);
            ignore_samples =
                reader.GetInteger("ReweightPmf", "ignoreSamples", 0);

            is_gawtm_eabf = reader.GetBoolean("ReweightPmf", "gawtm", false);
            auto namd_log_files_string =
                reader.Get("ReweightPmf", "namdLogFiles", "");
            namd_log_files =
                StringToVector(namd_log_files_string, std::string());

            std::vector<std::string> temp_output_path;
            pystring::rpartition(cvtrj_files[0], ".", temp_output_path);
            output_path = temp_output_path[0] + ".reweight";
        }

        // parse config file of ProjectPmf
        void ConfigProjectPmf(std::string& pmf_path,
                              std::vector<int>& final_columns,
                              double& temperature, std::string& output_path) {
            pmf_path = reader.Get("ProjectPmf", "Pmf", "");
            auto final_columns_string =
                reader.Get("ProjectPmf", "finalColumns", "");
            final_columns = StringToVector(final_columns_string, 0);
            temperature = reader.GetReal("ProjectPmf", "temperature", 300);
            std::vector<std::string> temp_output_path;
            pystring::rpartition(pmf_path, ".", temp_output_path);
            output_path = temp_output_path[0] + ".project";
        }

        // convert a string splited by ',' to a vector
        template <typename T>
        std::vector<T> StringToVector(std::string& input_string,
                                      T type_indicator) {
            std::vector<T> result;
            std::vector<std::string> splited_string;
            std::stringstream stream;
            T temp_variable;
            if (input_string != "") {
                pystring::split(input_string, splited_string, ",");
                for (auto& item : splited_string) {
                    stream << item;
                    stream >> temp_variable;
                    result.push_back(temp_variable);

                    // clear stream
                    std::stringstream().swap(stream);
                }
            }
            return result;
        }

        INIReader reader;

        // free-energy difference as a result
        double free_energy_difference;
    };
}  // namespace interface

#endif