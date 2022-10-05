#ifndef TOOL_PROJECTOR_HPP
#define TOOL_PROJECTOR_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "common_pmf_parser.hpp"

namespace projector {

    class Projector {
       public:
        Projector(const pmf_parser::Pmf<double>& pmf_data,
                  const std::vector<int>& final_axis, double temperature) {
            pmf_data_ = new pmf_parser::Pmf<double>(pmf_data);
            final_axis_ = final_axis;
            temperature_ = temperature;

            // the variable for PMF projected on
            std::vector<double> lowerboundary, upperboundary, width;
            for (size_t i = 0; i < final_axis_.size(); i++) {
                lowerboundary.push_back(
                    pmf_data_->GetLowerboundary()[final_axis_[i]]);
                width.push_back(pmf_data_->GetWidth()[final_axis_[i]]);
                upperboundary.push_back(
                    pmf_data_->GetUpperboundary()[final_axis_[i]]);
            }
            projected_pmf_ =
                new pmf_parser::Pmf<double>(lowerboundary, width, upperboundary);
            
            projection();
        }

        void WriteFile(const std::string& output_file) {
            projected_pmf_->WritePmfFile(output_file);
        }

        ~Projector() {
            delete projected_pmf_;
            delete pmf_data_;
        }

       private:
        void projection() {
            // axis to be averaged
            std::vector<int> other_axis;
            for (int i = 0; i < pmf_data_->GetDimension(); i++) {
                auto iter =
                    std::find(final_axis_.begin(), final_axis_.end(), i);
                if (iter == final_axis_.end()) {
                    other_axis.push_back(i);
                }
            }

            // shape of each axis
            std::vector<int> final_axis_shape, other_axis_shape;
            for (size_t i = 0; i < pmf_data_->GetDimension(); i++) {
                auto iter =
                    std::find(final_axis_.begin(), final_axis_.end(), i);
                if (iter == final_axis_.end()) {
                    other_axis_shape.push_back(pmf_data_->GetShape()[i]);
                } else {
                    final_axis_shape.push_back(pmf_data_->GetShape()[i]);
                }
            }

            // iterate over any dimension
            std::vector<int> loop_flag_outer(final_axis_shape.size(), 0);
            std::vector<int> loop_flag_inner(other_axis_shape.size(), 0);
            std::vector<int> current_position(pmf_data_->GetDimension(), 0);
            int n = loop_flag_outer.size();
            int m = loop_flag_inner.size();
            double numerator, denominator;
            while (n >= 0) {
                for (size_t i = 0, j = 0; i < pmf_data_->GetDimension(); i++) {
                    auto iter =
                        std::find(final_axis_.begin(), final_axis_.end(), i);
                    if (iter == final_axis_.end()) {
                        current_position[i] = loop_flag_outer[j];
                        j++;
                    }
                }

                numerator = 0;
                denominator = 0;
                // inner loop
                m = loop_flag_inner.size();
                while (m >= 0) {
                    // got the current position
                    for (size_t i = 0, j = 0; i < pmf_data_->GetDimension();
                         i++) {
                        auto iter = std::find(final_axis_.begin(),
                                              final_axis_.end(), i);
                        if (iter != final_axis_.end()) {
                            current_position[i] = loop_flag_inner[j];
                            j++;
                        }
                    }
                    numerator += (*pmf_data_)[current_position] *
                                 exp(-(*pmf_data_)[current_position] /
                                     (common_tools::kBoltzmann * temperature_));
                    denominator +=
                        exp(-(*pmf_data_)[current_position] /
                            (common_tools::kBoltzmann * temperature_));

                    // mimic an mD for loop
                    m = other_axis_shape.size() - 1;
                    while (m >= 0) {
                        loop_flag_inner[m] += 1;
                        if (loop_flag_inner[m] > other_axis_shape[m] - 1) {
                            loop_flag_inner[m] = 0;
                            m--;
                        } else {
                            break;
                        }
                    }
                }

                //std::cout << numerator << std::endl;
                (*projected_pmf_)[loop_flag_outer] = numerator / denominator;

                // mimic an nD for loop
                n = final_axis_shape.size() - 1;
                while (n >= 0) {
                    loop_flag_outer[n] += 1;
                    if (loop_flag_outer[n] > final_axis_shape[n] - 1) {
                        loop_flag_outer[n] = 0;
                        n--;
                    } else {
                        break;
                    }
                }
            }
        }

        // old pmf
        const pmf_parser::Pmf<double>* pmf_data_;

        // new pmf
        const pmf_parser::Pmf<double>* projected_pmf_;

        // the axis to be projected on
        std::vector<int> final_axis_;

        double temperature_;
    };

}  // namespace projector

#endif  // TOOL_PROJECTOR_HPP
