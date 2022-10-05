#ifndef COMMON_TOOLS_HPP
#define COMMON_TOOLS_HPP

#include <cmath>
#include <vector>

namespace common_tools {

    // the accuracy of float num calculation
    const double kAccuracy = 1e-8;
    const double kDecimalAccuracy = (-(log10(kAccuracy))) - 1;

    // boltzmann constant
    const double kBoltzmann = 0.001987191;

    inline double Round(double number, int digit) {
        return (std::round(number * pow(10, digit)) / pow(10, digit));
    }

    inline bool VectorEqual(const std::vector<double>& vec1,
                            const std::vector<double>& vec2) {
        if (vec1.size() != vec2.size()) {
            return false;
        }
        for (int i = 0; i < vec1.size(); i++) {
            if (abs(vec1[i] - vec2[i]) > kAccuracy) {
                return false;
            }
        }
        return true;
    }

    inline bool VectorInVectorOfVector(
        const std::vector<int>& vec,
        const std::vector<std::vector<int> >& vec_list) {
        for (auto& vec_item : vec_list) {
            if (vec == vec_item) {
                return true;
            }
        }
        return false;
    }

    inline bool VectorInVectorOfVector(
        const std::vector<double>& vec,
        const std::vector<std::vector<double> >& vec_list) {
        for (auto& vec_item : vec_list) {
            if (VectorEqual(vec, vec_item)) {
                return true;
            }
        }
        return false;
    }

}  // namespace common_tools

#endif  // header guard
