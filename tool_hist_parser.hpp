#ifndef TOOL_HIST_PARSER_HPP
#define TOOL_HIST_PARSER_HPP

#include <cassert>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "common_pmf_parser.hpp"

// parsing hist files
//
// Usage:
//   // read NAMD hist file
//   auto hist = = hist_parser::Hist<double>("test.hist.pmf")
//   // calculate RMSD with respect to the last frame in the hist file
//   auto result = his.CalcRMSD()
//   // calculate RMSD with respect to a reference pmf
//   auto result = his.CalcRMSD(pmf_parser::Pmf<double>("test.pmf"))
//   // get data
//   his.GetData()
//   his.GetDataForFrame(frame)
//   his.GetNumberOfFrames()

namespace hist_parser {

    template <typename T>
    class Hist {
       public:
        // current only NAMD hist file is supported
        Hist(const std::string& hist_file) {
            static_assert(
                std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be a kind of number");

            std::ifstream read_file;
            // tellg() only works when using std::ios::binary
            read_file.open(hist_file, std::ios::in | std::ios::binary);
            if (!read_file.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            // used for parsing
            std::string line;
            std::vector<std::string> splited_line;
            std::streampos pos;

            while (true) {
                // position of current stream
                pos = read_file.tellg();

                // whether end of the hist file
                getline(read_file, line);
                pystring::split(line, splited_line);
                if (splited_line.size() == 0) {
                    break;
                }
                read_file.seekg(pos);

                this->history_.push_back(new pmf_parser::Pmf<T>(read_file));
            }
            read_file.close();
        }

        // get data
        const std::vector<pmf_parser::Pmf<T>*>& GetData() const {
            return this->history_;
        }

        // get data for a given frame
        pmf_parser::Pmf<T>* GetDataForFrame(int frame) const {
            return this->history_[frame];
        }

        // get number of frames
        int GetNumberOfFrames() const { return this->history_.size(); }

        // calculate the time-evolution of RMSD with respect to a reference
        std::vector<double> CalcRMSD(const pmf_parser::Pmf<T>& ref_pmf) const {
            assert(this->history_.size() != 0);
            assert(this->history_[0]->GetShape() == ref_pmf.GetShape());

            // result
            std::vector<double> result;
            for (const auto& p : history_) {
                result.push_back(ref_pmf.CalcRMSD(*p));
            }
            return result;
        }

        // calculate the time-evolution of RMSD with respect to history[-1]
        std::vector<double> CalcRMSD() const {
            assert(this->history_.size() != 0);

            auto hist_size = history_.size();
            auto& ref_pmf = *(history_[hist_size - 1]);

            // result
            std::vector<double> result;
            for (const auto& p : history_) {
                result.push_back(ref_pmf.CalcRMSD(*p));
            }
            return result;
        }

        ~Hist() {
            for (const auto& item : history_) {
                delete item;
            }
        }

       private:
        // history is actually a set of PMFs
        std::vector<pmf_parser::Pmf<T>*> history_;
    };

}  // namespace hist_parser

#endif  // HISTPARSER_HPP
