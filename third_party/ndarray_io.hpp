#ifndef NDARRAYIO_HPP
#define NDARRAYIO_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ndarray.hpp"
#include "pystring.hpp"

// usage:
//   // read a dat file into an ndarray
//   // the second parameter determines the data type
//   auto a = ndarray::ReadDat(file, 0.0);
//   // write an ndarray into a external file
//   ndarray::WriteDat("file.txt", arr);
//
// note:
//   if this file is included, one must also include pystring
//

namespace ndarray {

    // read a datfile into a 2d NdArray
    template <typename T>
    inline NdArray<T> ReadDat(const std::string& file, T dummyVar) {
        static_assert(
            std::is_integral<T>::value || std::is_floating_point<T>::value,
            "T must be a number");

        std::ifstream read_file;
        read_file.open(file, std::ios::in);
        if (!read_file.is_open()) {
            std::cerr << "file cannot open!" << std::endl;
            exit(1);
        }

        // all lines in memory
        std::vector<std::vector<std::string> > all_lines;

        // used for parsing
        std::string line;
        std::vector<std::string> splited_line;

        // read the file into memory
        while (getline(read_file, line)) {
            if (pystring::startswith(line, "#")) {
                continue;
            }
            pystring::split(line, splited_line);
            if (splited_line.size() == 0) {
                continue;
            }
            all_lines.push_back(splited_line);
        }

        // the shape of the matrix
        // the col num is determined by the first row
        std::vector<int> shape = {int(all_lines.size()),
                                  int(all_lines[0].size())};
        NdArray<T> data(shape);

        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                if (j < all_lines[i].size()) {
                    data[{i, j}] = std::stod(all_lines[i][j]);
                } else {
                    data[{i, j}] = 0;
                }
            }
        }

        read_file.close();

        return data;
    }

    // write an 2d NdArray to a datFile
    template <typename T>
    inline void WriteDat(const std::string& file, NdArray<T> arr) {
        std::ofstream write_file;
        write_file.open(file, std::ios::out);
        if (!write_file.is_open()) {
            std::cerr << "file cannot open!" << std::endl;
            exit(1);
        }

        auto shape = arr.GetShape();
        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                write_file << arr[{i, j}] << " ";
            }
            write_file << "\n";
        }

        write_file.close();
    }
}  // namespace ndarray

#endif  // NDARRAYIO_HPP
