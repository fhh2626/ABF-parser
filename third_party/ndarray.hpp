#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <cassert>
// #define NDEBUG

#include <iostream>
#include <vector>

// a lightweighted n-dimensional array library
// Haohao Fu (fhh2626@gmail.com)
// version 0.12 beta
//
// usage:
//   // initialize by shape (vector<int>) and default value (int, default 0)
//   ndarray<int> arr({5,4}, 1);
//   // copy constructor
//   ndarray<double> arr2 = arr;
//   // calculation
//   std::cout << arr2 * 5 + 1;
//   // reshape
//   // one must understand the mechanism of the reshape function
//   // internal data are stored in a 1d array
//   // ndarray.Reshape() simply change the arrangement of the 1d array, for
//   instance
//   // [1,2,3,4,5,6] (shape{6}) -> [[1,2,3],[4,5,6]] (shape{2,3}) ->
//   [[1,2],[3,4],[5,6]] (shape{3,2}) 
//   arr.Reshape({4,5});
//   // return the max/min value
//   arr.MaxValue()
//   // get other information of the array
//   arr.GetShape()
//   arr.GetTotalSize()
//   arr.GetCArray()
//

namespace ndarray {

    // NdArray class
    template <typename T>
    class NdArray {
       public:
        // default constructor of NdArray
        NdArray(const std::vector<int>& shape, T default_value = 0) {
            static_assert(
                std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be a kind of number");

            // the shape of the tensor
            this->shape_ = shape;

            // calculate the length of the internal array
            this->total_size_ = 1;
            for (auto i : shape) {
                this->total_size_ *= i;
            }

            assert(this->total_size_ != 0);

            // get memory
            this->data_ = new T[total_size_];
            // initialize the new array
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] = default_value;
            }
        }

        // copy constructor
        NdArray(const NdArray& arr) {
            this->total_size_ = arr.GetTotalSize();
            this->shape_ = arr.GetShape();
            this->data_ = new T[this->total_size_];
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] = arr.GetCArray()[i];
            }
        }

        // type change constructor
        template <typename U>
        NdArray(const NdArray<U>& arr) {
            this->total_size_ = arr.GetTotalSize();
            this->shape_ = arr.GetShape();
            this->data_ = new T[this->total_size_];
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] = T(arr.GetCArray()[i]);
            }
        }

        // operator[] to get the desired item at a given pos
        T& operator[](const std::vector<int>& pos) {
            return data_[Position(pos)];
        }

        // operator[] to get the desired item at a given pos
        const T& operator[](const std::vector<int>& pos) const {
            return data_[Position(pos)];
        }

        // operator +
        NdArray operator+(const NdArray& arr) const {
            assert(this->shape_ == arr.shape_);
            NdArray sum_arr = arr;
            for (int i = 0; i < arr.total_size_; i++) {
                sum_arr.data_[i] += this->data_[i];
            }
            return sum_arr;
        }

        NdArray operator+(const T& num) const {
            NdArray sum_arr = *this;
            for (int i = 0; i < this->total_size_; i++) {
                sum_arr.data_[i] += num;
            }
            return sum_arr;
        }

        friend NdArray operator+(const T& num, const NdArray& arr) {
            return arr + num;
        }

        // operator +=
        NdArray operator+=(const NdArray& arr) {
            assert(this->shape_ == arr.shape_);
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] += arr.data_[i];
            }
            return *this;
        }

        NdArray& operator+=(const T& num) {
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] += num;
            }
            return *this;
        }

        // operator -
        NdArray operator-(const NdArray& arr) const {
            assert(this->shape_ == arr.shape_);
            NdArray sum_arr = *this;
            for (int i = 0; i < arr.total_size_; i++) {
                sum_arr.data_[i] -= arr.data_[i];
            }
            return sum_arr;
        }

        NdArray operator-(const T& num) const {
            NdArray sum_arr = *this;
            for (int i = 0; i < this->total_size_; i++) {
                sum_arr.data_[i] -= num;
            }
            return sum_arr;
        }

        // operator -=
        NdArray operator-=(const NdArray& arr) {
            assert(this->shape_ == arr.shape_);
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] -= arr.data_[i];
            }
            return *this;
        }

        NdArray& operator-=(const T& num) {
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] -= num;
            }
            return *this;
        }

        // operator *
        NdArray operator*(const NdArray& arr) const {
            assert(this->shape_ == arr.shape_);
            NdArray sum_arr = *this;
            for (int i = 0; i < arr.total_size_; i++) {
                sum_arr.data_[i] *= arr.data_[i];
            }
            return sum_arr;
        }

        NdArray operator*(const T& num) const {
            NdArray sum_arr = *this;
            for (int i = 0; i < this->total_size_; i++) {
                sum_arr.data_[i] *= num;
            }
            return sum_arr;
        }

        friend NdArray operator*(const T& num, const NdArray& arr) {
            return arr * num;
        }

        // operator *=
        NdArray operator*=(const NdArray& arr) {
            assert(this->shape_ == arr.shape_);
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] *= arr.data_[i];
            }
            return *this;
        }

        NdArray& operator*=(const T& num) {
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] *= num;
            }
            return *this;
        }

        // operator /
        NdArray operator/(const NdArray& arr) const {
            assert(this->shape_ == arr.shape_);
            NdArray sum_arr = *this;
            for (int i = 0; i < arr.total_size_; i++) {
                sum_arr.data_[i] /= arr.data_[i];
            }
            return sum_arr;
        }

        NdArray operator/(const T& num) const {
            NdArray sum_arr = *this;
            for (int i = 0; i < this->total_size_; i++) {
                sum_arr.data_[i] /= num;
            }
            return sum_arr;
        }

        // operator /=
        NdArray operator/=(const NdArray& arr) {
            assert(this->shape_ == arr.shape_);
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] /= arr.data_[i];
            }
            return *this;
        }

        NdArray& operator/=(const T& num) {
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] /= num;
            }
            return *this;
        }

        // operator %
        NdArray operator%(const NdArray& arr) const {
            assert(this->shape_ == arr.shape_);
            NdArray sum_arr = *this;
            for (int i = 0; i < arr.total_size_; i++) {
                sum_arr.data_[i] %= arr.data_[i];
            }
            return sum_arr;
        }

        NdArray operator%(const T& num) const {
            NdArray sum_arr = *this;
            for (int i = 0; i < this->total_size_; i++) {
                sum_arr.data_[i] %= num;
            }
            return sum_arr;
        }

        // operator %=
        NdArray operator%=(const NdArray& arr) {
            assert(this->shape_ == arr.shape_);
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] %= arr.data_[i];
            }
            return *this;
        }

        NdArray& operator%=(const T& num) {
            for (int i = 0; i < this->total_size_; i++) {
                this->data_[i] %= num;
            }
            return *this;
        }

        // operator <<
        friend std::ostream& operator<<(std::ostream& out, const NdArray& arr) {
            // iterate over any dimension
            int n = 0;
            std::vector<int> loop_flag(arr.shape_.size(), 0);
            while (n >= 0) {
                out << arr.data_[arr.Position(loop_flag)] << " ";

                // mimic an nD for loop
                n = arr.shape_.size() - 1;
                while (n >= 0) {
                    loop_flag[n] += 1;
                    if (loop_flag[n] > arr.shape_[n] - 1) {
                        loop_flag[n] = 0;
                        n--;
                        out << "\n";
                    } else {
                        break;
                    }
                }
            }
            return out;
        }

        // return the max value
        T MaxValue() const {
            T max_n = this->data_[0];
            for (int i = 0; i < this->total_size_; i++) {
                if (this->data_[i] > max_n) {
                    max_n = data_[i];
                }
            }
            return max_n;
        }

        // return the min value
        // return the max value
        T MinValue() const {
            T min_n = this->data_[0];
            for (int i = 0; i < this->total_size_; i++) {
                if (this->data_[i] < min_n) {
                    min_n = data_[i];
                }
            }
            return min_n;
        }

        // reshape the nd array
        void Reshape(const std::vector<int>& new_shape) {
            // calculate the length corresponding the new shape
            int total_size = 1;
            for (auto i : new_shape) {
                total_size *= i;
            }

            assert(this->total_size_ = total_size);

            this->shape_ = new_shape;
        }

        // return the totol size of the ndArray
        int GetTotalSize() const { return this->total_size_; }

        // return the shape of the ndArray
        const std::vector<int>& GetShape() const { return this->shape_; }

        // get the C-style Array
        const T* const GetCArray() const { return this->data_; }

        // default destructor
        ~NdArray() { delete[] data_; }

       private:
        // find the real position of an item based on a vector
        int Position(const std::vector<int>& pos) const {
            assert(pos.size() == this->shape_.size());

            int real_pos = 0;
            int temp_sum = 0;
            for (int i = 0; i < pos.size(); i++) {
                temp_sum = pos[i];
                for (int j = i + 1; j < pos.size(); j++) {
                    temp_sum *= this->shape_[j];
                }
                real_pos += temp_sum;
            }
            return real_pos;
        }

        // pointers to data
        T* data_;
        // the total number of items
        int total_size_;
        // the shape of the NdArray
        std::vector<int> shape_;
    };
}  // namespace NdArray

#endif  // NDARRAY_HPP
