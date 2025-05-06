#pragma once

#include <cmath>

template<class T>
class ArrayTree {
public:
    ArrayTree(uint64_t count_array, uint64_t size) 
        : count_array_(count_array)
        , size_array_(size)
        , double_size_array_(std::pow(2, std::ceil(std::log2(size_array_)) + 1) - 1)
        , start_array_(std::pow(2, std::ceil(std::log2(size_array_))) - 1)
        , full_size_(count_array * double_size_array_)
        , data_(new T[full_size_])
        , end_(data_ + full_size_) {}

    ~ArrayTree() {
        delete[] data_;
    }

    T& operator[](uint64_t index) {
        return data_[getIndex(index)];
    }

    const T& operator[](uint64_t index) const {
        return data_[getIndex(index)];
    }

    void debug() const {
        for (uint64_t array = 0; array < count_array_; ++array) {
            for (uint64_t element = 0; element < double_size_array_; ++element) {
                std::cout << data_[array * double_size_array_ + element] << " ";
            }
            std::cout << std::endl;
        }
    }

    uint64_t fastChoose(uint64_t array, T total_data_, double *random_number) {
        double rn = *random_number * total_data_;
        uint64_t index = 0;
        uint64_t offset = array * double_size_array_;
        T total = 0;

        while (index < start_array_) {
            index = index * 2 + 1;
            if (total + data_[offset + index] < rn) {
                total += data_[offset + index];
                ++index;
            }
        }
        
        *random_number = (rn - total) / data_[offset + index];
        return index - start_array_;
    }

    // uint64_t fastChooseSkip(uint64_t array, T total_data_, double *random_number, uint64_t skip) {
    //     double rn = *random_number * total_data_;
    //     uint64_t index = start;
    //     if (skip == 0) {
    //         ++index;
    //     }
    //     T total = data_[index];
    //     while (total < rn) {
    //         ++index;
    //         if (index != skip) {
    //             total += data_[index];
    //         }
    //     }
    //     *random_number = (rn - (total - data_[index])) / data_[index];
    //     return index;
    // }

    void updateElement(uint64_t array, uint64_t element) {
        uint64_t offset = array * double_size_array_;
        uint64_t real_position = element + start_array_;
        bool left_right = real_position % 2 == 0; // false is left node, true is right node
        while (real_position != 0) {
            data_[offset + (real_position - 1) / 2] = data_[offset + real_position] + data_[offset + real_position + (left_right ? -1 : 1)];
            real_position = (real_position - 1) / 2;
            left_right = real_position % 2 == 0;
        }
    }

    void updateArray(uint64_t array) {
        uint64_t offset = (array + 1) * double_size_array_ - 2;
        int64_t parent_index = offset - double_size_array_ / 2;
        for (int64_t index = offset; index > (int64_t)(array * double_size_array_); index = index - 2, --parent_index) {
            data_[parent_index] = data_[index] + data_[index + 1];
        }
    }

    T getSum(uint64_t array) {
        return data_[array * double_size_array_];
    }

private:
    inline uint64_t getIndex(uint64_t index) const {
        return index / size_array_ * double_size_array_ + index % size_array_ + start_array_;
    }

    uint64_t count_array_;
    uint64_t size_array_;
    uint64_t double_size_array_;
    uint64_t start_array_;
    uint64_t full_size_;
    T* data_;
    T* end_;
};