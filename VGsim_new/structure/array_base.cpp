#pragma once

#include <iostream>

template<class T>
class ArrayBase {
public:
    ArrayBase(uint64_t count_array, uint64_t size) 
        : count_array_(count_array)
        , size_array_(size)
        , full_size_(count_array_ * size_array_)
        , data_(new T[full_size_])
        , end_(data_ + full_size_) {}

    ~ArrayBase() {
        delete[] data_;
    }

    T& operator[](uint64_t index) {
        return data_[getIndex(index)];
    }

    const T& operator[](uint64_t index) const {
        return data_[getIndex(index)];
    }

    void debug() {
        std::cout << "Data address - " << data_ << std::endl;
        for (uint64_t array = 0; array < count_array_; ++array) {
            for (uint64_t element = 0; element < size_array_; ++element) {
                std::cout << data_[array * size_array_ + element] << " ";
            }
            std::cout << std::endl;
        }
    }

    uint64_t fastChoose(uint64_t array, T total_data_, double *random_number) {
        double rn = *random_number * total_data_;
        uint64_t index = 0;
        uint64_t offset = array * size_array_;
        T total = 0;

        while (total + data_[offset + index] < rn) {
            total += data_[offset + index];
            ++index;
        }

        *random_number = (rn - total) / data_[offset + index];
        return index;
    }

    // uint64_t fastChooseSkip(uint64_t start, T total_data_, double *random_number, uint64_t skip) {
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
        return;
    }

    void updateArray(uint64_t array) {
        return;
    }

    T getSum(uint64_t array) {
        T sum = 0;
        for (uint64_t index = array * size_array_; index < (array + 1) * size_array_; ++index) {
            sum += data_[index];
        }
        return sum;
    }

private:
    inline uint64_t getIndex(uint64_t index) const {
        return index;
    }

    uint64_t count_array_;
    uint64_t size_array_;
    uint64_t full_size_;
    T* data_;
    T* end_;
};