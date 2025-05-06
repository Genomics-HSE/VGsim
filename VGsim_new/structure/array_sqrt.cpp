#pragma once

#include <cmath>

template<class T>
class ArraySqrt {
public:
    ArraySqrt(uint64_t count_array, uint64_t size) 
        : count_array_(count_array)
        , size_array_(size)
        , count_part_array_(std::round(std::log2(size_array_)))
        , base_part_size_(size_array_ / count_part_array_)
        , double_size_array_(1 + count_part_array_ + size_array_)
        , start_array_(1 + count_part_array_)
        , full_size_(count_array_ * double_size_array_)
        , count_add_one_(size_array_ % count_part_array_)
        , r_size_(size_array_ / count_part_array_)
        , l_size_(r_size_ + (count_add_one_ > 0 ? 1 : 0))
        , data_(new T[full_size_])
        , end_(data_ + full_size_) {}

    ~ArraySqrt() {
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
        uint64_t index = 1;
        uint64_t offset = array * double_size_array_;
        T total = 0;

        while (total + data_[offset + index] < rn) {
            total += data_[offset + index];
            ++index;
        }
        index = start_array_ + getStartPart(index - 1);
        while (total + data_[offset + index] < rn) {
            total += data_[offset + index];
            ++index;
        }

        *random_number = (rn - total) / data_[offset + index];
        return index - start_array_;
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
        updateElementPrivate(array, getPartIndex(element));
        updateRoot(array);
    }

    void updateArray(uint64_t array) {
        for (uint64_t index = 0; index < count_part_array_; ++index) {
            updateElementPrivate(array, index);
        }
        updateRoot(array);
    }

    T getSum(uint64_t array) {
        return data_[array * double_size_array_];
    }

private:
    inline uint64_t getIndex(uint64_t index) const {
        return index / size_array_ * double_size_array_ + index % size_array_ + start_array_;
    }

    void updateElementPrivate(uint64_t array, uint64_t part_array) {
        uint64_t offset = array * double_size_array_;
        
        uint64_t part_index = offset + 1 + part_array;
        data_[part_index] = 0;

        uint64_t index_start = offset + start_array_ + getStartPart(part_array);
        for (uint64_t index = index_start; index < index_start + getSizePart(part_array); ++index) {
            data_[part_index] += data_[index];
        }
    }

    void updateRoot(uint64_t array) {
        uint64_t offset = array * double_size_array_;
        data_[offset] = 0;
        for (uint64_t index = offset + 1; index < offset + 1 + count_part_array_; ++index) {
            data_[offset] += data_[index];
        }
    }

    uint64_t getPartIndex(uint64_t element) {
        if (l_size_ * count_add_one_ < element) {
            return count_add_one_ + (element - l_size_ * count_add_one_) / r_size_;
        }
        return element / l_size_;
    }

    uint64_t getStartPart(uint64_t part_index) {
        if (count_add_one_ < part_index) {
            return l_size_ * count_add_one_ + r_size_ * (part_index - count_add_one_);
        }
        return l_size_ * part_index;
    }

    uint64_t getSizePart(uint64_t part_index) {
        if (count_add_one_ <= part_index) {
            return r_size_;
        }
        return l_size_;
    }

    uint64_t count_array_;
    uint64_t size_array_;
    uint64_t count_part_array_;
    uint64_t base_part_size_;
    uint64_t double_size_array_;
    uint64_t start_array_;
    uint64_t full_size_;

    uint64_t count_add_one_;
    uint64_t r_size_;
    uint64_t l_size_;

    T* data_;
    T* end_;
};