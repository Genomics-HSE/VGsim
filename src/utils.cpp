#pragma once

#include "utils.h"

template<class T>
uint64_t fastChoose(T* weight, T total_weight, double *random_number) {
    double rn = *random_number * total_weight;
    uint64_t index = 0;
    T total = weight[index];
    while (total < rn) {
        ++index;
        total += weight[index];
    }
    *random_number = (rn - (total - weight[index])) / weight[index];
    return index;
}

template<class T>
uint64_t fastChooseSkip(T* weight, T total_weight, double *random_number, uint64_t skip) {
    double rn = *random_number * total_weight;
    uint64_t index = 0;
    if (skip == 0) {
        ++index;
    }
    T total = weight[index];
    while (total < rn) {
        ++index;
        if (index != skip) {
            total += weight[index];
        }
    }
    *random_number = (rn - (total - weight[index])) / weight[index];
    return index;
}

template<class T>
void PrintArray1nd(std::string text, T* array, uint64_t size) {
    uint64_t index = 0;
    std::cout << text;
    for (uint64_t _ = 0; _ < size; ++_) {
        std::cout << " " << array[index++];
    }
    std::cout << std::endl;
}

template<class T>
void PrintArray2nd(std::string text, T* array, uint64_t size1, uint64_t size2) {
    uint64_t index = 0;
    std::cout << text << std::endl;
    for (uint64_t _ = 0; _ < size1; ++_) {
        for (uint64_t _ = 0; _ < size2; ++_) {
            std::cout << array[index++] << " ";
        }
        std::cout << std::endl;
    }
}

template<class T>
void PrintArray3nd(std::string text, T* array, uint64_t size1, uint64_t size2, uint64_t size3) {
    uint64_t index = 0;
    std::cout << text << std::endl;
    for (uint64_t _ = 0; _ < size1; ++_) {
        for (uint64_t _ = 0; _ < size2; ++_) {
            for (uint64_t _ = 0; _ < size3; ++_) {
                std::cout << array[index++] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

template<class T>
void PrintArray4nd(std::string text, T* array, uint64_t size1, uint64_t size2, uint64_t size3, uint64_t size4) {
    uint64_t index = 0;
    std::cout << text << std::endl;
    for (uint64_t _ = 0; _ < size1; ++_) {
        for (uint64_t _ = 0; _ < size2; ++_) {
            for (uint64_t _ = 0; _ < size3; ++_) {
                for (uint64_t _ = 0; _ < size4; ++_) {
                    std::cout << array[index++] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void PrintMutation(const Mutation& mutation) {
    std::cout << "Node: " << mutation.node
              << ", site: " << mutation.site
              << ", derived: " << mutation.derivedState
              << ", ancestral: " << mutation.ancestralState
              << ", time: " << mutation.time
              << "\n";
}

void PrintMigration(uint64_t node, const Migration& migration) {
    std::cout << "Node: " << node
              << ", old: " << migration.oldPopulation
              << ", new: " << migration.newPopulation
              << ", time: " << migration.time
              << "\n";
}

uint64_t GetNewHaplotype(uint64_t haplotype, uint64_t site, uint64_t DS, uint64_t sites) {
    uint64_t digit4 = 0;
    uint64_t AS = 0;

    digit4 = std::pow(4, sites - site - 1);
    AS = haplotype / digit4 % 4;
    if (DS >= AS) {
        ++DS; 
    }
    return haplotype + (DS - AS) * digit4;
}

std::vector<uint64_t> GetIndexes(int64_t index, uint64_t max_index) {
    std::vector<uint64_t> indexes;
    if (index == -1) {
        indexes.reserve(max_index);
        for (uint64_t i = 0; i < max_index; ++i) {
            indexes.push_back(i);
        }
    } else {
        indexes.push_back(static_cast<uint64_t>(index));
    }
    return indexes;
}

void CheckValue(double rate, std::string smth) {
    if (rate < 0) {
        throw std::runtime_error("Incorrect value of " + smth + ". Value should be more or equal 0.");
    }
}

// def check_value(self, value, smth, edge=None, none=False):
//     if none:
//         if isinstance(value, (int, float)) == False and value is not None:
//             raise TypeError('Incorrect type of ' + smth + '. Type should be int or float or None.')
//     else:
//         if isinstance(value, (int, float)) == False:
//             raise TypeError('Incorrect type of ' + smth + '. Type should be int or float.')
            
//     if edge is None:
//         if isinstance(value, (int, float)):
//             if value < 0:
//                 raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0.')
//     else:
//         if isinstance(value, (int, float)):
//             if value < 0 or value > edge:
//                 raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0 and equal or less ' + str(edge) + '.')

void CheckIndex(int64_t index, uint64_t max_index, std::string smth) {
    if (!(index == -1 || (index >= 0 && static_cast<uint64_t>(index) < max_index))) {
        throw std::runtime_error("There are no such " + smth + "!");
    }
}

// def check_index(self, index, edge, smth, hap=False, none=True):
//     if none == False and index is None:
//         raise TypeError('Incorrect type of ' + smth + '. Type should be int.')
//     elif isinstance(index, int):
//         if index < 0 or index >= edge:
//             raise IndexError('There are no such ' + smth + '!')
//     elif isinstance(index, str) and hap:
//         if index.count("A") + index.count("T") + index.count("C") + index.count("G") + index.count("*") \
//         != self.sites:
//             raise ValueError('Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"*\" and length of haplotype should be equal number of mutations sites.')
//     elif index is not None:
//         if hap:
//             raise TypeError('Incorrect type of haplotype. Type should be int or str or None.')
//         else:
//             raise TypeError('Incorrect type of ' + smth + '. Type should be int or None.')

