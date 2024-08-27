#include <iostream>

#include "simulator.cpp"

// clang++ main.cpp -std=c++17 -O2 -Wall -Werror -Wsign-compare -o fast_solution
// g++ main.cpp -std=c++17 -O2 -Wall -Werror -Wsign-compare -o fast_solution
// clang++ main.cpp -fsanitize=address,undefined -fno-sanitize-recover=all -std=c++17 -Wall -Werror -Wsign-compare -o debug_solution

// clang++ wrap.cpp -std=c++17 -O2 -Wall -Werror -Wsign-compare -o fast_solution

int main() {
    Simulator simulator(2, 1, 3, 1239);
    simulator.Simulate(100);
    simulator.Genealogy();
    // simulator.Debug();
    return 0;
}