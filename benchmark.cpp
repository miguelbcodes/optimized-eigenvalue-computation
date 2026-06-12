#include "include/matrix.h"
#include "include/hessenberg.h"
#include "include/qr.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>

static long get_current_rss_kb() {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            long kb = 0;
            for (char c : line)
                if (c >= '0' && c <= '9') kb = kb * 10 + (c - '0');
            return kb;
        }
    }
    return 0;
}

int main() {
    const int sizes[] = {10, 50, 100, 250, 500, 1000};
    const int num_sizes = 6;
    const double tolerance = 1e-10;

    std::cout << "algorithm,N,time_ms,memory_kb" << std::endl;

    for (int s = 0; s < num_sizes; ++s) {
        int n = sizes[s];

        long mem_before = get_current_rss_kb();

        Matrix A(n);
        Matrix::fill_random(A);

        auto start = std::chrono::high_resolution_clock::now();
        Matrix H = reduce_to_hessenberg(A);
        Matrix schur = qr_iterate(H, tolerance);
        auto end = std::chrono::high_resolution_clock::now();

        long mem_after = get_current_rss_kb();
        long mem_used = mem_after - mem_before;
        if (mem_used < 0) mem_used = 0;

        double elapsed = std::chrono::duration<double, std::milli>(end - start).count();

        std::cout << "hessenberg_qr," << n << ","
                  << std::fixed << std::setprecision(2) << elapsed << ","
                  << mem_used << std::endl;

        std::cerr << "N=" << n << " done (" << std::fixed << std::setprecision(2) << elapsed << " ms)" << std::endl;
    }

    return 0;
}
