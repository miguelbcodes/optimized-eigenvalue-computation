#include "../include/matrix.h"
#include "../include/hessenberg.h"
#include "../include/qr.h"
#include "../include/tridiagonal.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <malloc.h>

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

static void emit_row(const char* algorithm, int n, double elapsed_ms, long memory_kb) {
    std::cout << algorithm << "," << n << ","
              << std::fixed << std::setprecision(2) << elapsed_ms << ","
              << memory_kb << std::endl;
}

int main(int argc, char* argv[]) {
    const int sizes[] = {10, 50, 100, 250, 500, 1000};
    const int num_sizes = 6;
    const double tolerance = 1e-10;
    const int samples = (argc > 1) ? std::atoi(argv[1]) : 1;

    std::cout << "algorithm,N,time_ms,memory_kb" << std::endl;

    for (int sample = 0; sample < samples; ++sample) {
    for (int s = 0; s < num_sizes; ++s) {
        int n = sizes[s];

        // The same symmetric matrix feeds both algorithms for a fair comparison.
        // It is allocated before either measurement, so it sits in the shared
        // baseline and cancels out of each algorithm's memory delta.
        Matrix A(n);
        Matrix::fill_random_symmetric(A);

        // General algorithm: dense Hessenberg reduction + QR iteration, O(N^3).
        // H and schur (both O(N^2)) are alive when memory is sampled.
        {
            long mem_before = get_current_rss_kb();
            auto start = std::chrono::high_resolution_clock::now();
            Matrix H = reduce_to_hessenberg(A);
            Matrix schur = qr_iterate(H, tolerance);
            auto end = std::chrono::high_resolution_clock::now();
            long mem_used = get_current_rss_kb() - mem_before;
            if (mem_used < 0) mem_used = 0;
            double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
            emit_row("hessenberg_qr", n, elapsed, mem_used);
        }
        // Return the dense temporaries to the OS so the next baseline is clean.
        malloc_trim(0);

        // Optimized algorithm: symmetric tridiagonal reduction + implicit-shift
        // QL, with O(N) bands and an O(N^2) QR phase. Only O(N) state (T, eig) is
        // alive when memory is sampled.
        {
            long mem_before = get_current_rss_kb();
            auto start = std::chrono::high_resolution_clock::now();
            SymmetricTridiagonal T = reduce_to_tridiagonal(A);
            std::vector<double> eig = qr_iterate_tridiagonal(T, tolerance);
            auto end = std::chrono::high_resolution_clock::now();
            long mem_used = get_current_rss_kb() - mem_before;
            if (mem_used < 0) mem_used = 0;
            double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
            emit_row("tridiagonal_qr", n, elapsed, mem_used);
        }
        malloc_trim(0);

        std::cerr << "sample=" << (sample+1) << "/" << samples << " N=" << n << " done" << std::endl;
    }
    }

    return 0;
}
