#include "ParallelHDKP1.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cmath>
#include <K2tree.h>
#include "seq/HDKP1.h"
#include "seq/HDKP2.h"

int main(int argc, char **argv) {
    std::string fileA = (argc > 1) ? argv[1] : "set1Conjunto1T.csv.kt";
    std::string fileB = (argc > 2) ? argv[2] : "set2Conjunto1T.csv.kt";

    std::cout << "========================================================" << std::endl;
    std::cout << "  Parallel_HSD_K2: HDKP1 Benchmark (Sequential vs SYCL)  " << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Loading dataset A: " << fileA << std::endl;
    std::cout << "Loading dataset B: " << fileB << std::endl;

    MREP2 *KA = loadK2treeFromFile(fileA.c_str());
    MREP2 *KB = loadK2treeFromFile(fileB.c_str());

    if (!KA || !KB) {
        std::cerr << "Error: Failed to load one or both k2-tree files." << std::endl;
        return 1;
    }

    // Initialize SYCL queue targeting CPU
    sycl::queue q{sycl::cpu_selector_v};
    std::cout << "SYCL Target Device: " << q.get_device().get_info<sycl::info::device::name>() << std::endl;
    std::cout << "Warming up SYCL runtime outside execution timer..." << std::endl;
    warmupSYCL(q);

    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << "1. Directed HDKP1: Sequential vs Parallel SYCL" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    // --- Sequential HDKP1 ---
    auto t1_seq = std::chrono::high_resolution_clock::now();
    double h_seq = hdkp1(KA, KB);
    auto t2_seq = std::chrono::high_resolution_clock::now();
    double ms_seq = std::chrono::duration<double, std::milli>(t2_seq - t1_seq).count();

    std::cout << "   -> Sequential HDKP1: h(A, B) = " << std::fixed << std::setprecision(5) << h_seq
              << " [Host time: " << std::setprecision(3) << ms_seq << " ms]" << std::endl;

    // --- Parallel HDKP1 (SYCL) ---
    auto t1_par = std::chrono::high_resolution_clock::now();
    double h_par = parallelHDKP1(KA, KB, &q);
    auto t2_par = std::chrono::high_resolution_clock::now();
    double ms_par = std::chrono::duration<double, std::milli>(t2_par - t1_par).count();

    std::cout << "   -> Parallel   HDKP1: h(A, B) = " << std::fixed << std::setprecision(5) << h_par
              << " [Host time: " << std::setprecision(3) << ms_par << " ms]" << std::endl;

    // Verify equality
    double diff_directed = std::abs(h_seq - h_par);
    bool directed_equal = (diff_directed < 1e-6);
    double speedup_directed = (ms_par > 0.0) ? (ms_seq / ms_par) : 0.0;

    std::cout << "   -> Verification: " 
              << (directed_equal ? "PASSED (results are equal, diff = " : "FAILED (diff = ")
              << diff_directed << ")" << std::endl;
    std::cout << "   -> Speedup: " << std::setprecision(2) << speedup_directed << "x" << std::endl;

    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << "2. Symmetric HDKP1: max(h(A, B), h(B, A))" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    // --- Symmetric Sequential HDKP1 ---
    t1_seq = std::chrono::high_resolution_clock::now();
    double hSym_seq = symmetricHDKP1(KA, KB);
    t2_seq = std::chrono::high_resolution_clock::now();
    double ms_sym_seq = std::chrono::duration<double, std::milli>(t2_seq - t1_seq).count();

    std::cout << "   -> Sequential Symmetric HDKP1 = " << std::fixed << std::setprecision(5) << hSym_seq
              << " [Host time: " << std::setprecision(3) << ms_sym_seq << " ms]" << std::endl;

    // --- Symmetric Parallel HDKP1 ---
    t1_par = std::chrono::high_resolution_clock::now();
    double hSym_par = symmetricParallelHDKP1(KA, KB, &q);
    t2_par = std::chrono::high_resolution_clock::now();
    double ms_sym_par = std::chrono::duration<double, std::milli>(t2_par - t1_par).count();

    std::cout << "   -> Parallel   Symmetric HDKP1 = " << std::fixed << std::setprecision(5) << hSym_par
              << " [Host time: " << std::setprecision(3) << ms_sym_par << " ms]" << std::endl;

    // Verify equality
    double diff_sym = std::abs(hSym_seq - hSym_par);
    bool sym_equal = (diff_sym < 1e-6);
    double speedup_sym = (ms_sym_par > 0.0) ? (ms_sym_seq / ms_sym_par) : 0.0;

    std::cout << "   -> Verification: "
              << (sym_equal ? "PASSED (results are equal, diff = " : "FAILED (diff = ")
              << diff_sym << ")" << std::endl;
    std::cout << "   -> Speedup: " << std::setprecision(2) << speedup_sym << "x" << std::endl;

    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << "3. Sequential HDKP2 Baseline (for reference)" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    t1_seq = std::chrono::high_resolution_clock::now();
    double h2 = hdkp2(KA, KB);
    t2_seq = std::chrono::high_resolution_clock::now();
    double ms_h2 = std::chrono::duration<double, std::milli>(t2_seq - t1_seq).count();

    std::cout << "   -> HDKP2 h(A, B) = " << std::fixed << std::setprecision(5) << h2
              << " [Host time: " << std::setprecision(3) << ms_h2 << " ms]" << std::endl;

    destroyK2tree(KA);
    destroyK2tree(KB);

    std::cout << "\n========================================================" << std::endl;
    if (directed_equal && sym_equal) {
        std::cout << "  ALL VERIFICATION CHECKS PASSED SUCCESSFULLY           " << std::endl;
    } else {
        std::cout << "  WARNING: SOME VERIFICATION CHECKS FAILED              " << std::endl;
    }
    std::cout << "========================================================" << std::endl;

    return (directed_equal && sym_equal) ? 0 : 1;
}
