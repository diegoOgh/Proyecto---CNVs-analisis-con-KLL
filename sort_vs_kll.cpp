#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>

using hr_clock = std::chrono::high_resolution_clock;

int main() {

    std::ifstream in("bin_experiment.csv");
    if (!in) return 1;

    std::ofstream out("kll_vs_sort_comparison.csv");
    out << "bin_size,num_bins,"
           "kll_time_sec,kll_memory_bytes,"
           "sort_time_sec,sort_memory_bytes\n";

    std::string line;
    std::getline(in, line); // header

    while (std::getline(in, line)) {

        std::stringstream ss(line);
        std::string token;

        int bin_size;
        size_t num_bins;
        double kll_time;
        size_t kll_mem;

        std::getline(ss, token, ','); bin_size = std::stoi(token);
        std::getline(ss, token, ','); num_bins = std::stoull(token);

        // saltar p25,p50,p75,p95
        for (int i = 0; i < 4; ++i)
            std::getline(ss, token, ',');

        std::getline(ss, token, ','); // kll_items
        std::getline(ss, token, ','); // kll_k
        std::getline(ss, token, ','); kll_time = std::stod(token);
        std::getline(ss, token, ','); kll_mem  = std::stoull(token);

        // --- SORT ---
        std::vector<uint32_t> data(num_bins, 100);

        auto t1 = hr_clock::now();
        std::sort(data.begin(), data.end());
        auto t2 = hr_clock::now();

        double sort_time =
            std::chrono::duration<double>(t2 - t1).count();

        size_t sort_mem = num_bins * sizeof(uint32_t);

        out << bin_size << ","
            << num_bins << ","
            << kll_time << ","
            << kll_mem << ","
            << sort_time << ","
            << sort_mem << "\n";
    }

    return 0;
}
