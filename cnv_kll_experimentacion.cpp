#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <string>
#include <chrono>

#include <htslib/sam.h>
#include "datasketches-cpp/kll/include/kll_sketch.hpp"

using namespace datasketches;
using hr_clock = std::chrono::high_resolution_clock;

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <archivo.bam>\n";
        return 1;
    }

    const char* bam_file = argv[1];

    std::vector<int> bin_sizes = {100 ,200, 500, 1000, 2000, 5000, 10000};
    const int K = 200;

    std::ofstream csv("bin_experiment.csv");
    csv << "bin_size,num_bins,p25,p50,p75,p95,"
           "kll_items,kll_k,kll_time_sec,kll_memory_bytes\n";

    for (int bin_size : bin_sizes) {

        std::cout << "\n=== BIN SIZE: " << bin_size << " bp ===\n";

        // --- Abrir BAM ---
        samFile* bam_fp = sam_open(bam_file, "r");
        if (!bam_fp) {
            std::cerr << "Error abriendo BAM\n";
            return 1;
        }

        sam_hdr_t* header = sam_hdr_read(bam_fp);
        bam1_t* aln = bam_init1();

        kll_sketch<float> coverage_sketch(K);
        std::unordered_map<uint32_t, uint32_t> bins;
        bins.reserve(200000);

        std::string current_chr;

        while (sam_read1(bam_fp, header, aln) >= 0) {

            if (aln->core.flag & (BAM_FUNMAP |
                                  BAM_FSECONDARY |
                                  BAM_FSUPPLEMENTARY))
                continue;

            std::string chr = header->target_name[aln->core.tid];

            if (current_chr.empty())
                current_chr = chr;

            if (chr != current_chr) {
                for (const auto& e : bins)
                    coverage_sketch.update((float)e.second);
                bins.clear();
                current_chr = chr;
            }

            int start = aln->core.pos;
            int end   = bam_endpos(aln);

            for (int p = start; p < end; p += bin_size) {
                bins[p / bin_size]++;
            }
        }

        // --- Timing KLL ---
        auto t1 = hr_clock::now();
        for (const auto& e : bins)
            coverage_sketch.update((float)e.second);
        auto t2 = hr_clock::now();

        double kll_time =
            std::chrono::duration<double>(t2 - t1).count();

        float p25 = coverage_sketch.get_quantile(0.25);
        float p50 = coverage_sketch.get_quantile(0.50);
        float p75 = coverage_sketch.get_quantile(0.75);
        float p95 = coverage_sketch.get_quantile(0.95);

        size_t kll_items = coverage_sketch.get_num_retained();
        size_t kll_mem   = coverage_sketch.get_serialized_size_bytes();

        csv << bin_size << ","
            << bins.size() << ","
            << p25 << ","
            << p50 << ","
            << p75 << ","
            << p95 << ","
            << kll_items << ","
            << K << ","
            << kll_time << ","
            << kll_mem << "\n";

        std::cout << "Mediana: " << p50 << "×\n";
        std::cout << "KLL items: " << kll_items << "\n";
        std::cout << "Memoria KLL: " << kll_mem / 1024.0 << " KB\n";
        std::cout << "Tiempo KLL: " << kll_time << " s\n";

        bam_destroy1(aln);
        sam_hdr_destroy(header);
        sam_close(bam_fp);
    }

    csv.close();
    std::cout << "\nExperimento terminado → bin_experiment.csv\n";
    return 0;
}
