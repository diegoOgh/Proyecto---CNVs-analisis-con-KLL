#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <string>
#include <chrono>
#include <fstream>

#include <htslib/sam.h>
#include "datasketches-cpp/kll/include/kll_sketch.hpp"

using namespace datasketches;

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Uso: " << argv[0]
                  << " <archivo.bam> <bin_size> <output.csv>\n";
        return 1;
    }

    const char* bam_file = argv[1];
    int bin_size = std::atoi(argv[2]);
    const char* csv_file = argv[3];

    // --- Abrir BAM ---
    samFile* bam_fp = sam_open(bam_file, "r");
    if (!bam_fp) {
        std::cerr << "Error al abrir BAM\n";
        return 1;
    }

    sam_hdr_t* header = sam_hdr_read(bam_fp);
    if (!header) {
        std::cerr << "Error leyendo header\n";
        sam_close(bam_fp);
        return 1;
    }

    bam1_t* aln = bam_init1();

    // --- KLL ---
    constexpr int K = 200;
    kll_sketch<float> coverage_sketch(K);

    using clock = std::chrono::steady_clock;
    std::chrono::duration<double> kll_time{0};

    std::string current_chr;
    std::unordered_map<uint32_t, uint32_t> bins;
    bins.reserve(100000);

    uint64_t total_reads = 0;
    uint64_t total_bins = 0;

    for (int i = 0; i < header->n_targets; ++i) {
    uint64_t chr_len = header->target_len[i];
    total_bins += (chr_len + bin_size - 1) / bin_size;
}



    // --- Lectura BAM ---
    while (sam_read1(bam_fp, header, aln) >= 0) {

        total_reads++;

        if (aln->core.flag & (BAM_FUNMAP |
                              BAM_FSECONDARY |
                              BAM_FSUPPLEMENTARY |
                              BAM_FDUP))
            continue;

        std::string chr = header->target_name[aln->core.tid];

        if (current_chr.empty())
            current_chr = chr;

        if (chr != current_chr) {
            auto t1 = clock::now();
            for (const auto& e : bins)
                coverage_sketch.update(static_cast<float>(e.second));
            auto t2 = clock::now();
            kll_time += (t2 - t1);

            bins.clear();
            current_chr = chr;
        }

        int start = aln->core.pos;
        int end   = bam_endpos(aln);

        for (int p = start; p < end; p += bin_size) {
            uint32_t bin = p / bin_size;
            bins[bin]++;
        }
    }

    // Flush final
    auto t1 = clock::now();
    for (const auto& e : bins)
        coverage_sketch.update(static_cast<float>(e.second));
    auto t2 = clock::now();
    kll_time += (t2 - t1);

    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(bam_fp);

    // --- Percentiles EXACTOS ---
    float p1  = coverage_sketch.get_quantile(0.01);
    float p5  = coverage_sketch.get_quantile(0.05);
    float p25 = coverage_sketch.get_quantile(0.25);
    float p50 = coverage_sketch.get_quantile(0.50);
    float p75 = coverage_sketch.get_quantile(0.75);
    float p95 = coverage_sketch.get_quantile(0.95);
    float p99 = coverage_sketch.get_quantile(0.99);

    float minv = coverage_sketch.get_min_item();
    float maxv = coverage_sketch.get_max_item();

    float iqr = p75 - p25;
    float deletion_threshold    = p50 * 0.5f;
    float duplication_threshold = p50 * 1.5f;

    size_t kll_items = coverage_sketch.get_num_retained();

    // --- CSV ---
    bool write_header = false;
    std::ifstream check(csv_file);
    if (!check.good())
        write_header = true;
    check.close();

    std::ofstream out(csv_file, std::ios::app);
    out << std::fixed << std::setprecision(6);

    if (write_header) {
        out << "bin_size,"
            << "total_bins," 
            << "p1,p5,p25,p50,p75,p95,p99,"
            << "min,max,iqr,"
            << "deletion_threshold,duplication_threshold,"
            << "kll_items,kll_k,kll_time_sec\n";
    }

    out << bin_size << ","
        << total_bins << ","
        << p1 << "," << p5 << "," << p25 << "," << p50 << ","
        << p75 << "," << p95 << "," << p99 << ","
        << minv << "," << maxv << ","
        << iqr << ","
        << deletion_threshold << ","
        << duplication_threshold << ","
        << kll_items << ","
        << K << ","
        << kll_time.count() << "\n";

    out.close();

    return 0;
}
