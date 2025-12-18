#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <string>
#include <chrono>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <htslib/sam.h>
#include "kll_sketch.hpp"


using namespace datasketches;
using hr_clock = std::chrono::steady_clock;

/* ===============================
   Utilidades baseline exacto
   =============================== */

float exact_quantile(const std::vector<uint32_t>& v, double q) {
    size_t idx = static_cast<size_t>(q * (v.size() - 1));
    return static_cast<float>(v[idx]);
}

double rank_error(
    const std::vector<uint32_t>& v,
    float value,
    double q
) {
    auto it = std::lower_bound(v.begin(), v.end(), value);
    size_t rank_kll = std::distance(v.begin(), it);
    size_t rank_exact = static_cast<size_t>(q * v.size());
    return std::abs((double)rank_kll - (double)rank_exact) / v.size();
}

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Uso: " << argv[0]
                  << " <archivo.bam> <bin_size> <output.csv>\n";
        return 1;
    }

    const char* bam_file = argv[1];
    int bin_size = std::atoi(argv[2]);
    const char* csv_file = argv[3];

    /* ===============================
       1️⃣ BASELINE EXACTO
       =============================== */

    samFile* bam_fp = sam_open(bam_file, "r");
    sam_hdr_t* header = sam_hdr_read(bam_fp);
    bam1_t* aln = bam_init1();

    std::string current_chr;
    std::unordered_map<uint32_t, uint32_t> bins;
    bins.reserve(100000);

    std::vector<uint32_t> exact_values;

    while (sam_read1(bam_fp, header, aln) >= 0) {

        if (aln->core.flag & (BAM_FUNMAP |
                              BAM_FSECONDARY |
                              BAM_FSUPPLEMENTARY |
                              BAM_FDUP))
            continue;

        std::string chr = header->target_name[aln->core.tid];

        if (current_chr.empty())
            current_chr = chr;

        if (chr != current_chr) {
            for (const auto& e : bins)
                exact_values.push_back(e.second);
            bins.clear();
            current_chr = chr;
        }

        int start = aln->core.pos;
        int end   = bam_endpos(aln);

        for (int p = start; p < end; p += bin_size)
            bins[p / bin_size]++;
    }

    for (const auto& e : bins)
        exact_values.push_back(e.second);

    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(bam_fp);

    std::sort(exact_values.begin(), exact_values.end());
    size_t N = exact_values.size();

    /* ===============================
       CSV header
       =============================== */

    std::ofstream out(csv_file);
    out << std::fixed << std::setprecision(6);
    out << "bin_size,K,"
        << "p50_kll,p50_exact,p50_rank_error,"
        << "p95_kll,p95_exact,p95_rank_error,"
        << "kll_time_sec,kll_bytes\n";

    /* ===============================
       2️⃣ EXPERIMENTO KLL
       =============================== */

    for (int K : {100, 200, 300, 400, 500, 1000, 2000}) {

        kll_sketch<float> sketch(K);
        auto t1 = hr_clock::now();

        for (uint32_t v : exact_values)
            sketch.update(static_cast<float>(v));

        auto t2 = hr_clock::now();
        double kll_time = std::chrono::duration<double>(t2 - t1).count();

        float p5_kll = sketch.get_quantile(0.05);
        float p50_kll = sketch.get_quantile(0.50);
        float p95_kll = sketch.get_quantile(0.95);

        float p5_exact = exact_quantile(exact_values, 0.05);
        float p50_exact = exact_quantile(exact_values, 0.50);
        float p95_exact = exact_quantile(exact_values, 0.95);

        double p5_err = rank_error(exact_values, p5_kll, 0.05);
        double p50_err = rank_error(exact_values, p50_kll, 0.50);
        double p95_err = rank_error(exact_values, p95_kll, 0.95);
        

        size_t bytes = sketch.get_serialized_size_bytes();

        out << bin_size << "," << K << ","
            << p5_kll << "," << p5_exact << "," << p5_err << ","
            << p50_kll << "," << p50_exact << "," << p50_err << ","
            << p95_kll << "," << p95_exact << "," << p95_err << ","
            << kll_time << "," << bytes << "\n";
    }

    out.close();
    return 0;
}
