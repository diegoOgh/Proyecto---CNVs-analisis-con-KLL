#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <algorithm>

#include <htslib/sam.h>

/*
 * ============================
 *   Estructuras de datos
 * ============================
 */

struct BaselineStats {
    int bin_size;
    float p25;
    float p50;
    float p75;
    float iqr;
    float deletion_threshold;
    float duplication_threshold;
};

struct CNV {
    std::string chr;
    uint64_t start;
    uint64_t end;
    std::string type;  // "DEL" o "DUP"
    float mean_coverage;
    uint64_t num_bins;
};

/*
 * ============================
 *   Utilidades
 * ============================
 */

BaselineStats load_baseline(const std::string& csv_file, int bin_size) {
    std::ifstream in(csv_file);
    if (!in.is_open())
        throw std::runtime_error("No se pudo abrir baseline CSV");

    std::string line;
    std::getline(in, line); // header

    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string field;
        BaselineStats b{};

        // bin_size
        std::getline(ss, field, ',');
        b.bin_size = std::stoi(field);
        if (b.bin_size != bin_size)
            continue;

        // total_bins (skip)
        std::getline(ss, field, ',');

        // p1, p5 (skip)
        std::getline(ss, field, ',');
        std::getline(ss, field, ',');

        // p25, p50, p75
        std::getline(ss, field, ','); b.p25 = std::stof(field);
        std::getline(ss, field, ','); b.p50 = std::stof(field);
        std::getline(ss, field, ','); b.p75 = std::stof(field);

        // p95, p99 (skip)
        std::getline(ss, field, ',');
        std::getline(ss, field, ',');

        // min, max (skip)
        std::getline(ss, field, ',');
        std::getline(ss, field, ',');

        // iqr
        std::getline(ss, field, ','); b.iqr = std::stof(field);

        // deletion_threshold, duplication_threshold
        std::getline(ss, field, ','); b.deletion_threshold = std::stof(field);
        std::getline(ss, field, ','); b.duplication_threshold = std::stof(field);

        return b;
    }

    throw std::runtime_error("Bin size no encontrado en baseline CSV");
}


/*
 * ============================
 *   Detección de CNVs
 * ============================
 */

void detect_cnvs_for_chr(
    const std::string& chr,
    const std::unordered_map<uint32_t, uint32_t>& bins,
    const BaselineStats& base,
    std::vector<CNV>& cnvs
) {
    uint32_t run_start_bin = 0;
    uint32_t run_len = 0;
    uint64_t run_sum = 0;
    std::string run_type;

    auto flush_run = [&](uint32_t end_bin) {
        if (run_len == 0) return;

        CNV cnv;
        cnv.chr = chr;
        cnv.start = uint64_t(run_start_bin) * base.bin_size;
        cnv.end   = uint64_t(end_bin + 1) * base.bin_size;
        cnv.type  = run_type;
        cnv.num_bins = run_len;
        cnv.mean_coverage = float(run_sum) / run_len;

        cnvs.push_back(cnv);

        run_len = 0;
        run_sum = 0;
        run_type.clear();
    };

    std::vector<std::pair<uint32_t, uint32_t>> sorted_bins(
    bins.begin(), bins.end()
);

    std::sort(sorted_bins.begin(), sorted_bins.end(),
            [](auto& a, auto& b) {
                return a.first < b.first;
            });


    for (const auto& e : bins) {
        uint32_t bin = e.first;
        uint32_t cov = e.second;

        std::string current_type;

        if (cov < base.deletion_threshold)
            current_type = "DEL";
        else if (cov > base.duplication_threshold)
            current_type = "DUP";
        else {
            flush_run(bin - 1);
            continue;
        }

        if (run_len == 0) {
            run_start_bin = bin;
            run_len = 1;
            run_sum = cov;
            run_type = current_type;
        }
        else if (current_type == run_type && bin == run_start_bin + run_len) {
            run_len++;
            run_sum += cov;
        }
        else {
            flush_run(bin - 1);
            run_start_bin = bin;
            run_len = 1;
            run_sum = cov;
            run_type = current_type;
        }
    }

    if (run_len > 0)
        flush_run(run_start_bin + run_len - 1);
}

/*
 * ============================
 *   MAIN
 * ============================
 */

int main(int argc, char* argv[]) {

    if (argc != 6) {
        std::cerr << "Uso: " << argv[0]
                  << " <archivo.bam> <bin_size> <baseline.csv> <output_cnvs.csv> <min_bins>\n";
        return 1;
    }

    const char* bam_file = argv[1];
    int bin_size = std::atoi(argv[2]);
    std::string baseline_csv = argv[3];
    std::string output_csv = argv[4];
    int min_bins = std::atoi(argv[5]);

    // --- Leer baseline ---
    BaselineStats base = load_baseline(baseline_csv, bin_size);

    // --- Abrir BAM ---
    samFile* bam_fp = sam_open(bam_file, "r");
    if (!bam_fp) {
        std::cerr << "Error abriendo BAM\n";
        return 1;
    }

    sam_hdr_t* header = sam_hdr_read(bam_fp);
    bam1_t* aln = bam_init1();

    std::string current_chr;
    std::unordered_map<uint32_t, uint32_t> bins;
    bins.reserve(100000);

    std::vector<CNV> cnvs;

    // --- Lectura BAM ---
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
            detect_cnvs_for_chr(current_chr, bins, base, cnvs);
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
    detect_cnvs_for_chr(current_chr, bins, base, cnvs);

    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(bam_fp);

    // --- Guardar CNVs ---
    std::ofstream out(output_csv);
    out << "chr,start,end,type,mean_coverage,num_bins\n";

    for (const auto& c : cnvs) {
        if (c.num_bins < (uint64_t)min_bins)
            continue;

        out << c.chr << ","
            << c.start << ","
            << c.end << ","
            << c.type << ","
            << std::fixed << std::setprecision(2)
            << c.mean_coverage << ","
            << c.num_bins << "\n";
    }

    

    out.close();

    return 0;
}
