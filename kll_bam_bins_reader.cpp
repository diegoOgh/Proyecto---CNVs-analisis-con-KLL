#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include "datasketches-cpp/kll/include/kll_sketch.hpp"
#include <htslib/sam.h>

using namespace datasketches;

// struct para los binsS
struct GenomicBin {
    std::string chrom;
    int start;
    int end;
    uint32_t coverage;
    
    GenomicBin(std::string c, int s, int e) 
        : chrom(c), start(s), end(e), coverage(0) {}
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Uso: " << argv[0] << " <archivo.bam> <bin_size>" << std::endl;
        return 1;
    }
    
    const char* bam_file = argv[1];
    int bin_size = std::atoi(argv[2]); // ejemplo:  10000 = 10kb bins
    
    std::cout << "Analizando cobertura en bins de " << bin_size << " bp..." << std::endl;
    
    // Abrir BAM
    samFile* bam_fp = sam_open(bam_file, "r");
    if (!bam_fp) {
        std::cerr << "Error al abrir " << bam_file << std::endl;
        return 1;
    }
    
    sam_hdr_t* header = sam_hdr_read(bam_fp);
    if (!header) {
        std::cerr << "Error al leer header" << std::endl;
        sam_close(bam_fp);
        return 1;
    }
    
    // Mapa para almacenar cobertura por bin
    // Key: "chr1:0", Value: coverage count
    std::map<std::string, uint32_t> bin_coverage;
    
    // Para generar las keys de bins
    auto get_bin_key = [&](const std::string& chr, int pos) -> std::string {
        int bin_start = (pos / bin_size) * bin_size;
        return chr + ":" + std::to_string(bin_start);
    };
    
    // Alineamiento
    bam1_t* aln = bam_init1();
    uint64_t total_reads = 0;
    
    std::cout << "Procesando reads..." << std::endl;
    
    // Leer BAM y contar cobertura por bin
    while (sam_read1(bam_fp, header, aln) >= 0) {
        total_reads++;
        
        if (total_reads % 5000000 == 0) {
            std::cout << "Procesados: " << total_reads << " reads..." << std::endl;
        }
        
        // Skip unmapped
        if (aln->core.flag & BAM_FUNMAP) continue;
        
        // Skip secondary/supplementary
        if (aln->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        
        // Obtener cromosoma y posición
        std::string chrom = header->target_name[aln->core.tid];
        int pos = aln->core.pos;
        int end_pos = bam_endpos(aln);
        
        // Incrementar cobertura en todos los bins que cubre este read
        for (int p = pos; p < end_pos; p += bin_size) {
            std::string bin_key = get_bin_key(chrom, p);
            bin_coverage[bin_key]++;
        }
    }
    
    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(bam_fp);
    
    std::cout << "Total reads procesados: " << total_reads << std::endl;
    std::cout << "Total bins con cobertura: " << bin_coverage.size() << std::endl;
    
    // --- Insertar coberturas en KLL sketch ---
    kll_sketch<float> coverage_sketch(200);
    
    std::cout << "\nInsertando coberturas en KLL sketch..." << std::endl;
    
    for (const auto& entry : bin_coverage) {
        coverage_sketch.update(static_cast<float>(entry.second));
    }
    
    // --- Análisis de distribución de cobertura ---
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n=== DISTRIBUCIÓN DE COBERTURA (por bin de " << bin_size << " bp) ===" << std::endl;
    
    float p1 = coverage_sketch.get_quantile(0.01);
    float p5 = coverage_sketch.get_quantile(0.05);
    float p25 = coverage_sketch.get_quantile(0.25);
    float p50 = coverage_sketch.get_quantile(0.50);
    float p75 = coverage_sketch.get_quantile(0.75);
    float p95 = coverage_sketch.get_quantile(0.95);
    float p99 = coverage_sketch.get_quantile(0.99);
    
    std::cout << "P1  (percentil 1):  " << p1 << "×" << std::endl;
    std::cout << "P5  (percentil 5):  " << p5 << "×" << std::endl;
    std::cout << "P25 (Q1):           " << p25 << "×" << std::endl;
    std::cout << "P50 (Mediana):      " << p50 << "× ← Cobertura típica" << std::endl;
    std::cout << "P75 (Q3):           " << p75 << "×" << std::endl;
    std::cout << "P95:                " << p95 << "×" << std::endl;
    std::cout << "P99:                " << p99 << "×" << std::endl;
    std::cout << "Mínimo:             " << coverage_sketch.get_min_item() << "×" << std::endl;
    std::cout << "Máximo:             " << coverage_sketch.get_max_item() << "×" << std::endl;
    
    // --- Detección de anomalías (heurística simple) ---
    std::cout << "\n=== DETECCIÓN DE POSIBLES CNVs (heurística) ===" << std::endl;
    
    float iqr = p75 - p25;
    float low_threshold = p25 - 1.5 * iqr;   // Posibles deleciones
    float high_threshold = p75 + 1.5 * iqr;  // Posibles duplicaciones
    
    // Alternativamente, usar 2× y 0.5× de la mediana
    float duplication_threshold = p50 * 1.5;  // 1.5× mediana
    float deletion_threshold = p50 * 0.5;     // 0.5× mediana
    
    std::cout << "Cobertura normal: " << p25 << "× - " << p75 << "×" << std::endl;
    std::cout << "Umbral deleción (< 0.5× mediana): " << deletion_threshold << "×" << std::endl;
    std::cout << "Umbral duplicación (> 1.5× mediana): " << duplication_threshold << "×" << std::endl;
    
    // Contar bins anómalos
    uint32_t low_bins = 0, high_bins = 0;
    
    for (const auto& entry : bin_coverage) {
        if (entry.second < deletion_threshold) low_bins++;
        if (entry.second > duplication_threshold) high_bins++;
    }
    
    std::cout << "\nBins con baja cobertura (posibles deleciones): " << low_bins 
              << " (" << (100.0 * low_bins / bin_coverage.size()) << "%)" << std::endl;
    std::cout << "Bins con alta cobertura (posibles duplicaciones): " << high_bins 
              << " (" << (100.0 * high_bins / bin_coverage.size()) << "%)" << std::endl;
    
    // --- Métricas adicionales ---
    std::cout << "\n=== MÉTRICAS DE COMPRESIÓN ===" << std::endl;
    std::cout << "Bins totales: " << bin_coverage.size() << std::endl;
    std::cout << "Items en sketch: " << coverage_sketch.get_num_retained() << std::endl;
    std::cout << "Compresión: " << (float)bin_coverage.size() / coverage_sketch.get_num_retained() << "x" << std::endl;
    
    return 0;
}