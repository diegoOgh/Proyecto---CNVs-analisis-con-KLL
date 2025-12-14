#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <string>
#include <chrono>

#include <htslib/sam.h>
#include "datasketches-cpp/kll/include/kll_sketch.hpp"

using namespace datasketches;

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << "Uso: " << argv[0] << " <archivo.bam> <bin_size>\n";
        return 1;
    }

    const char* bam_file = argv[1];
    int bin_size = std::atoi(argv[2]);

    std::cout << "Analizando cobertura en bins de "
              << bin_size << " bp\n";

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

    // --- Alineamientos ---
    bam1_t* aln = bam_init1();

    // --- KLL Sketch ---
    kll_sketch<float> coverage_sketch(200);
    using clock = std::chrono::steady_clock;
    std::chrono::duration<double> kll_time{0};

    // --- Estado por cromosoma ---
    std::string current_chr = "";
    std::unordered_map<uint32_t, uint32_t> bins;
    bins.reserve(100000); // evita rehash

    uint64_t total_reads = 0;

    std::cout << "Procesando reads...\n";

    // Leer BAM y contar cobertura por bin
    while (sam_read1(bam_fp, header, aln) >= 0) {

        total_reads++;

        //cada 5 millones de reads procesados, imprimir progreso
        if (total_reads % 100'000'000 == 0)
            std::cout << "Procesados: " << total_reads << "\n";

        // Filtros que saltan reads que no deben contar para cobertura
        // Donde BAM_FUNMAP = read no mapeado
        //       BAM_FSECONDARY = read secundario
        //       BAM_FSUPPLEMENTARY = read suplementario
        //       BAM_FDUP = read duplicado
        if (aln->core.flag & (BAM_FUNMAP |
                              BAM_FSECONDARY |
                              BAM_FSUPPLEMENTARY |
                              BAM_FDUP))
            continue;

        // Obtener cromosoma actual
        std::string chr = header->target_name[aln->core.tid];

        // Inicialización
        if (current_chr.empty())
            current_chr = chr;

        // Cambio de cromosoma → flush
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

        //Recorrer posiciones del read y actualizar bins
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


    // Limpieza
    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(bam_fp);

    std::cout << "Total reads procesados: "
              << total_reads << "\n";

    // --- Análisis de distribución de cobertura ---
std::cout << std::fixed << std::setprecision(2);
std::cout << "\n=== DISTRIBUCIÓN DE COBERTURA (por bin de "
          << bin_size << " bp) ===\n";

float p1  = coverage_sketch.get_quantile(0.01);
float p5  = coverage_sketch.get_quantile(0.05);
float p25 = coverage_sketch.get_quantile(0.25);
float p50 = coverage_sketch.get_quantile(0.50);
float p75 = coverage_sketch.get_quantile(0.75);
float p95 = coverage_sketch.get_quantile(0.95);
float p99 = coverage_sketch.get_quantile(0.99);

std::cout << "P1  (percentil 1):  " << p1  << "×\n";
std::cout << "P5  (percentil 5):  " << p5  << "×\n";
std::cout << "P25 (Q1):           " << p25 << "×\n";
std::cout << "P50 (Mediana):      " << p50 << "× ← Cobertura típica\n";
std::cout << "P75 (Q3):           " << p75 << "×\n";
std::cout << "P95:                " << p95 << "×\n";
std::cout << "P99:                " << p99 << "×\n";
std::cout << "Mínimo:             " << coverage_sketch.get_min_item() << "×\n";
std::cout << "Máximo:             " << coverage_sketch.get_max_item() << "×\n";

// --- Detección de anomalías (heurística simple) ---
std::cout << "\n=== DETECCIÓN DE POSIBLES CNVs (heurística) ===\n";

float iqr = p75 - p25;
float low_threshold_iqr  = p25 - 1.5f * iqr;
float high_threshold_iqr = p75 + 1.5f * iqr;

// Umbrales basados en la mediana
float deletion_threshold    = p50 * 0.5f;
float duplication_threshold = p50 * 1.5f;

std::cout << "Cobertura normal (IQR): "
          << p25 << "× - " << p75 << "×\n";
std::cout << "Umbral deleción (< 0.5× mediana): "
          << deletion_threshold << "×\n";
std::cout << "Umbral duplicación (> 1.5× mediana): "
          << duplication_threshold << "×\n";

// --- Métricas de compresión ---
std::cout << "\n=== MÉTRICAS DE COMPRESIÓN ===\n";
std::cout << "Items retenidos en KLL: "
          << coverage_sketch.get_num_retained() << "\n";
std::cout << "Parámetro k: 200\n";

// --- Tiempo KLL ---
std::cout << "\n=== RENDIMIENTO KLL ===\n";
std::cout << "Tiempo total KLL (resumen de cobertura): "
          << kll_time.count() << " segundos\n";

    return 0;
}