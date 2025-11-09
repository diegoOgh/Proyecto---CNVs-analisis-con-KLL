#include <iostream>
#include <iomanip>
#include <string>
#include "datasketches-cpp/kll/include/kll_sketch.hpp"
#include <htslib/sam.h>

using namespace datasketches;

int main(int argc, char* argv[]) {
    // Verificar argumentos
    if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <archivo.bam>" << std::endl;
        return 1;
    }
    
    const char* bam_file = argv[1];
    
    // Crear sketches para diferentes métricas
    kll_sketch<float> mapping_quality_sketch(200);
    kll_sketch<float> read_length_sketch(200);
    kll_sketch<float> insert_size_sketch(200);
    
    // Abrir archivo BAM
    samFile* bam_fp = sam_open(bam_file, "r");
    if (!bam_fp) {
        std::cerr << "Error: No se pudo abrir " << bam_file << std::endl;
        return 1;
    }
    
    // Leer header
    sam_hdr_t* header = sam_hdr_read(bam_fp);
    if (!header) {
        std::cerr << "Error: No se pudo leer el header del BAM" << std::endl;
        sam_close(bam_fp);
        return 1;
    }
    
    // Alineamiento
    bam1_t* aln = bam_init1();
    
    std::cout << "Procesando archivo BAM en streaming..." << std::endl;
    
    uint64_t total_reads = 0;
    uint64_t mapped_reads = 0;
    uint64_t paired_reads = 0;
    
    // Leer alineamientos uno por uno (streaming)
    while (sam_read1(bam_fp, header, aln) >= 0) {
        total_reads++;
        
        // Mostrar progreso cada 1M reads
        if (total_reads % 1000000 == 0) {
            std::cout << "Procesados: " << total_reads << " reads..." << std::endl;
        }
        
        // Calidad de mapeo (MAPQ)
        int mapq = aln->core.qual;
        mapping_quality_sketch.update(static_cast<float>(mapq));
        
        // Longitud del read
        int read_len = aln->core.l_qseq;
        read_length_sketch.update(static_cast<float>(read_len));
        
        // Verificar si está mapeado
        if (!(aln->core.flag & BAM_FUNMAP)) {
            mapped_reads++;
        }
        
        // Insert size (solo para paired-end properly paired)
        if ((aln->core.flag & BAM_FPAIRED) && 
            (aln->core.flag & BAM_FPROPER_PAIR) &&
            aln->core.isize > 0) {
            paired_reads++;
            insert_size_sketch.update(static_cast<float>(aln->core.isize));
        }
    }
    
    // Limpiar
    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(bam_fp);
    
    std::cout << "\n¡Procesamiento completado!" << std::endl;
    
    // --- Resultados ---
    std::cout << std::fixed << std::setprecision(2);
    
    std::cout << "\n=== ESTADÍSTICAS GENERALES ===" << std::endl;
    std::cout << "Total de reads: " << total_reads << std::endl;
    std::cout << "Reads mapeados: " << mapped_reads 
              << " (" << (100.0 * mapped_reads / total_reads) << "%)" << std::endl;
    std::cout << "Reads paired proper: " << paired_reads << std::endl;
    
    // Calidad de mapeo
    std::cout << "\n=== CALIDAD DE MAPEO (MAPQ) ===" << std::endl;
    std::cout << "Mediana: " << mapping_quality_sketch.get_quantile(0.5) << std::endl;
    std::cout << "P25: " << mapping_quality_sketch.get_quantile(0.25) << std::endl;
    std::cout << "P75: " << mapping_quality_sketch.get_quantile(0.75) << std::endl;
    std::cout << "P90: " << mapping_quality_sketch.get_quantile(0.90) << std::endl;
    std::cout << "Mínimo: " << mapping_quality_sketch.get_min_item() << std::endl;
    std::cout << "Máximo: " << mapping_quality_sketch.get_max_item() << std::endl;
    std::cout << "Compresión: " << (float)mapping_quality_sketch.get_n() 
              / mapping_quality_sketch.get_num_retained() << "x" << std::endl;
    
    // Longitud de reads
    std::cout << "\n=== LONGITUD DE READS ===" << std::endl;
    std::cout << "Mediana: " << read_length_sketch.get_quantile(0.5) << std::endl;
    std::cout << "P25: " << read_length_sketch.get_quantile(0.25) << std::endl;
    std::cout << "P75: " << read_length_sketch.get_quantile(0.75) << std::endl;
    std::cout << "Mínimo: " << read_length_sketch.get_min_item() << std::endl;
    std::cout << "Máximo: " << read_length_sketch.get_max_item() << std::endl;
    
    // Insert size (solo si hay paired-end data)
    if (paired_reads > 0) {
        std::cout << "\n=== INSERT SIZE (Paired-end) ===" << std::endl;
        std::cout << "Mediana: " << insert_size_sketch.get_quantile(0.5) << std::endl;
        std::cout << "P25: " << insert_size_sketch.get_quantile(0.25) << std::endl;
        std::cout << "P75: " << insert_size_sketch.get_quantile(0.75) << std::endl;
        std::cout << "P90: " << insert_size_sketch.get_quantile(0.90) << std::endl;
        std::cout << "P99: " << insert_size_sketch.get_quantile(0.99) << std::endl;
        std::cout << "Compresión: " << (float)insert_size_sketch.get_n() 
                  / insert_size_sketch.get_num_retained() << "x" << std::endl;
    }
    
    return 0;
}