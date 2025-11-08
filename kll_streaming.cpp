#include <iostream>
#include <random>
#include <iomanip>
#include "datasketches-cpp/kll/include/kll_sketch.hpp"

using namespace datasketches;

int main() {
    kll_sketch<float> kll_sketch(200);
    
    std::cout << "Procesando 10,000,000 de puntos en streaming..." << std::endl;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float> dist(0.0f, 1.0f);
    
    // Streaming puro: NO guardamos los datos
    for (int i = 0; i < 10'000'000; ++i) {
        kll_sketch.update(dist(gen));
    }
    
    std::cout << "¡Sketch actualizado!" << std::endl;
    
    // Consultar resultados
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n--- Resultados del KLL Sketch ---" << std::endl;
    std::cout << "Mediana (p50): " << kll_sketch.get_quantile(0.5) << std::endl;
    std::cout << "Percentil 90: " << kll_sketch.get_quantile(0.90) << std::endl;
    std::cout << "Percentil 99: " << kll_sketch.get_quantile(0.99) << std::endl;
    std::cout << "Items procesados: " << kll_sketch.get_n() << std::endl;
    std::cout << "Items retenidos: " << kll_sketch.get_num_retained() << std::endl;
    std::cout << "Compresión: " << (float)kll_sketch.get_n() / kll_sketch.get_num_retained() << "x" << std::endl;
    
    return 0;
}