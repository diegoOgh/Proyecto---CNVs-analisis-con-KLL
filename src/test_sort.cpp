#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

// --- Implementación manual de QuickSort ---
void quicksort(std::vector<int>& arr, int left, int right) {
    if (left >= right) return;

    int pivot = arr[(left + right) / 2];
    int i = left;
    int j = right;

    while (i <= j) {
        while (arr[i] < pivot) i++;
        while (arr[j] > pivot) j--;
        if (i <= j) {
            std::swap(arr[i], arr[j]);
            i++;
            j--;
        }
    }

    if (left < j) quicksort(arr, left, j);
    if (i < right) quicksort(arr, i, right);
}

int main() {
    const size_t N = 200'000'000;

    // Generar datos aleatorios
    std::vector<int> datos(N);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, 1'000'000'000);
    for (auto &x : datos)
        x = dist(rng);

    // --- Copia para std::sort ---
    std::vector<int> datos_std = datos;

    // --- Quicksort manual ---
    auto inicio1 = std::chrono::high_resolution_clock::now();
    quicksort(datos, 0, datos.size() - 1);
    auto fin1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracion1 = fin1 - inicio1;
    std::cout << "Quicksort manual: " << duracion1.count() << " segundos\n";

    // --- std::sort para comparar ---
    auto inicio2 = std::chrono::high_resolution_clock::now();
    std::sort(datos_std.begin(), datos_std.end());
    auto fin2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracion2 = fin2 - inicio2;
    std::cout << "std::sort: " << duracion2.count() << " segundos\n";

    return 0;
}
