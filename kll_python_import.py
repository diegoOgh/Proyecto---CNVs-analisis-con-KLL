import numpy as np
from datasketches import kll_floats_sketch as KLL

kll_sketch = KLL(k=200)

print("Actualizando el sketch con 10,000,000 de puntos...")
data_stream = np.random.randn(10_000_000)

total_items = len(data_stream)

kll_sketch.update(data_stream)
print("¡Sketch actualizado!")

# --- 3. Output: Consultar el Sketch ---
median_approx = kll_sketch.get_quantile(0.5)
p90_approx = kll_sketch.get_quantile(0.90)
p99_approx = kll_sketch.get_quantile(0.99)

# CORRECCIÓN FINAL: Mezcla de propiedades y métodos
num_items_retenidos = kll_sketch.num_retained     # Propiedad (sin paréntesis)
min_val_aprox = kll_sketch.get_min_value()        # Método (con paréntesis)
max_val_aprox = kll_sketch.get_max_value()        # Método (con paréntesis)

# --- 4. Resultados ---
print("\n--- Resultados del KLL Sketch ---")
print(f"Mediana (p50) aproximada: {median_approx:.4f}")
print(f"Percentil 90 (p90) aproximado: {p90_approx:.4f}")
print(f"Percentil 99 (p99) aproximado: {p99_approx:.4f}")
print(f"Total de items procesados: {total_items}")
print(f"Items retenidos en memoria: {num_items_retenidos}")
print(f"Compresión: {total_items/num_items_retenidos:.1f}x")
print(f"Valor mínimo: {min_val_aprox:.4f}")
print(f"Valor máximo: {max_val_aprox:.4f}")

# --- Comparación ---
print("\n--- Comparación con valores reales (NumPy) ---")
print(f"Mediana real: {np.quantile(data_stream, 0.5):.4f}")
print(f"Percentil 90 real: {np.quantile(data_stream, 0.90):.4f}")
print(f"Percentil 99 real: {np.quantile(data_stream, 0.99):.4f}")