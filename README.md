# Distancia de Hausdorff Paralela "Naive", basado en HDKP1, implementado con SYCL

Implementación paralela "embarazosa" del algoritmo **HDKP1** para el cálculo de distancia de Hausdorff (dirigida y simétrica) sobre conjuntos de puntos 2D representados mediante $k^2$-trees (`MREP2`).

---

## Concepto de la Implementación 

Esta versión corresponde al reemplazo directo del ciclo secuencial `for (p in A)` por una ejecución paralela:

1. **Extracción de Puntos**: Se extraen los puntos de $K_A$ a un arreglo/vector plano mediante `extractPointK2tree(A, 0)`.
2. **Paralelización del Ciclo**: Se divide el rango de puntos $[0, N)$ estáticamente entre los workers (`std::async`).
3. **`NNMAX` en distintos workers**: Se invoca directamente la función secuencial original `nnMax` (`src/seq/NNMAX.cpp`)
5. **Cota Global Atómica `cmax`**: Variable en memoria compartida (USM de SYCL) actualizada asincrónicamente mediante `sycl::atomic_ref<double, sycl::memory_order::relaxed>::fetch_max`.

---

## Compilación y Ejecución

### Requisitos
- **Compilador**: Intel oneAPI DPC++/C++ Compiler (`icpx`).
- **Target**: CPU (`-fsycl -fsycl-targets=spir64_x86_64 -O3`).
- **CMake**: >= 3.20.

Cargar el entorno Intel oneAPI:
```bash
source /opt/intel/oneapi/setvars.sh
```

### Compilar
```bash
cd /path/to/HD_K2_Par
mkdir -p build && cd build
cmake .. -DCMAKE_CXX_COMPILER=icpx -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)
```

### Ejecutar
```bash
./build/Parallel_HSD_K2 data/set1Conjunto7R.csv.kt data/set2Conjunto7R.csv.kt
```

### Argumentos de línea de comandos:
```bash
./build/Parallel_HSD_K2 [datasetA.csv.kt] [datasetB.csv.kt]
```
Si se omiten los argumentos, se ejecutan por defecto `set1Conjunto1T.csv.kt` y `set2Conjunto1T.csv.kt`.

---

## Estructura del Repositorio

```text
HD_K2_Par/
├── CMakeLists.txt          # Configuración CMake con soporte icpx y SYCL
├── README.md               # Documentación e instrucciones
├── src/
│   ├── main.cpp            # Driver de benchmarking y verificación (Seq vs Par)
│   ├── ParallelHDKP1.h     # Declaraciones de la implementación paralela ingenua
│   ├── ParallelHDKP1.cpp   # Implementación paralela ingenua sobre puntos
│   └── seq/                # Algoritmos secuenciales base (HDKP1, HDKP2, NNMAX original)
├── K2tree/                 # Estructura compacta k2-tree
├── Util/                   # Primitivas geométricas (Point, SpatialArea)
├── libcds/                 # Biblioteca de estructuras compactas
└── data/                   # Datasets reales de prueba (.csv.kt)
```
