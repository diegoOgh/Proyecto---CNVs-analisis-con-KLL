Detección de CNVs en WGS usando KLL Sketch

Proyecto semestral realizado para el electivo de postgrado
“Tópicos en Manejo de Grandes Volúmenes de Datos”.

Este proyecto evalúa el uso de KLL Sketch para estimar cuantiles de cobertura en datos de Whole Genome Sequencing (WGS) y aplicar estos cuantiles a la detección de Copy Number Variations (CNVs), reduciendo el uso de memoria y tiempo de cómputo frente a métodos exactos basados en ordenamiento.

Descripción general

La detección de CNVs en datos WGS requiere analizar la cobertura de millones de bins genómicos. El cálculo exacto de cuantiles mediante ordenamiento resulta costoso en tiempo y memoria.
En este proyecto se utiliza KLL Sketch, un algoritmo probabilístico de estimación de cuantiles que permite controlar explícitamente el error de rank mediante el parámetro K, manteniendo precisión suficiente para la detección de CNVs.

Estructura del repositorio
├── src/                    # Códigos fuente en C++
├── datasketches-cpp/       # Apache DataSketches (KLL)
├── graficos/               # Notebooks para generación de gráficos
├── *.csv                   # Outputs de los experimentos
└── README.md

Dependencias
Requisitos generales

Sistema operativo Linux o macOS

Compilador g++ con soporte C++17

HTSlib (lectura de archivos BAM)

HTSlib se utiliza para leer archivos BAM directamente desde C++.

Ubuntu / Debian

sudo apt update
sudo apt install -y libhts-dev


macOS (Homebrew)

brew install htslib

Samtools (preprocesamiento de archivos BAM)

Samtools se utiliza para:

Extraer cromosomas específicos desde un BAM completo

Indexar archivos BAM

Ubuntu / Debian

sudo apt update
sudo apt install -y samtools


macOS (Homebrew)

brew install samtools


Verificación:

samtools --version

Apache DataSketches (KLL – C++)

Clonar el repositorio dentro de la carpeta raíz del proyecto:

git clone https://github.com/apache/datasketches-cpp.git


El proyecto utiliza directamente los headers:

datasketches-cpp/kll/include

datasketches-cpp/common/include

No es necesario compilar la librería completa.

Datos utilizados

Este proyecto utiliza datos reales de Whole Genome Sequencing (WGS) del consorcio Genome in a Bottle (GIAB), correspondientes al individuo HG002 (NA24385) del trío Ashkenazim.

Fuente de los datos

Repositorio índice:

https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_Illumina300X_wgs_novoalign_GRCh37_GRCh38_NHGRI_07282015.HG002


Archivo BAM completo (GRCh38, ~300X):

ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam


Archivo índice:

ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam.bai


⚠️ El archivo BAM completo tiene un tamaño aproximado de 500 GB.

Preparación del archivo BAM usado en el proyecto

Para los experimentos no se utiliza el genoma completo, sino únicamente los cromosomas:

chr1, chr2, chr3, chr4, chr5


El archivo esperado por el proyecto es:

HG002.chr1-5.bam


A continuación se muestran dos formas alternativas de obtener este archivo.

Opción A: Descargar el BAM completo y extraer los cromosomas 1–5

⚠️ Esta opción requiere aproximadamente 500 GB de espacio en disco.

Descargar el BAM completo y su índice

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam.bai


Extraer los cromosomas 1–5

samtools view -b HG002.GRCh38.300x.bam chr1 chr2 chr3 chr4 chr5 > HG002.chr1-5.bam


Indexar el nuevo archivo BAM

samtools index HG002.chr1-5.bam

Opción B: Descargar solo los cromosomas 1–5 (recomendado)

Esta opción evita descargar el BAM completo y reduce drásticamente el espacio necesario.

Requiere que el servidor FTP soporte acceso aleatorio mediante el índice .bai (lo cual se cumple para los datos de GIAB).

Descargar únicamente el índice

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam.bai


Extraer directamente los cromosomas 1–5 desde el BAM remoto (190 GB)

samtools view -b \
  ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam \
  chr1 chr2 chr3 chr4 chr5 > HG002.chr1-5.bam


Indexar el BAM resultante

samtools index HG002.chr1-5.bam

Consideraciones

Opción A es más simple conceptualmente, pero costosa en almacenamiento.

Opción B es la forma recomendada para reproducir este proyecto.

Todos los programas del repositorio asumen que el archivo resultante se llama:

HG002.chr1-5.bam


Si se utiliza un nombre distinto, debe ajustarse manualmente al ejecutar los programas.

Compilación y ejecución

⚠️ Todos los comandos deben ejecutarse desde la raíz del proyecto.

1. cnv_kll_experimentacion.cpp

Genera bins de cobertura y realiza experimentación inicial con KLL.

Compilación

g++ -O3 -std=c++17 src/cnv_kll_experimentacion.cpp \
    -I./datasketches-cpp/kll/include \
    -I./datasketches-cpp/common/include \
    -lhts \
    -o cnv_kll_experimentacion


Ejecución

./cnv_kll_experimentacion HG002.chr1-5.bam


Output

bin_experiment.csv

2. sort_vs_kll.cpp

Compara KLL Sketch con el método exacto basado en ordenamiento.

Compilación

g++ -O3 -std=c++17 src/sort_vs_kll.cpp \
    -I./datasketches-cpp/kll/include \
    -I./datasketches-cpp/common/include \
    -lhts \
    -o sort_vs_kll


Ejecución

./sort_vs_kll


Output

kll_vs_sort_comparison.csv

3. k_experimentacion.cpp

Evalúa el impacto del parámetro K en error, tiempo y memoria.

Compilación

g++ -O3 -std=c++17 src/k_experimentacion.cpp \
    -I./datasketches-cpp/kll/include \
    -I./datasketches-cpp/common/include \
    -lhts \
    -o k_experimentacion


Ejecución

./k_experimentacion HG002.chr1-5.bam 1000 k_experimentacion.csv

4. bam_reader_mejorado.cpp (baseline exacto)

Genera el baseline exacto de cobertura por bin.

Compilación

g++ -O3 -std=c++17 src/bam_reader_mejorado.cpp \
    -I./datasketches-cpp/kll/include \
    -I./datasketches-cpp/common/include \
    -lhts \
    -o kll_bam_reader


Ejecución

./kll_bam_reader HG002.chr1-5.bam 1000 cnv_1000.csv

5. cnv_pasada.cpp (detección de CNVs)

Detecta CNVs usando cuantiles de cobertura.

Compilación

g++ -O3 -std=c++17 src/cnv_pasada.cpp \
    -I./datasketches-cpp/kll/include \
    -I./datasketches-cpp/common/include \
    -lhts \
    -o cnv_pasada


Ejecución

./cnv_pasada HG002.chr1-5.bam 1000 cnv_1000.csv cnv_detection.csv 5

Gráficos

La carpeta graficos/ contiene notebooks de Jupyter para generar los gráficos del análisis.
Los nombres de los archivos .csv deben coincidir con los definidos en el notebook.

Conclusión

KLL Sketch permite estimar cuantiles de cobertura para la detección de CNVs en datos WGS con bajo uso de memoria y tiempo de cómputo, manteniendo errores de rank pequeños y controlables mediante el parámetro K, lo que lo convierte en una alternativa escalable frente a métodos exactos.