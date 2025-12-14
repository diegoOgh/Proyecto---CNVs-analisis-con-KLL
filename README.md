# Proyecto---CNVs-analisis-con-KLL
Proyecto semestral realizado para el electivo de postgrado "TOPICOS EN MANEJO DE GRANDES VOLUMENES DE DATOS"


//Provisorio:

Para correr bam_reader_mejorado.cpp deben hacer git clone del KLL de apache dentro de la carpeta del proyecto, luego compilar de la siguiente manera:

g++ -O3 -std=c++17 bam_reader_mejorado.cpp     -I./datasketches-cpp/kll/include     -I./datasketches-cpp/common/include     -o cnv_kll     -lhts

Y luego correrlo con :

./cnv_kll HG002.chr1-5.bam 1000