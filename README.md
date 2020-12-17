# mihola2021
Analyses for Mihola et al paper (rat DSB hotspots) 2021

PREREQUISITES: 
The pipeline is built in nextflow (nextflow.io) using singularity. It has been tested using nextflow/20.10.0 and singularity/3.6.4. Earlier versions of nextflow will not work.

ACCESSORY DATA: 
Aligned BAM files for all samples are needed to run the pipeline. These can be obtained from:
https://hpc.nih.gov/~brickkm/mihola2021/mihola2021.tar.gz

Extract this tar in the folder that contains the mihola.nf script

RUN: 
Clone the repo. From the source folder, run:

nextflow run -c config.nf -profile local mihola2021.nf
