# Singel-cell-RNAseq-Analysis
The pipeline aims to streamline processing raw fastq files from 10X platfrom from accounting for batch effects, sequence alignment, clustering, cell type annotation, and downstream differential analysis, generating publication-ready plots. Most of the analyses are built on top of existing tools and/or packages publicly available. Please find the reference as indicated. Here, the pipeline demonstration consisits of two parts, QC and batch effect adjustment in Bash script (part I) and clustering annotation in R script (Part II):


Part I. Cell Ranger (Bash)
-----
- [Cell Ranger](https://github.com/10XGenomics/cellranger): Process of raw bcl data from sequencer including demultiplexing, reads alignment, and normalization of library size

Part II. Clustering and differnetial analysis (R)
-----
- [SingleCellExperiment](https://github.com/drisso/SingleCellExperiment): Data storage tool for data types required for single cell data analysis
- [Seurat](https://github.com/satijalab/seurat): Toolkit for single cell genomics analysis and figure plotting
- [scater](https://github.com/jimhester/scater): Single cell data pre-processing before downstream analysis
- [scuttle](https://rdrr.io/github/LTLA/scuttle/): Sequencing QC, reads normalization, and data transformation
- [scran](https://rdrr.io/bioc/scran/): For cell phase assignment and feature selection of highly variable genes 
- [celldex](https://rdrr.io/github/LTLA/celldex/): Import hematopoietic cell type assignment refernce profile from ImmGen
- [SingleR](https://github.com/dviraran/SingleR): Cell type annotation of single cell based on expression profile of the reference data
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html): Differential expression analysis of gene expression profile with biological replicates 
- [ggplot2](https://github.com/tidyverse/ggplot2): For basic QC plotting
- [dplyr](https://github.com/tidyverse/dplyr): Organize count data for downstream analysis
