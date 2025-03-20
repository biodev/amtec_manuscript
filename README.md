## AMTEC Manuscript Workflow

Workflow code supporting the manuscript 'Utilizing cohort-level and individual networks to predict best response in patients with metastatic triple negative breast cancer'

## Required files

Retrieve the below files and place them in the corresponding directory

-   `data/annotation`
    -   [gencode.v24.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz)
    -   [h.all.v7.5.1.symbols.gmt](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt)
    -   [Cell_marker_Human.xlsx](http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt)
    -   `10780432ccr173509-sup-192911_3_supp_4675335_p6rxmz.xlsx`
        -   This is Table S1 from [Tamborero et al 2018](https://aacrjournals.org/clincancerres/article/24/15/3717/80876/A-Pan-cancer-Landscape-of-Interactions-between)
-   `data/burstein`
    -   `TNBC_Ding_77_gene_signatures.txt` derived from [Table S1](https://pmc.ncbi.nlm.nih.gov/articles/instance/6349443/bin/oncotarget-10-198-s002.docx)
-   `data/tcga`
    -   [TCGA_BRCA_counts.tsv.gz](https://osf.io/gqrz9/files/osfstorage/57855f67594d9001fbb9c53e)
    -   [TCGA_BRCA_tpm.tsv.gz](https://osf.io/gqrz9/files/osfstorage/57855fdd9ad5a10200c739be)
    -   [TCGA_ID_MAP.csv](https://osf.io/gqrz9/files/osfstorage/578e44e29ad5a101fe95330e)
    -   [brca_tcga_pan_can_atlas_2018_clinical_data.tsv](https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018)
-   `data/zhang`
    -   bcell_wrkbk.xlsx derived from [Zhang et al 2021](https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00499-2?elqTrackId=b1cd684119f04471ab546f6b63bdf5e8)
        -   Copy the gene list embedded in the 'Survival analysis' section of the 'Method details' into an Excel file with two columns: **signature** and **genes**
-   `data/metabric`
    -   Download the following files from <https://www.cbioportal.org/study/summary?id=brca_metabric>:
        -   `data_clinical_patient.txt`
        -   `data_clinical_sample.txt`
        -   `data_mrna_agilent_microarray.txt`
-   `data/amtec` and `data/wgcna`
        - [Request](mailto:mcweeney@ohsu.edu) access to these files

## Installation

``` r
#CRAN:
install.packages(c("data.table", "targets", "tidymodels", 
"openxlsx", "ggplot2", "stringr", "ggrepel", "WGCNA",
"caret", "partykit", "patchwork", "viridis", "ggsurvfit"))

#BioConductor
##May need to install "BiocManager" first from CRAN
BiocManager::install(c("limma", "edgeR", "sva", "GSVA", "ComplexHeatmap"))
```

## Pipeline Execution

``` bash
mkdir output 
mkdir figures
```

``` r
targets::tar_make()
```

## Licensing

[GNU General Public License v3.0](LICENSE)

This code was developed by Daniel Bottomly, a member of the McWeeney Lab and is protected under copyright by the Oregon Health and Science University Knight Cancer Institute, 2024.
