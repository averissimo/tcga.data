# Import datasets from TCGA

André Veríssimo *and* Marta Lopes

updated July 2017

-   [Package information](#package-information)
    -   [How to use the dataset](#how-to-use-the-dataset)
        -   [Example](#example)
    -   [Description of the data set](#description-of-the-data-set)
        -   [Explaining the TCGA codes](#explaining-the-tcga-codes)
-   [Acknowledgements](#acknowledgements)
    -    [Funding](#funding)


# Package information

The Cancer Genome Atlas data collection is part of a larger effort to build a research community focused on connecting cancer phenotypes to genotypes by providing clinical images matched to subjects from The Cancer Genome Atlas ([TCGA](https://cancergenome.nih.gov/)).

The data is publicly available (<https://gdc-portal.nci.nih.gov/>).

## How to use the dataset

1.  Download a released package (brca.data or tcga.data)

2.  Install it using `install.packages('<path to package>', repos = NULL, type="source")`
    -   Alternatively use the link directly with devtools package `devtools::install_url('<link to package zip>')`

3.  Load the library

4.  Load the required datasets (one or more of the following)
    -   `clinical`
    -   `fpkm.per.tissue`
    -   `fpkm.per.tissue.barcode`
    -   `mutation`
    -   `gdc`

### Example using the direct link

``` r
# load library or use directly if insta
install.packages('devtools')
# The library can also be loaded and use the function install_git without 'devtools::' prefix
devtools::install_url('https://github.com/averissimo/tcga.data/releases/download/2016.12.15-brca/brca.data_1.0.tar.gz')
#
# Load the brca.data package
library(brca.data)
# start using the data, for example the tissue data
data(fpkm.per.tissue)
# tissue is now in the enviromnet and will be loaded on the first
#  time it is used. For example:
names(fpkm.per.tissue)
```

## Description of the data set

The Cancer Genome Atlas data collection is part of a larger effort to build a research community focused on connecting cancer phenotypes to genotypes by providing clinical images matched to subjects from The Cancer Genome Atlas ([TCGA](https://cancergenome.nih.gov/)).

The data is publicly available (<https://gdc-portal.nci.nih.gov/>) and this package only has available the following data:

1.  the gene expression data, decomposed on the origin of the sample (i.e. normal tissue, primary solid tumor, etc..);

2.  the clinical data also decomposed in the same way as gene expression data. (additional information is also available in the gdc data variable, such as follow\_up, drug, radiation, ...);

3.  mutation information presented in a sparse matrix with the count of mutations per patient/ensemble gene.

### Explaining the TCGA codes

The following links explain (1) the individuals' barcode and (2) the sample type code:

1.  [Link to individuals' barcode from tcga](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode?desktop=true&macroName=unmigrated-inline-wiki-markup)

1. [Link to sample type code from tcga](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)

Explanation of TCGA barcode using example `TCGA-02-0001-01C-01D-0182-01`

![Details on TCGA Barcode](https://wiki.nci.nih.gov/download/attachments/39294833/barcode.png?version=1&modificationDate=1300400318000&api=v2)

| Label       | Identifier for     | Value | Value Description | Possible Values |
|-------------|--------------------|-------|-------------------|-----------------|
| Project     | Project name       | TCGA | TCGA project | TCGA |
| TSS         | Tissue source site | 02 | GBM (brain tumor) sample from MD Anderson | See [Code Tables Report](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes) |
| Participant | Study participant | 0001 | The first participant from MD Anderson for GBM study | Any alpha-numeric value |
| Sample      | Sample type       | 01 | A solid tumor | Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See [Code Tables Report](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) for a complete list of sample codes |
| Vial        | Order of sample in a sequence of samples       | C | The third vial | A to Z |
| Portion     | Order of portion in a sequence of 100 - 120 mg sample portions | 01 | The first portion of the sample | 01-99 |
| Analyte     | Molecular type of analyte for analysis         | D | The analyte is a DNA sample | See [Code Tables Report](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes) |
| Plate       | Order of plate in a sequence of 96-well plates | 0182 | The 182nd plate | 4-digit alphanumeric value |
| Center      | Sequencing or characterization center that will receive the aliquot for analysis | 01 | The Broad Institute GCC | See [Code Tables Report](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes) |


# Acknowledgements

Marta Lopes was a fundamental person on the development of this package framework and methods, providing feedback forcing this to be more than just a personal package.

Susana Vinga for supervising and guidance of this work.

## Funding

### Projects

-    EU H2020 Personalizing Health and Care Program - [https://ec.europa.eu/commission/](https://ec.europa.eu/commission/index_en)
    -    **SOUND** – Statistical multi-Omics UNDerstanding of Patient Samples – Horizon 2020 - 633974 – [http://www.sound-biomed.eu/](http://www.sound-biomed.eu/)

-    Fundação para a Ciência e Tecnologia (FCT) - [http://www.fct.pt/](http://www.fct.pt/)
    -    **PERSEIDS** – Personalizing cancer therapy through integrated modeling and decision – FCT Project PTDC/EMS-SIS/0642/2014

### Fellowships

-    **Fundação para a Ciência e Tecnologia (FCT)** - [http://www.fct.pt/](http://www.fct.pt/)
    -    André Veríssimo's PhD grant - SFRH/BD/97415/2013
    -    Susana Vinga' Investigador FCT grant - IF/00653/2012
