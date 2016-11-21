Import BRCA dataset from TCGA
================
André Veríssimo *and* Marta Lopes
November 2016

-   [Package information](#package-information)
    -   [How to use the dataset](#how-to-use-the-dataset)
        -   [Example](#example)
    -   [Description of the data set](#description-of-the-data-set)
        -   [Explaining the TCGA codes](#explaining-the-tcga-codes)
-   [Script to generate the package data](#script-to-generate-the-package-data)
    -   [Loading data from TCGA](#loading-data-from-tcga)
        -   [Importing using TCGA2STAT package](#importing-using-tcga2stat-package)
        -   [Extracting all the barcodes](#extracting-all-the-barcodes)
    -   [Processing the data](#processing-the-data)
        -   [Mapping all types of sample (tumor, normal, metastases, control...)](#mapping-all-types-of-sample-tumor-normal-metastases-control...)
        -   [Size of data from each tissue](#size-of-data-from-each-tissue)
        -   [Store all patient's barcode from `tissue` in `tissue.barcode`](#store-all-patients-barcode-from-tissue-in-tissue.barcode)
        -   [Store patient's clinical data](#store-patients-clinical-data)
    -   [Exported data](#exported-data)

Package information
===================

How to use the dataset
----------------------

1.  Install brca.data by using `devtools` package.

2.  Load the library

3.  Load the required datasets (one or more of the following)
    -   `clinical`
    -   `tissue`
    -   `tissue.all`
    -   `tissue.barcode`
    -   `tissue.ix`

### Example

``` r
# load library or use directly if insta
install.packages('devtools')
# The library can also be loaded and use the function install_git without 'devtools::' prefix
devtools::install_git('http://sels.tecnico.ulisboa.pt/gitlab/averissimo/rpackage-brca.git')
#
# Load the brca.data package
library(brca.data)
# start using the data, for example the tissue data
data(tissue)
# tissue is now in the enviromnet and will be loaded on the first
#  time it is used. For example:
names(tissue)
```

Description of the data set
---------------------------

The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) data collection is part of a larger effort to build a research community focused on connecting cancer phenotypes to genotypes by providing clinical images matched to subjects from The Cancer Genome Atlas ([TCGA](https://cancergenome.nih.gov/)).

The BRCA data is publicly available (<https://gdc-portal.nci.nih.gov/>) and is decomposed into two data sets:

1.  the gene expression data, composed of `20501` variables for a total of `1205` samples with `1093` individuals. From those samples, `1093` with primary solid tumor and `112` with normal tissue;

2.  the clinical data is composed of 18 variables obtained from the same individuals.

### Explaining the TCGA codes

The following links explain (1) the individuals' barcode and (2) the sample type code:

1.  [Link to individuals' barcode from tcga](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode?desktop=true&macroName=unmigrated-inline-wiki-markup)

2.  [Link to sample type code from tcga](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)

Explanation of TCGA barcode using example `TCGA-02-0001-01C-01D-0182-01`

![Details on TCGA Barcode](https://wiki.nci.nih.gov/download/attachments/39294833/barcode.png?version=1&modificationDate=1300400318000&api=v2)

<table>
<colgroup>
<col width="17%" />
<col width="25%" />
<col width="9%" />
<col width="24%" />
<col width="22%" />
</colgroup>
<thead>
<tr class="header">
<th>Label</th>
<th>Identifier for</th>
<th>Value</th>
<th>Value Description</th>
<th>Possible Values</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Project</td>
<td>Project name</td>
<td>TCGA</td>
<td>TCGA project</td>
<td>TCGA</td>
</tr>
<tr class="even">
<td>TSS</td>
<td>Tissue source site</td>
<td>02</td>
<td>GBM (brain tumor) sample from MD Anderson</td>
<td>See <a href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes">Code Tables Report</a></td>
</tr>
<tr class="odd">
<td>Participant</td>
<td>Study participant</td>
<td>0001</td>
<td>The first participant from MD Anderson for GBM study</td>
<td>Any alpha-numeric value</td>
</tr>
<tr class="even">
<td>Sample</td>
<td>Sample type</td>
<td>01</td>
<td>A solid tumor</td>
<td>Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See <a href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes">Code Tables Report</a> for a complete list of sample codes</td>
</tr>
<tr class="odd">
<td>Vial</td>
<td>Order of sample in a sequence of samples</td>
<td>C</td>
<td>The third vial</td>
<td>A to Z</td>
</tr>
<tr class="even">
<td>Portion</td>
<td>Order of portion in a sequence of 100 - 120 mg sample portions</td>
<td>01</td>
<td>The first portion of the sample</td>
<td>01-99</td>
</tr>
<tr class="odd">
<td>Analyte</td>
<td>Molecular type of analyte for analysis</td>
<td>D</td>
<td>The analyte is a DNA sample</td>
<td>See <a href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes">Code Tables Report</a></td>
</tr>
<tr class="even">
<td>Plate</td>
<td>Order of plate in a sequence of 96-well plates</td>
<td>0182</td>
<td>The 182nd plate</td>
<td>4-digit alphanumeric value</td>
</tr>
<tr class="odd">
<td>Center</td>
<td>Sequencing or characterization center that will receive the aliquot for analysis</td>
<td>01</td>
<td>The Broad Institute GCC</td>
<td>See <a href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes">Code Tables Report</a></td>
</tr>
</tbody>
</table>

Script to generate the package data
===================================

``` r
#
#
#
#
# This now documents the generation of the dataset.
#
# No longer a generic documentation of the package!
#
# 
#
```

Loading data from TCGA
----------------------

### Importing using TCGA2STAT package

The BRCA data set can also be extracted through the [TCGA2STAT R package](https://cran.r-project.org/web/packages/TCGA2STAT/index.html)

``` r
# Library that downloads the data
library(TCGA2STAT)
# Library to print some information
library(futile.logger)
# Cached filed is used internally for the generation of this document. Not included in package.
#  As only the post-processed variables are of interest.
if (file.exists('cache.RData')) {
  load('cache.RData')
} else {
  brca <- TCGA2STAT::getTCGA(disease = "BRCA", data.type = "RNASeq2", type = "RPKM", clinical = TRUE)
  save(brca, file = 'cache.RData')
}
```

### Extracting all the barcodes

Pattern to retrieve TCGA patient and sample type from extended barcode.

``` r
# $1 gets the TCGA participant portion
# $2 gets the sample type
tcga.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4})-([0-9]{2}).+'
#
# patient barcode for each column in brca$dat
tissue.barcode <- list()
tissue.barcode$all <- gsub(tcga.barcode.pattern, '\\1', colnames(brca$dat))
#
# getting all tisues
tissue <- list()
tissue$all <- brca$dat
#
# clinical data
clinical <- list()
clinical$all <- brca$clinical
#
#
getParticipant <- function(fullBarcode) {
  source('R/functions.R')
  getParticipantCode(fullBarcode)
}
```

Processing the data
-------------------

### Mapping all types of sample (tumor, normal, metastases, control...)

See [Code Tables Report](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) for a complete list of sample codes.

``` r
# sample type per sample
tissue.type  <- as.numeric(gsub(tcga.barcode.pattern,'\\2', colnames(brca$dat)))
#
# mapping for tissue type
tissue.mapping <- list( 
  'primary.solid.tumor' = 01,
  'recurrent.solid.tumor' = 02,
  'primary.blood.derived.cancer.peripheral.blood' = 03,
  'recurrent.blood.derived.cancer.bone.marrow' = 04,
  'additional.new.primary' = 05,
  'metastatic' = 06,
  'additional.metastatic' = 07,
  'human.tumor.original.cells' = 08,
  'primary.blood.derived.cancer.bone.marrow' = 09,
  'blood.derived.normal' = 10,
  'solid.tissue.normal' = 11,
  'buccal.cell.normal' = 12,
  'ebv.immortalized.normal' = 13,
  'bone.marrow.normal' = 14,
  'sample.type.15' = 15,
  'sample.type.16' = 16,
  'control.analyte' = 20,
  'recurrent.blood.derived.cancer.peripheral.blood' = 40,
  'cell.lines' = 50,
  'primary.xenograft.tissue' = 60,
  'cell.line.derived.xenograft.tissue' = 61,
  'sample.type.99' = 99
)
#
# Types of tissues
tissue.ix <- list()
for (el in names(tissue.mapping)) {
  tissue.ix[[el]] <- (tissue.type == tissue.mapping[[el]])
  if (any(tissue.ix[[el]]))
    tissue[[el]] <- tissue$all[,tissue.ix[[el]]]
}
```

### Size of data from each tissue

``` r
sample.size <- c()
for (el in names(tissue)) {
  sample.size <- rbind(sample.size, c(ncol(tissue[[el]]), nrow(tissue[[el]])))
}
rownames(sample.size) <- names(tissue)
colnames(sample.size) <- c('# of Samples', '# of Genes')
futile.logger::flog.info('Tissue information per tissue type:', sample.size, capture = TRUE)
```

    ## INFO [2016-11-21 19:42:24] Tissue information per tissue type:
    ## 
    ##                     # of Samples # of Genes
    ## all                         1212      20501
    ## primary.solid.tumor         1093      20501
    ## metastatic                     7      20501
    ## solid.tissue.normal          112      20501

### Store all patient's barcode from `tissue` in `tissue.barcode`

``` r
for (el in names(tissue)) {
  tissue.barcode[[el]] <- getParticipant(colnames(tissue[[el]]))
}
```

### Store patient's clinical data

``` r
for (el in names(tissue)) {
  clinical[[el]] <- clinical$all[tissue.barcode[[el]],]
}
sample.size <- c()
for (el in names(clinical)) {
  sample.size <- rbind(sample.size, c(nrow(clinical[[el]]), ncol(clinical[[el]])))
}
rownames(sample.size) <- names(tissue)
colnames(sample.size) <- c('# of Samples', '# of Features')
futile.logger::flog.info('Clinical information per tissue type:', sample.size, capture = TRUE)
```

    ## INFO [2016-11-21 19:42:24] Clinical information per tissue type:
    ## 
    ##                     # of Samples # of Features
    ## all                         1212            18
    ## primary.solid.tumor         1093            18
    ## metastatic                     7            18
    ## solid.tissue.normal          112            18

Exported data
-------------

-   `tissue.ix`: indexes of all tissues;
    -   Has all types of possible tissue samples;
-   `tissue`: gene expression data for all types of tissues, see `names(tissue)`;
-   `tissue.barcode`: Patient's participation data code (TCGA-XX-XXXX) per type of tissue;
-   `clinical`: Clinical data per tissue type. Has the same structure as tissue.

``` r
# Extract all expression data to a different variable to remove redundancy in `tissue`
tissue.all     <- tissue$all
tissue$all     <- 'Stored in tissue.all in brca.data package'
#
#
# Uncomment to update data variables
#
# devtools::use_data(tissue.ix, tissue, tissue.barcode, clinical, tissue.all, overwrite = TRUE)
```
