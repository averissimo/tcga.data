TCGA.DATA R Package
================

-   [Package information](#package-information)
    -   [How to use the dataset](#how-to-use-the-dataset)
-   [How to build own data package](#how-to-build-own-data-package)
-   [Ackowledgements](#ackowledgements)

This R Package allows to retrieve Gene Expression, Mutation and clinical data from [TCGA database](http://gdc-portal.nci.nih.gov/) (The Cancer Genome Atlas). It retrieves a single type of cancer at a time.

We publish diferent package in the [releases page](https://github.com/averissimo/tcga.data/releases) that allow to quickly use the datasets.

The genome expression datasets are already in a matrix format ready to be used. The data is in FPKM (Fragments Per Kilobase Million) format. Any additional normalization to use in models must be performed

Package information
-------------------

### How to use the dataset

1.  Install `brca.data` by using `devtools` package. (`brca.data`, `prad.data` or `skcm.data`)

2.  Load the library

3.  Load the required datasets (one or more of the following)
    -   `multiAssay`
    -   `gdc.original`

In older versions of this package, prior to September 2018, the dataset was named `fpkm.per.tissue` or `mutation`, but we since improved the storage using a MultiAssayExperiment object from bioconductor.

To recover the datasets in the old matrix format use the following

``` r
data('multiAssay')
fpkm.data <- build.matrix('RNASeqFPKM', multiAssay)
fpkm.per.tissue <- fpkm.data$data
fpkm.clinical   <- fpkm.data$clinical
```

#### Example for BRCA package

``` r
# The library can also be loaded and use the function install_git without 'devtools::' prefix
BiocManager::install('https://github.com/averissimo/tcga.data/releases/download/2016.12.15-brca/brca.data_1.0.tar.gz')
#
# Load the brca.data package
library(brca.data)
# start using the data, for example the tissue data
data(fpkm.per.tissue)
# tissue is now in the enviromnet and will be loaded on the first
#  time it is used. For example:
names(fpkm.per.tissue)
```

How to build own data package
-----------------------------

1.  Open vignettes/build\_data.Rmd
2.  Change in the header of the Rmd *(beginning of the document)* the project param to the target TCGA project
3.  Open DESCRITION and change the name of the package to the desired name

-   we use a convention of \#\#\#\#.data where \#\#\#\# is the tcga project name in lowercase

1.  Run the vignettes/build\_data.Rmd to build the cache of the data
2.  Run `devtools::document()` to create documentation
3.  Run `devtools::build()` to build the actual package

Ackowledgements
---------------

This package was developed primarily by *[André Veríssimo](http://web.tecnico.ulisboa.pt/andre.verissimo/)* with support from *Marta Lopes*, *Eunice Carrasquinha* and *[Susana Vinga](http://web.tecnico.ulisboa.pt/susanavinga/)*

This work was supported by:

-   [FCT](www.fct.pt), through IDMEC, under LAETA, projects *(UID/EMS/50022/2013)*;
-   Susana Vinga acknowledges support by program Investigador FCT *(IF/00653/2012)* from [FCT](www.fct.pt), co-funded by the European Social Fund *(ESF)* through the Operational Program Human Potential *(POPH)*;
-   André Veríssimo acknowledges support from [FCT](www.fct.pt) *(SFRH/BD/97415/2013)*.
