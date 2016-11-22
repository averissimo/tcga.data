#' brca.data: Breast Cancer data from TCGA
#'
#' To use this package load the following datasets:
#'
#' - `tissue.ix`: indexes of all tissues;
#' - `tissue`: gene expression data for all types of tissues, see `names(tissue)`;
#' - `tissue.barcode`: Patient's participation data code (TCGA-XX-XXXX) per type of tissue;
#' - `clinical`: Clinical data per tissue type. Has the same structure as tissue;
#' - `gene.ranges`: Genomic ranges with description for genes in tissue[.all] matrices.
#'
#' To see how the data was retrieved, check the `README.Rmd` in the package source.
#'
#' @name brca.data
#' @docType package
NULL
