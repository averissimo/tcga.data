#' Gene expression per tissue.
#'
#' A dataset containing a list of tissue names from TCGA breast cancer dataset (BRCA).
#'
#' Each element has a normalized gene expression (RPKM) matrix with genes (rows) x patients (columns).
#'
#' For all gene expression data matrix, see `tissue.all` dataset
#'
#' @format A list with RPKM gene expression data.
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"tissue"

#' Gene expression per tissue.
#'
#' A dataset containing a list of tissue names from TCGA breast cancer dataset (BRCA).
#'
#' Each element has a normalized gene expression (RPKM) matrix with genes (rows) x patients (columns).
#'
#' It also has an element with all the dataset.
#'
#' @format A list with RPKM gene expression data.
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"tissue.all"

#' Index of patients per tissue.
#'
#' A dataset containing a list of index tissue names from TCGA breast cancer dataset (BRCA).
#'
#' Each element has a vector of logical index for the tissue data.
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"tissue.ix"

#' Patient's barcode per tissue type
#'
#' A dataset containing a list of vectors with TCGA barcodes per tissue type.
#'
#' @format A list with RPKM gene expression data.
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"tissue.barcode"

#' Clinical data per tissue type
#'
#' A dataset containing a list of dataframes with clinical data per tissue type.
#'
#' @format A list with RPKM gene expression data.
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"clinical"
