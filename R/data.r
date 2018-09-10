#' Original data from GDC
#'
#' Contains:
#' - All clinical data (clinical, follow-up, etc..)
#' - Gene expression metadata (FPKM and Counts)
#' - Mutation data
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"gdc.original"


#' MultiAssayExperiment data
#'
#' To get RNA-Seq assay data, please run loadGDCRnaSeq()
#'
#' It is not included in the data originaly to avoid redundant
#'  data in the package.
#'
#' Contains:
#' - All clinical data (clinical, follow-up, etc..)
#' - Gene expression metadata (assay can be obtained by function 'loadGDCRnaSeq')
#' - Mutation data
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"multiAssay"