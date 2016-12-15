#' Gene expression levels per tissue type.
#'
#' List of matrices with RNASeq v2 FPKM (rows: genomic ranges) x (columns: patients)
#'
#' FPKM: fragments per kilobase per million
#' = [# of fragments]/[length of transcript in kilo base]/[million mapped reads]
#'
#' For all FPKM expression data matrix, run `joinRNASeqData` dataset
#'
#' @format A list with FPKM gene expression data.
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"fpkm.per.tissue"

#' Patient's barcode per tissue type
#'
#' A dataset containing a list of vectors with TCGA barcodes per tissue type.
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"fpkm.per.tissue.barcode"

#' Clinical data per tissue type
#'
#' A dataset containing a list of dataframes with clinical data per tissue type.
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"clinical"

#' Genomic ranges with description for genes in tissue matrices
#'
#' Contains ensembl id and external gene names.
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"gene.ranges"

#' List with mutation data
#'
#' Contains:
#' - Count of mutations per each pair of gene/case
#' - Original data from gdc (filtered out duplicate rows)
#' - Duplicated rows that were filtered
#'
#' @source \url{https://gdc-portal.nci.nih.gov/projects/TCGA-BRCA}
"mutation"

#' Original data from GDC
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
"gdc"
