#' Get Participant part of TCGA Barcode
#'
#' TCGA-XX-XXXX
#'
#' @param fullBarcode string, or  vector of strings
#'
#' @return participant's barcode
#' @export
#'
#' @examples
#' getParticipantCode('TCGA-3C-AAAU-01A-11R-A41B-07')
#' getParticipantCode(c('TCGA-3C-AAAU-01A-11R-A41B-07', 'TCGA-3C-AALI-01A-11R-A41B-07'))
getParticipantCode <- function(fullBarcode) {
  tcga.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4}).*'
  return(gsub(tcga.barcode.pattern, '\\1', fullBarcode))
}


#' Get Sample Type part of TCGA Barcode
#'
#' XX - where 01 is solid tumor, 06 is metasteses, etc..
#'
#' See [Code Tables Report](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) for a complete list of sample codes.
#'
#' @param fullBarcode number, or  vector of numbers
#'
#' @return sample code
#' @export
#'
#' @examples
#' getParticipantCode('TCGA-3C-AAAU-01A-11R-A41B-07')
#' getParticipantCode(c('TCGA-3C-AAAU-01A-11R-A41B-07', 'TCGA-3C-AALI-01A-11R-A41B-07'))
getSampleTypeCode <- function(fullBarcode) {
  tcga.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4})-([0-9]{2}).*'
  return(gsub(tcga.barcode.pattern, '\\2', fullBarcode))
}

#' Joins all expression data in a single matrix
#'
#' This is not stored directly as it would be duplicate information
#'  which takes a lot of storage
#'
#' @return a matrix with gene expression levels and all tissue samples
#' @export
#'
#' @examples
joinRNASeqData <- function() {
  # load tissue data
  data(fpkm.per.tissue)

  # iterate on all and join in a single matrix
  arrays.in.dat <- sapply(seq(length(tissue)), function(ix){ is.array(tissue[[ix]]) })

  out.dat <- c()
  for ( ix in which(arrays.in.dat)) {
    out.dat <- cbind(out.dat, tissue[[ix]])
  }
  return(out.dat)
}

#' Get rnaseq assay data to GDC object
#'
#' This cannot be cached in the package as it takes too much space and
#'  would be redundant with fpkm.per.tissue data
#'
#' @param project
#' @param workflow.type
#'
#' @return
#' @export
#'
#' @examples
loadGDCRnaSeq <- function() {
  data(gdc)

  temp.dat <- joinRNASeqData()
  rownames(temp.dat) <- rownames(gdc$rnaseq)
  temp.dat <- temp.dat[,colnames(gdc$rnaseq)]

  gdc$rnaseq@assays[[1]] <- temp.dat

  return(gdc)
}
