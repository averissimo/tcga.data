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
