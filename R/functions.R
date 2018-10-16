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
  return(TCGAutils::TCGAbarcode(fullBarcode))
}

#' Build list of matrices with different types of tissues
#'
#' @param source.name data type
#' @param assay.data data from MultiAssayExperiment
#'
#' @return data divided by tissue type, with clinical and source.name assay
#' @export
#'
#' @examples
#' data('multiAssay')
#' fpkm.data <- build.matrix('RNASeqFPKM', multiAssay)
#' fpkm.per.tissue <- fpkm.data$data
#' fpkm.clinical   <- fpkm.data$clinical
#' names(fpkm.per.tissue)
build.matrix <- function(source.name, assay.data) {
  
  ret.list  <- list()
  cli.list  <- list()
  full.code <- list()
  
  if (source.name %in% names(assay.data@ExperimentList)) {
    if (source.name == 'Mutation') {
      return(assay.data[[source.name]]@assays[['counts']])
    } else if (source.name %in% c('RNASeqFPKM', 'RNASeqCounts')) {
      suppressMessages(
        new.assay <- assay.data[,,source.name]
      )
      for (sample_id in unique(new.assay[[source.name]]@phenoData$definition)) {
        sample.id <- gsub(' ', '.', sample_id) %>% tolower()
        # keep only one type of sample
        suppressMessages(
          tmp.assay <- new.assay[,new.assay[[source.name]]$definition == sample_id,source.name]
          )
        ret.list[[sample.id]] <- Biobase::exprs(tmp.assay[[source.name]])
        cli.list[[sample.id]] <- tmp.assay@colData
        
        full.code[[sample.id]]          <- colnames(ret.list[[sample.id]])
        colnames(ret.list[[sample.id]]) <- strtrim(colnames(ret.list[[sample.id]]), 12)
      }
      
      futile.logger::flog.info('Individuals per sample type for %s', source.name)
      for (ix.name in names(ret.list)) {
        futile.logger::flog.info('  * %s: %d', ix.name, ncol(ret.list[[ix.name]]))
      }
      
      cli.list$all <- new.assay@colData
      return(list(clinical = cli.list, data = ret.list, original.codes = full.code))
    } else {
      stop('Source must be one of the assays in multiAssay data object')  
    }
  } else {
    stop('Source must be one of the assays in multiAssay data object')
  }
}
