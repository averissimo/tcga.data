#' Create MultiAssayExperiment object from data
#'
#' @param clinical use custom clinical (that can be pre-processed)
#'
#' @return a MultiAssayExperiment object
#' @export
#'
#' @examples
#' assay <- build.assay()
#' assay[['RNASeq']]
#' assar$vital_status
build.assay <- function(clinical.custom = NULL, 
                        gdc.custom = NULL, 
                        mutation.custom = NULL,
                        rnaseq.custom = NULL) {
  # get clinical data
  if (is.null(clinical.custom)) {
    data(clinical)
    clin <- clinical$all  
  } else {
    clin <- clinical.custom
  }
  
  futile.logger::flog.info('Loading \'Biospecimen\' data...')
  if (is.null(gdc.custom)) {
    data(gdc)
    gdc.custom <- gdc
  }
  
  # get all RNASeq data
  futile.logger::flog.info('Joining \'RNASeq\' data...')
  if (is.null(rnaseq.custom)) {
    rnaseq.custom <- joinRNASeqData()
  }
      
  futile.logger::flog.info('Loading \'Mutation\' data...')
  if (is.null(mutation.custom)) {
    data(mutation)
    mutation.custom <- mutation$count
  }
  
  #
  # Expression data
  
  # map expression data with clinical
  es.map <- data.frame(master = strtrim(colnames(rnaseq.custom), 12), 
                       assay = colnames(rnaseq.custom), 
                       stringsAsFactors = FALSE)
  
  # filter only valid date.. i.e expression that have clinical data
  valid.ix <- es.map$master %in%  clin$bcr_patient_barcode
  valid.dat <- rnaseq.custom[, valid.ix]
  
  sample.barcode <- strtrim(colnames(valid.dat), 16)
  valid.codes <- sample.barcode[sample.barcode %in% gdc$bio.sample$bcr_sample_barcode]
  
  temp.df <- Biobase::AnnotatedDataFrame(gdc$bio.sample[valid.codes,])
  rownames(temp.df) <- colnames(valid.dat)
  
  # build expression set
  es  <- Biobase::ExpressionSet(assayData = valid.dat, phenoData = temp.df)
  
  #
  # Mutation data
  mutation.colnames <- colnames(mutation.custom)  
  valid.ix <- colnames(mutation.custom) %in% clin$bcr_patient_barcode
  
  mut.map <- data.frame(master = mutation.colnames[valid.ix], assay = mutation.colnames[valid.ix])
  
  mut <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mutation.custom))
  
  #
  # Setup to create MultiAssayExperiment object
  
  futile.logger::flog.info('Building Assay...')
  listmap <- list(es.map, mut.map)
  names(listmap) <- c("RNASeq", "Mutation")
  
  dfmap <- MultiAssayExperiment::listToMap(listmap)
  objlist <- list("RNASeq" = es, "Mutation" = mut)
  my.assay <- MultiAssayExperiment::MultiAssayExperiment(objlist, clin, dfmap)
  
  return(my.assay)
}