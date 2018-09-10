#' Create MultiAssayExperiment object from data
#'
#' @param clinical use custom clinical (that can be pre-processed)
#'
#' @return a MultiAssayExperiment object
#' @export
build.assay <- function(clinical.custom, 
                        gdc.custom, 
                        mutation.custom,
                        rnaseq.fpkm.custom,
                        rnaseq.counts.custom) {

  #
  # Expression data (FPKM)
  
  # map expression data with clinical
  es.map.fpkm <- data.frame(master = strtrim(colnames(rnaseq.fpkm.custom), 12), 
                            assay  = colnames(rnaseq.fpkm.custom), 
                            stringsAsFactors = FALSE)
  
  # filter only valid date.. i.e expression that have clinical data
  valid.ix.fpkm  <- es.map.fpkm$master %in%  clinical.custom$bcr_patient_barcode
  valid.dat.fpkm <- rnaseq.fpkm.custom[, valid.ix.fpkm]
  
  sample.barcode.fpkm <- strtrim(colnames(valid.dat.fpkm), 16)
  valid.codes.fpkm    <- sample.barcode.fpkm[sample.barcode.fpkm %in% gdc$bio.sample$bcr_sample_barcode]
  
  fpkm.df           <- Biobase::AnnotatedDataFrame(gdc$bio.sample[valid.codes.fpkm,])
  rownames(fpkm.df) <- colnames(valid.dat.fpkm)
  
  # build expression set
  es.fpkm <- Biobase::ExpressionSet(assayData = valid.dat.fpkm, phenoData = fpkm.df)
  
  #
  # Expression data (Counts)
  
  # map expression data with clinical
  es.map.counts <- data.frame(master = strtrim(colnames(rnaseq.counts.custom), 12), 
                       assay = colnames(rnaseq.counts.custom), 
                       stringsAsFactors = FALSE)
  
  # filter only valid date.. i.e expression that have clinical data
  valid.ix.counts <- es.map.counts$master %in%  clinical.custom$bcr_patient_barcode
  valid.dat.counts <- rnaseq.counts.custom[, valid.ix.counts]
  
  sample.barcode.counts <- strtrim(colnames(valid.dat.counts), 16)
  valid.codes.counts    <- sample.barcode.counts[sample.barcode.counts %in% gdc$bio.sample$bcr_sample_barcode]
  
  counts.df <- Biobase::AnnotatedDataFrame(gdc$bio.sample[valid.codes.counts,])
  rownames(counts.df) <- colnames(valid.dat.counts)
  
  # build expression set
  es.counts  <- Biobase::ExpressionSet(assayData = valid.dat.counts, phenoData = counts.df)
  
  #
  # Mutation data
  mutation.colnames <- colnames(mutation.custom)  
  valid.ix          <- colnames(mutation.custom) %in% clinical$bcr_patient_barcode
  
  mut.map <- data.frame(master = mutation.colnames[valid.ix], assay = mutation.colnames[valid.ix])
  mut     <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mutation.custom))
  
  #
  # Setup to create MultiAssayExperiment object
  
  futile.logger::flog.info('Building Assay...')
  listmap <- list(es.map.fpkm, es.map.counts, mut.map)
  names(listmap) <- c("RNASeqFPKM", "RNASeqCounts", "Mutation")
  
  dfmap <- MultiAssayExperiment::listToMap(listmap)
  objlist <- list("RNASeqFPKM" = es.fpkm, "RNASeqCounts" = es.fpkm, "Mutation" = mut)
  my.assay <- MultiAssayExperiment::MultiAssayExperiment(objlist, clinical.custom, dfmap)
  
  return(my.assay)
}