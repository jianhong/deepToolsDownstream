#' write computeMatrix output
#' @description Write a SummarizedExperiment object into computeMatrix output. 
#' @param se a SummarizedExperiment object from \link{importCount}.
#' @param filename output filename
#' @return zipped output filename
#' @export
#' @importFrom R.utils gzip
#' @importFrom rjson toJSON
#' @importFrom SummarizedExperiment rowRanges assays
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @importFrom utils write.table
#' @import GenomicRanges
#' @import BiocGenerics
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' se <- importCount(file)
#' exportCount(se)
exportCount <- function(se, filename = tempfile()){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  header <- metadata(se)
  if(length(header$group_labels)==1) header$group_labels <- list(header$group_labels)
  if(length(header$sample_labels)==1) header$sample_labels <- list(header$sample_labels)
  headerL <- toJSON(header)
  anno <- as.data.frame(rowRanges(se))
  out <- do.call(cbind, assays(se))
  out[is.na(out)] <- 0 ## remove NA.
  out <- cbind(anno[, c("seqnames", "start", "end")], name=rownames(anno), 
               anno[, c("score", "strand")], round(out, digits = 6))
  filen <- sub(".gz", "", filename)
  writeLines(paste0("@", headerL), filen)
  write.table(out, file = filen, append = TRUE, quote=FALSE,
              sep="\t", row.names = FALSE, col.names = FALSE)
  gzip(filename = filen, overwrite=TRUE)
  return(paste0(filen, ".gz"))
}
