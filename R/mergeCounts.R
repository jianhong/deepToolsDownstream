#' merge profile from multiple counts file
#' @description Merge profile from multiple counts file by given function
#' @param files count filenames of computeMatrix output.
#' @param FUN function for each column in counts. Could set to 
#' colMeans, colSum.
#' @param na.rm logical. Should missing values (including NaN) be 
#' omitted from the calculations?
#' @return data.frame object.
#' @importFrom reshape2 melt
#' @export
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' mergeCounts(file)
mergeCounts <- function(files, FUN=colMeans, na.rm=TRUE){
  stopifnot(is.function(FUN))
  if(!as.character(substitute(FUN)) %in% c("colMeans", "colSum")){
    stop("FUN must be colMeans or colSum.")
  }
  cnts <- list()
  for(i in seq_along(files)){
    se <- importCount(files[i])
    counts <- assays(se)
    anno <- rowRanges(se)
    header <- metadata(se)
    ## split by group
    group <- anno$group
    anno <- split(anno, group)
    counts <- lapply(counts, function(.ele){
      split(.ele, group)
    })
    counts <- unlist(counts, recursive = FALSE)
    counts <- lapply(counts, FUN = FUN, na.rm=na.rm)
    cnts[[files[i]]] <- do.call(cbind, counts)
  }
  cn <- melt(cnts)
  colnames(cn) <- c("coord", "sample_group", "value", "L1")
  cn
}
