#' filter computeMatrix output
#' @description Filter output of importCount. 
#' @param se a SummarizedExperiment object from \link{importCount}.
#' @param subset filter condition for rows.
#' @return a SummarizedExperiment object
#' @export
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' se <- importCount(file)
#' library(SummarizedExperiment)
#' keep <- rowMeans(assays(se)[[1]], na.rm = TRUE) < 2 ## arbitory number
#' nrow(se)
#' se <- filterCount(se, subset=keep)
#' nrow(se)
#' table(keep)
filterCount <- function(se, subset){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  if(!missing(subset)){
    header <- metadata(se)
    group_boundaries <- rep(header$group_labels, diff(header$group_boundaries))
    group_boundaries <- group_boundaries[subset]
    header$group_labels <- unique(group_boundaries)
    group_boundaries <- table(group_boundaries)[header$group_labels]
    header$group_boundaries <- unname(cumsum(c(0, group_boundaries)))
    se <- subset(se, subset=subset)
    metadata(se) <- header
  }
  se
}
