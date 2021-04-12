#' change sample/group labels
#' @description Change sample/group labels of imported data.
#' @param se a SummarizedExperiment object from \link{importCount}.
#' @param attr attribute to be changed. Could be sample_labels or group_labels.
#' @param values values to be set. 
#' @export
#' @return a SummarizedExperiment object
#' @importFrom SummarizedExperiment assays<- rowRanges<-
#' @importFrom S4Vectors metadata<-
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' se <- importCount(file)
#' se <- relabels(se, "group_labels", c("gp1", "gp2"))
#' plotProfile(se)
relabels <- function(se, attr="sample_labels", values){
  stopifnot(is(se, "SummarizedExperiment"))
  attr <- attr[1]
  stopifnot(attr %in% c("sample_labels", "group_labels"))
  header <- metadata(se)
  if(length(values)!=length(header[[attr]])){
    stop("length of values is not identical to the data")
  }
  header[[attr]] <- values
  metadata(se) <- header
  if(attr=="sample_labels"){
    names(assays(se)) <- values
  }else{#attr=="group_labels"
    rowRanges(se)$group <- factor(rep(header$group_labels, diff(header$group_boundaries)),
                                  levels = header$group_labels)
  }
  se  
}