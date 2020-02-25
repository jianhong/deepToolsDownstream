#' plot profiles
#' @description Plot a profile for counts.
#' @param se a SummarizedExperiment object from \link{importCount}.
#' @param loessSmooth Use \link[stats:scatter.smooth]{loess.smooth} 
#' to smooth curve or not.
#' @param span smoothness parameter for loess.
#' @param facet group or sample. This will be passed to \link[ggplot2:facet_wrap]{facet_wrap}.
#' @param xaxis_breaks,xaxis_label xaxis breaks and labels. see \link[ggplot2:scale_continuous]{scale_x_continuous}.
#' @importFrom SummarizedExperiment rowRanges assays
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @importFrom stats loess.smooth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string geom_line theme_classic scale_x_continuous facet_wrap xlab ylab
#' @export
#' @return ggplot object
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' se <- importCount(file)
#' plotProfile(se)
plotProfile <- function(se, loessSmooth = TRUE, span = 1/25,
                        facet="group", 
                        xaxis_breaks=NULL,
                        xaxis_label=NULL){
  stopifnot(is(se, "SummarizedExperiment"))
  facet <- facet[1]
  if(!facet %in% c("sample", "group")){
    stop("facet must be sample or group")
  }
  counts <- assays(se)
  anno <- rowRanges(se)
  header <- metadata(se)
  ## split by group
  group <- anno$group
  anno <- split(anno, group)
  counts <- lapply(counts, function(.ele){
    split(.ele, group)
  })
  ## unlist counts and rename the samples by letters.
  old_name <- header$sample_labels
  if(length(counts)>26){
    stop("Can not handle samples more than 26.")
  }
  names(old_name) <- LETTERS[seq_along(old_name)]
  names(counts) <- LETTERS[seq_along(counts)]
  counts <- unlist(counts, recursive = FALSE)
  ## colSum
  dd_norm <- lapply(counts, colMeans, na.rm=FALSE)
  if(loessSmooth){
    dd_norm <- lapply(dd_norm, function(.ele){
      loess.smooth(x = seq_along(.ele),
                   y = .ele,
                   span = span,
                   family = "gaussian",
                   evaluation = length(.ele))$y
    })
  }
  dd_norm <- mapply(dd_norm, counts, FUN = function(a, b){
    names(a) <- colnames(b)
    a
  }, SIMPLIFY = FALSE)
  d_melt <- do.call(cbind, dd_norm)
  d_melt <- melt(d_melt)
  colnames(d_melt) <- c("coord", "sample_group", "value")
  d_melt$sample <- factor(old_name[sub("^(.).*$", "\\1", 
                                       as.character(d_melt$sample_group))])
  d_melt$group <- factor(sub("^..", "", as.character(d_melt$sample_group)))
  d_melt$x <- as.numeric(sub("^B|D", "", as.character(d_melt$coord)))
  lines <- unique(c(0, d_melt$x[which(grepl("^D", as.character(d_melt$coord)))[1]] - header$`bin size`[1]/2))
  x <- ifelse(facet=="group", "sample", "group")
  if(length(xaxis_breaks)==0){
    xaxis_breaks <- c(-1*header$upstream[1], lines, sum(lines)+header$downstream[1])
  }
  if(length(xaxis_label)!=length(xaxis_breaks)){
    xaxis_breaks <- c(-1*header$upstream[1], lines, sum(lines)+header$downstream[1])
    xaxis_label <- c(-1*header$upstream[1], c("TSS", "TES")[seq_along(lines)],
                     header$downstream[1])
  }
  p <- ggplot(d_melt, aes_string(x="x", y="value", color=x)) + 
    geom_line() + facet_wrap(facets = facet) +
    xlab("distance (bp)") + ylab(header$`bin avg type`) + 
    scale_x_continuous(breaks = xaxis_breaks, 
                       labels = xaxis_label) +
    theme_classic() 
}