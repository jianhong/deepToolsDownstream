#' plot heatmap
#' @description Plot a heatmap for counts.
#' @param se a SummarizedExperiment object from \link{importCount}.
#' @param facet This will be passed to 
#' \link[ggplot2:facet_grid]{facet_grid}:rows.
#' It can be set to NULL.
#' @param transformFUN The transform function of values.
#' @param pseudo Pseudo value will be add to values before apply transformFUN.
#' @param orderFUN function used to reorder Y axis of the heatmap
#' @param orderBy the sample name or index used to order Y axis
#' @param fill_gradient color fill gradient function. 
#' see \link[ggplot2:scale_colour_gradient]{scale_colour_gradient}.
#' @param boderColor The color of the heatmap cell border.
#' @param xaxis_breaks,xaxis_label xaxis breaks and labels.
#'  see \link[ggplot2:scale_continuous]{scale_x_continuous}.
#' @param yaxis_breaks,yaxis_label xaxis breaks and labels.
#'  see \link[ggplot2:scale_x_discrete]{scale_y_discrete}.
#' @importFrom SummarizedExperiment rowRanges assays
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @importFrom reshape2 melt dcast
#' @importFrom ggplot2 ggplot aes_string geom_line theme_classic 
#' scale_x_continuous facet_wrap xlab ylab scale_colour_gradient2
#' scale_fill_gradient scale_y_discrete theme element_blank geom_tile
#' facet_grid
#' @importFrom stats as.formula
#' @export
#' @return ggplot object
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' se <- importCount(file)
#' plotHeatmap(se, yaxis_breaks="100033817", yaxis_label="geneA")
plotHeatmap <- function(se,
                        facet=as.formula("group~sample"), 
                        transformFUN=log2,
                        pseudo=1,
                        orderFUN=rowMeans,
                        orderBy=1,
                        fill_gradient=scale_fill_gradient(
                          low = "blue", high = "red"),
                        boderColor=NA,
                        xaxis_breaks=NULL,
                        xaxis_label=NULL,
                        yaxis_breaks=NULL,
                        yaxis_label=NULL){
  stopifnot(is(se, "SummarizedExperiment"))
  if(length(facet)>0){
    if(is.character(facet)){
      facet <- facet[1]
      if(!facet %in% c("sample", "group")){
        stop("facet must be sample or group")
      }
    }
    setFacet <- TRUE
  }else{
    setFacet <- FALSE
  }
  counts <- assays(se)
  if(is.numeric(orderBy)){
    orderBy <- names(counts)[orderBy]
  }
  orderBy <- orderBy[orderBy %in% names(counts)]
  if(length(orderBy)==0){
    stop("orderBy is not proper names of samples!")
  }
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
  ## 
  d_melt <- melt(lapply(counts, as.matrix))
  colnames(d_melt) <- c("annoID", "coord", "value", "sample_group")
  d_melt$sample <- factor(old_name[sub("^(.).*$", "\\1", 
                                       as.character(d_melt$sample_group))])
  d_melt$group <- factor(sub("^..", "", as.character(d_melt$sample_group)))
  d_melt$x <- as.numeric(sub("^B|D", "", as.character(d_melt$coord)))
  if(is.function(transformFUN)){
    d_melt$value <- transformFUN(d_melt$value+pseudo)
  }
  d_melt$annoID <- as.character(d_melt$annoID)
  dataForReorder <- d_melt[d_melt$sample %in% orderBy,
                           c("annoID", "coord", "value", "sample")]
  dataForReorder <- split(dataForReorder[, c("annoID", "coord", "value")],
                          dataForReorder$sample, drop = TRUE)
  dataForReorder <- dataForReorder[orderBy]
  dataForReorder <- lapply(dataForReorder, function(.ele)
    dcast(.ele, formula = as.formula("annoID~coord")))
  annoID <- unique(unlist(lapply(dataForReorder, function(.ele) .ele$annoID)))
  dataForReorder <- lapply(dataForReorder, function(.ele){
    name <- .ele$annoID
    .ele <- .ele[, -1]
    .ele <- as.matrix(.ele)
    .ele[is.infinite(.ele)] <- NA
    if("na.rm" %in% names(formals(orderFUN))){
      values <- orderFUN(.ele, na.rm = TRUE)
    }else{
      values <- orderFUN(.ele)
    }
    names(values) <- name
    values[annoID]
  })
  d_melt$annoID <- 
    factor(d_melt$annoID,
           levels = annoID[do.call(function(...)
             order(..., decreasing = FALSE), dataForReorder)])
  
  lines <- unique(c(0, d_melt$x[which(grepl("^D", 
                                            as.character(d_melt$coord)))[1]] -
                      header$`bin size`[1]/2))
  
  if(length(xaxis_breaks)==0){
    xaxis_breaks <- c(-1*header$upstream[1], lines, 
                      sum(lines)+header$downstream[1])
  }
  if(length(xaxis_label)!=length(xaxis_breaks)){
    xaxis_breaks <- c(-1*header$upstream[1], lines,
                      sum(lines)+header$downstream[1])
    xaxis_label <- c(-1*header$upstream[1], c("TSS", "TES")[seq_along(lines)],
                     header$downstream[1])
  }
  
  p <- ggplot(d_melt, aes_string(x="x", y="annoID", fill="value")) + 
    geom_tile(color=boderColor) + fill_gradient
  if(setFacet){
    p <- p + facet_grid(rows = facet, scales="free")
  }
  p <- p +
    xlab("distance (bp)") + ylab(header$`bin avg type`) + 
    scale_x_continuous(breaks = xaxis_breaks, 
                       labels = xaxis_label) +
    scale_y_discrete(breaks = yaxis_breaks,
                     labels = yaxis_label) +
    theme_classic() + 
    theme(axis.line.y=element_blank())
  p
}
