#' read computeMatrix output
#' @description Read computeMatrix output into a SummarizedExperiment object
#' @param file count filename of computeMatrix output.
#' @return a SummarizedExperiment object
#' @export
#' @importFrom rjson fromJSON
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @import GenomicRanges
#' @import IRanges
#' @importFrom utils read.delim
#' @examples 
#' file <- system.file("extdata", "count.gz", package= "deepToolsDownstream")
#' se <- importCount(file)
importCount <- function(file){
  stopifnot(length(file)==1)
  f <- gzfile(file, open = "r")
  on.exit(close(f))
  headerL <- readLines(f, n=1)
  header <- fromJSON(substring(headerL, 2))
  d <- read.delim(f, header=FALSE, comment.char = "@")
  anno <- d[, seq.int(6)]
  d <- d[, -seq.int(6)]
  dd <- mapply(function(a, b, ups, body, dws, bin){
    .ele <- d[, a:b]
    U <- ups/bin
    B <- body/bin
    D <- dws/bin
    colnames(.ele) <- c(-1*rev(seq.int(U))*bin + bin/2,
                        paste0(rep(c("B", "D"), c(B, D)),
                               seq.int(B+D)*bin - bin/2))
    .ele
  }, header$sample_boundaries[-length(header$sample_boundaries)]+1, 
  header$sample_boundaries[-1],
  header$upstream,
  header$body,
  header$downstream,
  header$`bin size`,
  SIMPLIFY = FALSE)
  names(dd) <- header$sample_labels
  strand <- as.character(anno$V6)
  strand[strand=="."] <- "*"
  anno <- GRanges(anno$V1, ranges = IRanges(anno$V2, anno$V3, 
                                            names = make.names(anno$V4, unique = TRUE)),
                  score=anno$V5, strand = strand, 
                  group = factor(rep(header$group_labels, diff(header$group_boundaries)),
                                 levels = header$group_labels))
  dd <- lapply(dd, function(.ele){
    if(!identical(rownames(.ele), names(anno))){
      rownames(.ele) <- names(anno)
    }
    .ele
  })
  se <- SummarizedExperiment(assays = dd, rowRanges = anno, 
                             metadata = header)
  se
}
