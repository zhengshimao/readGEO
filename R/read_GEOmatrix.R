#' Read and Process GEO Expression Matrix
#'
#' This function reads a GEO-formatted expression matrix file, extracts expression
#' and phenotype data, and returns them in a structured list.
#'
#' @param f_matrix Character. Path to the GEO expression matrix file.
#'
#' @return A list containing:
#' \describe{
#'   \item{exp}{Data frame. Expression matrix with probe IDs in the first column.}
#'   \item{clinical_info}{Data frame. Phenotype data associated with the samples.}
#'   \item{eSet}{ExpressionSet object. The original GEOquery object containing the full data.}
#' }
#'
#' @details
#' This function uses `GEOquery::getGEO` to load the GEO file and retrieve expression and
#' clinical data. The expression matrix is sorted to match the order of samples in the
#' clinical information, ensuring consistency.
#'
#' @examplesIf FALSE
#' f_matrix <- system.file("extdata/GSE106172/matrix","GSE106172_series_matrix.txt.gz",package = "readGEO")
#' result <- read_GEOmatrix(f_matrix = f_matrix)
#'
#' @importFrom GEOquery getGEO
#' @export

read_GEOmatrix <- function(f_matrix){
  eSet <- GEOquery::getGEO(filename = f_matrix,
                           GSElimits = NULL,
                           AnnotGPL = F,
                           getGPL = F)

  exp_ori <- eSet@assayData$exprs
  # exp_ori0 <- Biobase::exprs(eSet)

  pd <- eSet@phenoData@data
  # pd <- Biobase::pData(eSet)

  # sort expression according to clinical infor
  if( ! identical(rownames(pd),colnames(exp_ori0)) ) {
    exp_ori0 = exp_ori0[,match(rownames(pd),colnames(exp_ori0))]
  }

  exp_ori <- data.frame(probe_id= rownames(exp_ori0),exp_ori0)
  stopifnot(identical(exp_ori$probe_id, exp_ori$probe_id))
  return(list(exp = exp_ori, clinical_info = pd, eSet = eSet))
}


