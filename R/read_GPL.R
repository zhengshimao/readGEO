#' Read and Process GEO Platform Annotation
#'
#' This function reads a GEO platform (GPL) annotation file in SOFT format, processes
#' probe-to-gene mappings, and optionally saves the result in `.rda` format.
#'
#' @param f_soft Character. Path to the GEO platform annotation file in SOFT format.
#' @param drop_na Logical. Whether to drop rows with missing gene symbols. Default is \code{TRUE}.
#' @param save_ids_only Logical. If \code{TRUE}, only the processed probe-to-gene mapping is returned.
#' Otherwise, additional raw and intermediate data are included in the output. Default is \code{TRUE}.
#' @param save_rda Character. Directory where the `.rda` file should be saved. Default is \code{NULL}, meaning no file will be saved.
#' @param rda_name Character. Name of the `.rda` file (without extension). If \code{NULL}, it defaults to the platform ID with a `_bioc` suffix.
#'
#' @return If \code{save_ids_only} is \code{TRUE}, a data frame with probe IDs and their associated gene symbols.
#' Otherwise, a list containing:
#' \describe{
#'   \item{ids}{Data frame. Processed probe-to-gene mappings.}
#'   \item{gpl_ori}{Data frame. Original platform annotation table.}
#'   \item{gpl_anno}{Data frame. Intermediate processed annotation table.}
#' }
#'
#' @details
#' The function parses the platform annotation table using `GEOquery::parseGEO` and extracts relevant columns for
#' probe-to-gene mappings. It applies further filtering and formatting to clean the data, ensuring consistency.
#'
#' If \code{save_rda} is specified, the processed data is saved as an `.rda` file in the given directory, mimicking
#' the storage format used in the \code{AnnoProbe} package.
#'
#' @examplesIf FALSE
#' f_soft <- "./inst/extdata/GSE106172/soft/GSE106172_family.soft.gz"
#' ids <- read_GPL(f_soft = f_soft, drop_na = TRUE, save_ids_only = TRUE)
#' full_data <- read_GPL(f_soft = f_soft, save_ids_only = FALSE)
#'
#' @importFrom GEOquery parseGEO
#' @importFrom dplyr select filter
#' @importFrom stringr str_extract str_remove_all
#' @importFrom stats na.omit
#' @export

read_GPL <- function(f_soft, drop_na = TRUE, save_ids_only = TRUE,save_rda = NULL, rda_name = NULL){
  # 读入探针注释
  gpl <- GEOquery::parseGEO(fname = f_soft,
                            GSElimits = NULL,
                            AnnotGPL = F,
                            getGPL = F)
  gpl_ori <- gpl@gpls[[1]]@dataTable@table # 获取平台注释表格
  # gpl_ori <- gpl@gpls[[1]]@dataTable %>% GEOquery::Table()
  # gpl_ori <- gpl@gpls$GPL16570@dataTable %>% GEOquery::Table()
  gpl_number <- gpl@header$platform_id # 获取平台编号

  if(F){
    txt = data.table::fread(f_soft, sep = "")[[1]]
    txt = txt[txt != ""]
    # return(.parseGPLTxt(txt))
    gpl_ori <- GEOquery:::.parseGPLTxt(txt) %>% GEOquery::Table()
  }
  # 进一步处理探针注释
  gpl_anno <- gpl_ori %>%
    dplyr::select(ID,gene_assignment) %>%
    dplyr::filter(gene_assignment != "---") #%>%
  # tidyr::separate(gene_assignment,c("drop","symbol"),sep="//") %>%
  # dplyr::select(-drop)

  gpl_anno["symbol"] <- stringr::str_extract(gpl_anno$gene_assignment ,pattern = "//.*?//") %>%
    stringr::str_remove_all(pattern = "//| ")

  ids <- gpl_anno %>% dplyr::select("ID","symbol")
  colnames(ids) <- c("probe_id","symbol")
  if(drop_na){
    ids <- stats::na.omit(ids)
  }

  # Save GPL probe data
  if(is.null(rda_name)){
    rda_name <- paste0(gpl_number,"_bioc")
  }
  assign(x = rda_name, ids)
  if(!is.null(save_rda)){
    if(!dir.exists(save_rda)){
      save_rda <- "./"
      message(paste0("./",rda_name))
    }
    # 保存Rdata #模仿AnnoProbe::idmap存储形式
    save(list = rda_name, file = paste0(save_rda,"/",rda_name,".rda"))
  }
  if(save_ids_only){
    return(ids)
  }else{
    return(list(ids = ids, gpl_ori = gpl_ori, gpl_anno = gpl_anno))
  }
}
