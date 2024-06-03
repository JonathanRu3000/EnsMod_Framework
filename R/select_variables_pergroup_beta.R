#' Select Variables Per Group for Habitat Suitability Modeling
#'
#' This function selects relevant environmental variables from a specified dataset for each group,
#' typically representing different ecological groups or conditions. It performs variable selection by
#' analyzing data using biserial correlation and checking for multicollinearity (VIF) utilizing the
#' SDMworkshop package by BlasBenito (https://github.com/BlasBenito/SDMworkshop).
#'
#' @param groups A list of raster stacks representing different groups of predictors.
#'               All raster objects must have the same CRS (Coordinate Reference System) and extent.
#' @param data An object returned by \code{\link{background_pt}}.
#'
#' @return A list where each element corresponds to a group with selected variables and their
#'         analysis metrics, including results of biserial correlation and VIF values.
#'
#' @examples
#' \dontrun{
#' my_groups <- handle_error(list(clim = clim, notclim = stack(topo, soil_new)))
#' selection <- select_variables_pergroup(groups = groups, data = backgr)
#' print(selection)
#' }
#'
#' @import raster
#' @import dplyr
#' @import stringr
#' @export
select_variables_pergroup <- function(groups = NULL, data = NULL) {
  
  # Create an empty list for the length of groups
  selection <- vector(mode = "list", length = length(groups))
  names(selection) <- names(groups)
  
  # Loop over all elements in groups
  # Extract values for all presences and absences for predictors of a group
  ext <- na.omit(data.frame(occ = data$pr_bg[, 3], raster::extract(groups[[1]], data$pr_bg[, 1:2])))
  
  # Helper function to detect bioclimatic indices from raster patterns
  get_bioclimatic_indices_by_patterns <- function(raster_stack, patterns) {
    if (!inherits(raster_stack, c("RasterStack", "RasterBrick"))) {
      stop("Input must be a RasterStack or RasterBrick")
    }
    layer_names <- names(raster_stack)
    combined_pattern <- stringr::str_c(patterns, collapse = "|")
    indices <- which(stringr::str_detect(layer_names, combined_pattern))
    return(setNames(indices, layer_names[indices]))
  }
  
  # Identify temperature and precipitation indices
  ind_temp <- get_bioclimatic_indices_by_patterns(groups$clim, c("bio1_", "bio2_", "bio3_", "bio4_", "bio5_", "bio6_", "bio7_", "bio8_", "bio9_", "bio10_", "bio11_"))
  ind_prec <- get_bioclimatic_indices_by_patterns(groups$clim, c("bio12_", "bio13_", "bio14_", "bio15_", "bio16_", "bio17_", "bio18_", "bio19_"))
  
  l <- list(temp = ind_temp, prec = ind_prec)
  prec_temp <- vector(mode = "list", length = length(l))
  names(prec_temp) <- names(l)
  
  # select temperature and precipitation related variables
  for (p in seq_along(l)) {
    ext2 <- data.frame(occ = ext$occ, ext[, l[[p]]])[,-2]
    names(ext2)[1] <- "presence"
    
    # run biserialCorrelation analysis and VIF analysis (try to keep the two first variables from bs analysis)
    bc <- biserialCorrelation(ext2, presence.column = "presence")
    av <- autoVIF(ext2[,-1], try.to.keep = bc$df$variable[1:2], verbose = TRUE)
    
    prec_temp[[p]]$bc <- bc
    prec_temp[[p]]$vifnames <- av
  }
  
  extfinal <- ext %>% dplyr::select(all_of(c(prec_temp[[1]]$vifnames, prec_temp[[2]]$vifnames)))
  final <- autoVIF(extfinal)
  
  selection[[1]]$prec_temp <- prec_temp
  selection[[1]]$vifnames <- final
  
  # select nonclimatic variables
  for (i in seq_along(groups)[-1]) {
    ext <- na.omit(data.frame(occ = data$pr_bg[, 3], raster::extract(groups[[i]], data$pr_bg[, 1:2])))
    if ("ID" %in% names(ext)) {
      ext <- ext %>% dplyr::select(!"ID")
    }
    names(ext)[1] <- "presence"
    
    bc <- biserialCorrelation(ext, presence.column = "presence")
    av <- autoVIF(ext[,-1], try.to.keep = bc$df$variable[1:2])
    
    selection[[i]]$bc <- bc
    selection[[i]]$vifnames <- av
  }
  
  selection$occs_which_na <- occs$stats$occs_final - nrow(ext %>% filter(presence == 1))
  
  return(selection)
}
