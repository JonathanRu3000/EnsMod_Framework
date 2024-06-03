#' Perform Tukey's HSD Test for selecting most significcant soil layers
#'
#' This function performs Tukey's Honest Significant Difference (HSD) test on soil layers for a given set of points and growth form. It identifies and removes non-significant soil layers based on the specified growth form or a custom pattern.
#' This function works only for soil layers from https://soilgrids.org/
#'
#' @param pattern A string specifying the soil depth layers to include. If NULL, a default pattern is selected based on the growth form. "0.5|5.15|15.30" selects soil layers from 0-5 cm, 5-15cm and 15-30 cm. Soil layers should be stored with the information for depth in the following way 0.5 for 0-5cm, 15.30 for 15-30 cm etc.
#' @param points An object created by \code{\link{background_pt}}.
#' @param soil A raster stack of soil layers from https://soilgrids.org/.
#' @param growthform A string specifying the growth form. Options are "small_annual", "annual", "shrub", or "tree".
#'
#' @return A modified raster stack with non-significant soil layers removed.
#'
#' @examples
#' \dontrun{
#' result <- tukey_soil(pattern = NULL, points = backgr, soil = soil, growthform = "small_annual")
#' print(result)
#' }
#'
#' @import dplyr
#' @import raster
#' @import reshape2
#' @export
#'
#' @details
#' This function uses Tukey's HSD test to compare the means of different soil layers and identify those that are not significantly different from each other. The test is conducted at a confidence level of 99% (conf.level = 0.99) and a significance level of 5% (p.adj > 0.05) for adjusted p-values.
#'
tukey_soil <- function(pattern = NULL,
                       points = NULL,
                       soil = NULL,
                       growthform = c(
                         "small_annual",
                         "annual",
                         "shrub",
                         "tree"
                       )) {


  # Match the growth form argument to the available options
  growthform <- match.arg(growthform)

  # Assign default pattern based on growth form if not provided
  if (is.null(pattern)) {
    pattern <- switch(growthform,
                      "small_annual" = "0.5|5.15|15.30",
                      "annual" = "0.5|5.15|15.30|30.60",
                      "shrub" = "0.5|5.15|15.30|30.60|60.100",
                      "tree" = "0.5|5.15|15.30|30.60|60.100|100.200")
  }

  # Select the desired soil layers
  dsoil <- soil[[grep(pattern, names(soil))]]

  # Extract values from raster and create data frame
  soilt <- data.frame(occ = points$pr_bg[, 3], raster::extract(dsoil, points$pr_bg[, -3]))

  # Get unique soil metrics
  names1 <- names(soilt)[-1]
  names2 <- unique(sapply(strsplit(names1, "_"), `[`, 2))

  # Perform Tukey's HSD test for each soil metric
  significant_layers <- sapply(names2, function(metric) {
    temp <- soilt %>% dplyr::select("occ", contains(metric))
    m <- melt(temp, id = "occ")[, -1]
    m$variable <- as.factor(m$variable)
    a <- aov(value ~ variable, m)
    tukey <- data.frame(TukeyHSD(x = a, conf.level = 0.99)$variable)
    nonsig <- tukey %>% filter(p.adj > 0.05)
    if (nrow(nonsig) > 0) {
      unique(sapply(strsplit(rownames(nonsig), "_"), `[`, 1:2)) %>%
        apply(2, paste, collapse = "_")
    } else {
      NULL
    }
  })

  # select the significant layers
  significant_layers <- unlist(significant_layers, use.names = FALSE)
  if (length(significant_layers) > 0) {
    soil_new <- dsoil %>% raster::dropLayer(setdiff(names(dsoil), significant_layers))
  } else {
    soil_new <- dsoil
  }

  return(soil_new)
}
