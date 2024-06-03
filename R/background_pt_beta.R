#' Generate Background Points for Habitat Suitability Modeling
#'
#' This function generates background points for habitat suitability modeling. Background points are taken from a calibration area (here defined by a convex hull with buffer)
#' (see, e.g., Hirzel & Le Lay, 2008). Additionally, this function does sampling bias correction using a so-called thickening approach (Barber et al., 2022)
#' by creating a density grid based on occurrences and selecting background points according to density weights from within the calibration area.
#' The calibration area is created using the \code{calib_area} function from the \code{flexsdm} package v1.3.3 (Velazco et al., 2022). A density grid is created using the \code{sp.kde} function from the \code{spatialEco} package (Evans et al. 2017).
#'
#' @param data An object created by \code{\link{prepare_occs}} with the cleaned occurrence data.
#' @param rast A raster object used for generating background points.
#' @param n_bg A numeric value indicating the proportion of background points to be generated relative to the non-NA cells in the raster.
#' @param width_cal_area A value indicating the buffer around occurrences for calculating the calibration area in km. (e.g., 100000 = 100 km)
#' @return A list with the following components:
#'   \describe{
#'     \item{bg}{A data frame of generated background points.}
#'     \item{pr_bg}{A data frame combining presence and background points.}
#'     \item{points_pr_bg}{A `SpatialPointsDataFrame` object for the combined presence and background points.}
#'     \item{area_M}{The calibration area as an `sf` object.}
#'     \item{isna}{A data frame of points with NA values.}
#'     \item{na_occs}{A data frame of occurrence points with NA values.}
#'     \item{na_bg}{A data frame of background points with NA values.}
#'   }
#'
#' @references Velazco, S. J. E., Rose, M. B., de Andrade, A. F. A., Minoli, I. & de Alvarenga, A. A. (2022). flexsdm: An r package for supporting a comprehensive and flexible species distribution modelling workflow. Methods in Ecology and Evolution, 13(8), 1661-1669. doi:10.1111/2041-210X.13874
#' @references Evans, J. S., Murphy, M. A., & Ram, K. (2021). spatialEco: Spatial Analysis and Modelling Utilities. R package version 1.3-1. Retrieved from http://archive.linux.duke.edu/cran/web/packages/spatialEco/spatialEco.pdf
#' @examples
#' \dontrun{
#' # Example usage
#' result <- background_pt(data = occs, rast = clim[[1]], n_bg = 0.1, width_cal_area = 100000)
#' print(result)
#' }
#'
#' @import dplyr
#' @import sf
#' @import terra
#' @import raster
#' @importFrom spatialEco sp.kde
#' @importFrom dismo randomPoints
#' @importFrom sp SpatialPointsDataFrame CRS
#' @export
background_pt <- function(data = NULL, rast = NULL, n_bg = 0.1, width_cal_area = 100000) {

  # Convert occurrences to simple feature
  sf_occs <- sf::st_as_sf(data$occs, coords = c("decimallongitude", "decimallatitude"),
                          crs = CRS("+init=epsg:4326"))

  # Calculate calibration area
  area_M <- flexsdm::calib_area(data$occs,
                                x = "decimallongitude",
                                y = "decimallatitude",
                                method = c("bmcp", width = width_cal_area),
                                crs = "+init=epsg:4326")

  # Crop raster layer to the extent of calibration area
  rs <- terra::rast(rast) %>%
    terra::crop(., area_M)

  # Calculate the number of background points
  nu <- terra::global(rs, fun = "notNA")$notNA * n_bg

  # Convert raster to `raster` package format for compatibility
  rs <- raster::raster(rs)

  # Calculate density grid
  tgb_kde <- spatialEco::sp.kde(x = sf::as_Spatial(sf_occs),
                                newdata = rs,
                                standardize = TRUE,
                                scale.factor = 10000)

  # Apply minimum density threshold and mask with calibration area
  tgb_kde[tgb_kde < 0.01] <- 0.01
  tgb_kde <- terra::rast(tgb_kde) %>%
    terra::mask(., area_M)

  # Generate background points from density grid
  bg <- data.frame(dismo::randomPoints(raster::raster(tgb_kde),
                                       n = nu,
                                       p = data$points,
                                       prob = TRUE),
                   occ = rep(0, nu))

  # Set column names to match presence data
  colnames(bg) <- colnames(data$occs)

  # Combine presences and absences in a common data table
  pr_bg <- rbind(bg, data$occs)

  # Extract values from raster and filter out NA values
  temp <- data.frame(pr_bg, value = raster::extract(rast, pr_bg[, c(1:2)]))
  final <- temp %>%
    filter(!is.na(value)) %>%
    dplyr::select(-value)

  # Separate NA values for occurrences and background points
  isna <- temp %>% filter(is.na(value))
  isna_occ <- isna %>% filter(occ == 1)
  isna_back <- isna %>% filter(occ == 0)

  # Create spatial points data frame for combined presences and absences
  spdf_pr_bg <- SpatialPointsDataFrame(
    coords = final[, c(1:2)],
    data = final,
    proj4string = CRS("+init=epsg:4326")
  )

  return(list(
    bg = bg,
    pr_bg = final,
    points_pr_bg = spdf_pr_bg,
    area_M = area_M,
    isna = isna,
    na_occs = isna_occ,
    na_bg = isna_back
  ))
}
