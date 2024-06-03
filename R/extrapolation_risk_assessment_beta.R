#' Extrapolation Risk Assessment for Habitat suitability modelling
#'
#' This function assesses the extrapolation risk of predictions based on the known calibration area.
#' It sources an external function 'exdet' from https://gist.github.com/johnbaums to perform the risk calculation.
#'
#' @param pred Raster with the environmental predictors used for modelling.
#' @param backgr An object returned by \code{\link{background_pt}} with information about the calibration area (`area_M`) to be used for masking.
#' @param st A spatial object containing the coordinate reference system of the projection area.
#' @param test_train A object returned by \code{\link{test_train_aoi}}.

#' @return A raster layer with extrapolation risk scores.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- extrapolation_risk_assessment(pred = pred, backgr = backgr, st = combined_stack)
#' plot(result)
#' }
extrapolation_risk_assessment <- function(pred, backgr, st, test_train) {

  # use the fast exdet function provided by johnbaums
  devtools::source_url("https://gist.githubusercontent.com/johnbaums/a3347b6e81e7e995e5f4027aa6178b8f/raw/321f053bcaf372c9e9ed5dbca416c7738fd0f220/exdet.R")

  # Crop raster layer to the extent of the defined box
  rs <- pred


  # Convert background area to sf object and set its coordinate reference system
  M <- sf::st_as_sf(backgr$area_M)
  sf::st_crs(M) <- sf::st_crs(st)

  # Mask the raster with the sf object
  M_mask <- raster::mask(rs, M)

  # Convert raster to dataframe and handle missing values
  df_G <- raster::as.data.frame(rs, xy = TRUE)
  df_M <- raster::as.data.frame(M_mask)

  # Run the extrapolation detection on non-missing data
  exdett <- exdet(na.omit(df_M), na.omit(df_G[,-c(1:2)]))
  x <- df_G[,3]

  # Replace non-missing values in x with scores from exdett
  scores <- exdett %>%
    dplyr::pull(score)

  x[!is.na(x)] <- scores

  # Combine the coordinates with the scores
  r <- cbind(df_G[,1:2], x)

  extent <- raster::extent(test_train$aoi)
  # Convert the dataframe back to a raster object
  rast_ext <- raster::rasterFromXYZ(r) %>%
    raster::crop(., extent)

  return(rast_ext)
}
