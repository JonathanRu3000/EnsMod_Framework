#' Create Training and Testing Datasets with Specified Area of Interest
#'
#' This function divides the input data into training and testing sets based on a specified partition ratio.
#' It extracts predictor variables for each set and adjusts the raster layers to the defined area of interest (AOI),
#' which is calculated using a minimum bounding box based on the extent of provided data points. The AOI specifies
#' the geographic area over which predictions should be applied.
#'
#' @param data An object returned by \code{\link{background_pt}}, containing spatial points data.
#' @param predictors A RasterStack or RasterBrick object representing the environmental variables used as predictors.
#' @param partition A numeric value between 0 and 1 indicating the proportion of the dataset to be used for training.
#'                  The default is 0.8, meaning 80% of the data is used for training and the remaining 20% for testing.
#' @param aoi_width A numeric value representing the distance in meters to create a bounding box around presence
#'                  and background points for projection. This bounding box defines the Area of Interest (AOI).
#'                  For example, 100000 represents a distance of 100 kilometers.
#'
#' @return A list with the following elements:
#'         \itemize{
#'           \item \code{training}: a sublist containing both the spatial points and the extracted predictor variables for the training set, confined to the AOI.
#'           \item \code{testing}: a sublist containing both the spatial points and the extracted predictor variables for the testing set, confined to the AOI.
#'           \item \code{comp}: a combined dataset of training and testing sets, both confined within the AOI.
#'           \item \code{aoi}: the cropped area of the raster layer adjusted to the extent of the bounding box defining the AOI.
#'         }
#' @export
#'
#' @examples
#' \dontrun{
#' result <- test_train_aoi(data = my_data, predictors = my_predictors, partition = 0.8, aoi_width = 100000)
#' print(result)
#' }
#' \dontrun{
#' result <- test_train_aoi(data = my_data, predictors = my_predictors, partition = 0.8, mbox = my_box)
#' print(result)
#' }
test_train_aoi <- function(data = NULL, predictors = NULL, partition = 0.8, aoi_width = NULL) {


  ## create minimum boundingbox based on the extent of a previously defined area
  area_M <- flexsdm::calib_area(data$points_pr_bg,
                                x = "decimallongitude",
                                y = "decimallatitude",
                                method = c("bmcp", width = aoi_width),
                                crs = "+init=epsg:4326"
  )

  box <- bbox(data$points_pr_bg)

  box[1, 1] <- terra::ext(area_M)[1]
  box[1, 2] <- terra::ext(area_M)[2]
  box[2, 1] <- terra::ext(area_M)[3]
  box[2, 2] <- terra::ext(area_M)[4]
  # Create training index
  train_index <- createDataPartition(y = as.factor(data$pr_bg$occ), p = partition, list = FALSE)

  # Subset the data into training and testing sets
  training_xy <- data$pr_bg[train_index, ]
  testing_xy <- data$pr_bg[-train_index, ]

  # Extract predictor variables for training and testing sets
  training_set <- na.omit(data.frame(occ = training_xy[, 3],
                                     raster::extract(predictors, training_xy[, -3])))
  testing_set <- na.omit(data.frame(occ = testing_xy[, 3],
                                    raster::extract(predictors, testing_xy[, -3])))

  # Rename columns if there is only one layer in predictors
  if (nlayers(predictors) == 1) {
    names(training_set) <- c("occ", names(predictors))
    names(testing_set) <- c("occ", names(predictors))
  }

  # Combine training and testing sets
  comp <- rbind(training_set, testing_set)

  # Crop raster layer to the extent of the bounding box
  aoi <- predictors %>%
    terra::crop(., box)

  # Create the result list
  tta <- list(
    training = list(training_xy = training_xy, training_set = training_set),
    testing = list(testing_xy = testing_xy, testing_set = testing_set),
    comp = comp,
    aoi = aoi
  )

  return(tta)
}
