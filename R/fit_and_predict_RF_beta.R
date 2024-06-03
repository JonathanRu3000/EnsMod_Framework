#' Fit Random Forest Models for Habitat suitability Predictions
#'
#' Constructs random forest models using the ranger package for both training and full data sets.
#' This function computes the AUC (Area Under the Curve) and Pearson correlation for models fitted on
#' the training data and the test data. Predictions for raster formats are also generated for both model types.
#' The methodology is based on the work provided by Valavi et al. (2022) and uses a downsampling approach to increase the
#' predictive power of the random forest model.
#'
#' @param train_test An object returned by \code{\link{test_train_aoi}}. A list containing training data, full data, and test data:
#'        $training$training_set: data frame of training data,
#'        $comp: data frame of full data,
#'        $testing$testing_set: data frame of test data,
#'        $aoi: Area of Interest (AOI) as a spatial object for raster prediction.
#' @param ntrees Integer, the number of trees in each forest. Default is 1000.
#' @param numbs An object returned by \code{\link{calc_numbs}}A list specifying numeric proportions and case weights for training and full data models:
#'        $train: list with `prNum` and `bgNum` for training model,
#'        $full: list with `prNum_full` and `bgNum_full` for full model.
#' @param rast_ext Raster object containing information of extrapolation risk returned by \code{\link{extrapolation_risk_assessment}}
#'
#' @return A structured list comprising two major components:
#'   - `$eval`: Contains the outcomes and evaluations based on the training dataset. It includes:
#'     - `test_model`: The model object returned from training.
#'     - `pred`: Predictions made by the model on the test dataset.
#'     - `roc`: The ROC curve object indicating the performance of the model at various threshold settings.
#'     - `AUC`: The area under the ROC curve, a scalar measure of model performance.
#'     - `COR`: The correlation coefficient between observed outcomes and model predictions, indicating prediction accuracy.
#'     - `raster_test`: Raster predictions based on the training model, masked by `rast_ext`.
#'   - `$full`: Contains the outcomes based on the complete dataset (training and testing combined). It includes:
#'     - `full_model`: The model object returned from training on the complete dataset.
#'     - `raster_full`: Raster predictions based on the full model, also masked by `rast_ext`.
#'
#' @references
#' Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022).
#' Predictive performance of presence-only species distribution models: A benchmark study with reproducible code.
#' Ecological Monographs, 92(1), e1486. https://doi.org/10.1002/ecm.1486
#' @export
fit_and_predict_rf <- function(train_test = NULL, ntrees = 1000, numbs = NULL, rast_ext = NULL) {

  # Unpack train_test list
  datatr <- train_test$training$training_set
  datacomp <- train_test$comp
  datatest <- train_test$testing$testing_set
  aoi <- train_test$aoi

  # Unpack numbs list
  numbstr <- numbs$train
  numbcomp <- numbs$full

  # Template structure for model output
  l <- list(
    eval = list(test_model = NULL, pred = NULL, r = NULL, AUC = NULL, COR = NULL, raster_test = NULL),
    full = list(full_model = NULL, raster_full = NULL)
  )

  # Run random forest model on training data
  mod_rf <- ranger(
    formula = as.factor(occ) ~ .,
    data = datatr,
    num.trees = ntrees,
    probability = TRUE,
    sample.fraction = numbstr$prNum / numbstr$bgNum,
    case.weights = numbstr$casewts,
    num.threads = 10
  )

  # Run random forest model on full data
  full_rf <- ranger(
    formula = as.factor(occ) ~ .,
    data = datacomp,
    num.trees = ntrees,
    probability = TRUE,
    sample.fraction = numbcomp$prNum_full / numbcomp$bgNum_full,
    case.weights = numbcomp$casewts_full,
    num.threads = 10
  )

  # Predict and evaluate on test data
  test_mod_rf <- predict(mod_rf, datatest, type = "response")
  r_rf <- roc(response = datatest$occ, predictor = as.numeric(test_mod_rf$predictions[, "1"]))
  auc_rf <- auc(r_rf)
  cor_rf <- cor(as.numeric(datatest$occ), test_mod_rf$predictions[, "1"])

  # Raster prediction for training data
  rast_pred_rf <- raster::mask(predict(aoi, mod_rf, fun = function(model, ...) predict(model, ...)$predictions[, "1"]),
                               rast_ext > 0 & rast_ext < 1, maskvalue = FALSE, updatevalue = 0)

  # Raster prediction for full data
  rast_pred_full_rf <- raster::mask(predict(aoi, full_rf, fun = function(model, ...) predict(model, ...)$predictions[, "1"]),
                                    rast_ext > 0 & rast_ext < 1, maskvalue = FALSE, updatevalue = 0)

  # Create random forest output
  rf <- l
  eval <- c("mod", "test_mod", "r", "auc", "cor", "rast_pred")
  full <- c("full", "rast_pred_full")
  for (i in 1:length(rf)) {
    if (i == 1) {
      for (k in 1:length(rf$eval)) {
        rf$eval[[k]] <- get(paste0(eval[k], "_rf"))
      }
    } else if (i == 2) {
      for (k in 1:length(rf$full)) {
        rf$full[[k]] <- get(paste0(full[k], "_rf"))
      }
    }
  }

  return(rf)
}
