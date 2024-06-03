#' Generalized Additive Models Using mgcv for Habitat suitability Predictions
#'
#' This function applies generalized additive models (GAM) using the `mgcv` package to training and full datasets,
#' and computes predictive models for species distribution. It supports parallel processing to improve performance.
#' The function calculates Area Under the Curve (AUC) and Pearson correlation (COR) for model evaluation,
#' and performs predictions on raster data for both training and full models. The code is based on the work provided by Valavi et al. (2022)
#'
#' @param train_test An object returned by \code{\link{test_train_aoi}}. A list containing training data, full data, and test data:
#'        $training$training_set: data frame of training data,
#'        $comp: data frame of full data,
#'        $testing$testing_set: data frame of test data,
#'        $aoi: Area of Interest (AOI) as a spatial object for raster prediction.
#' @param numbs An object returned by \code{\link{calc_numbs}}A list specifying numeric proportions and case weights for training and full data models:
#'        $train: list with `prNum` and `bgNum` for training model,
#'        $full: list with `prNum_full` and `bgNum_full` for full model.
#' @param rast_ext Raster object containing information of extrapolation risk returned by \code{\link{extrapolation_risk_assessment}}
#'
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
#'
#' @export
#' @examples
#' train_test <- list(training = list(training_set = train_data),
#'                    testing = list(testing_set = test_data),
#'                    comp = full_data, aoi = aoi_raster)
#' numbs <- list(train = list(casewts = weights_train), full = list(casewts_full = weights_full))
#' rast_ext <- extent_raster
#' results <- fit_gam(train_test, numbs, rast_ext)
#
fit_and_predcict_gam <- function(train_test = NULL, numbs = NULL, rast_ext = NULL) {

  # Unpack datasets and weights from input lists
  data.train <- train_test$training$training_set
  data.test <- train_test$testing$testing_set
  data.comp <- train_test$comp
  casewts.train <- numbs$train$casewts
  casewts.full <- numbs$full$casewts_full
  aoi <- train_test$aoi

  # Set up parallel processing
  n.cores <- parallel::detectCores() - 2 # leaving two cores aside for system stability
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)

  # Build template structure for model output
  l <- list(
    eval = list(test_model = NULL, pred = NULL, r = NULL, AUC = NULL, COR = NULL, raster_test = NULL),
    full = list(full_model = NULL, raster_full = NULL)
  )
  eval_keys <- c("mod", "test_mod", "r", "auc", "cor", "rast_pred")
  full_keys <- c("full", "rast_pred_full")

  # Prepare data and formula for GAM
  myform <- paste("occ ~", paste(paste0("s(", colnames(data.comp)[-1], ")"), collapse = " + "))
  training_set_gam <- data.frame("occ" = as.factor(data.train$occ), data.train[, -1, drop = FALSE])
  comp_set_gam <- data.frame("occ" = as.factor(data.comp$occ), data.comp[, -1, drop = FALSE])

  # Run GAM models
  mod_gam <- mgcv::bam(
    formula = as.formula(myform),
    data = training_set_gam,
    family = binomial(link = "logit"),
    weights = casewts.train,
    method = "REML",
    cluster = my.cluster
  )

  full_gam <- mgcv::bam(
    formula = as.formula(myform),
    data = comp_set_gam,
    family = binomial(link = "logit"),
    weights = casewts.full,
    method = "REML",
    cluster = my.cluster
  )

  parallel::stopCluster(my.cluster)

  # Evaluate the model using test data
  test_mod_gam <- predict(mod_gam, data.test, type = "response")
  r_gam <- roc(response = data.test$occ, predictor = as.numeric(test_mod_gam))
  auc_gam <- auc(r_gam)
  cor_gam <- cor(as.numeric(data.test$occ), test_mod_gam)

  # Predict to raster based on model outputs
  rast_pred_gam <- raster::mask(raster::predict(aoi, mod_gam, type = "response"),
                                rast_ext > 0 & rast_ext < 1, maskvalue = FALSE, updatevalue = 0)
  rast_pred_full_gam <- raster::mask(raster::predict(aoi, full_gam, type = "response"),
                                     rast_ext > 0 & rast_ext < 1, maskvalue = FALSE, updatevalue = 0)

  # Compile model outputs into the list
  gam <- l
  for (i in 1:length(gam)) {
    if (i == 1) {
      for (k in 1:length(gam$eval)) {
        gam$eval[[k]] <- get(paste0(eval_keys[k], "_gam"))
      }
    } else if (i == 2) {
      for (k in 1:length(gam$full)) {
        gam$full[[k]] <- get(paste0(full_keys[k], "_gam"))
      }
    }
  }

  return(gam)
}
