#' Fit a Boosted Regression Trees (BRT) Model for Habitat Suitability Modeling
#'
#' This function fits a Boosted Regression Trees (BRT) model to ecological data for habitat suitability analysis.
#' It processes both training and complete datasets, performs predictions, and evaluates the model using metrics
#' such as ROC, AUC, and correlation (COR). Additionally, it allows for the application of model predictions to
#' raster layers. The methodology is based on the approach described by Valavi et al. (2022).
#'
#' @param training_test A list containing separate sub-lists for training and testing datasets along with
#'        any additional parameters needed for the model.
#' @param numbs A list containing weights for model fitting on training and complete datasets.
#' @param rast_ext A raster object to mask the areas where model predictions are applied, enhancing the
#'        relevance of predictions to specified geographic extents.
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
#' @example
#' fit_and_predict_brt(
#' training_test=tta,
#' numbs = numbs,
#' rast_ext = rast_ext
#' )
#'
#' @importFrom gbm gbm.step
#' @importFrom raster predict mask
#' @importFrom pROC roc auc
#' @importFrom stats cor
#' @export
#'
fit_and_predict_brt <- function(training_test=NULL,
                    numbs=numbs,
                    rast_ext = NULL) {

  data.train = training_test$training$training_set
  data.test = training_test$testing$testing_set
  data.comp = training_test$comp
  casewts.train = numbs$train$casewts
  casewts.full = numbs$full$casewts_full
  aoi = training_test$aoi


  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # building template structure for model output ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  l <- list(
    eval = list(
      test_model = NULL, pred = NULL, r = NULL, AUC = NULL, COR = NULL,
      raster_test = NULL
    ),
    full = list(full_model = NULL, raster_full = NULL)
  )

  eval <- c("mod", "test_mod", "r", "auc", "cor", "rast_pred")
  full <- c("full", "rast_pred_full")

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # run models ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # model based on training data
  mod_brt <- gbm.step(
    data = data.train,
    gbm.x = 2:ncol(data.train),
    gbm.y = 1,
    family = "bernoulli",
    site.weights = casewts.train,
    tree.complexity = 3,
    learning.rate = 0.001,
    n.trees = 50,
    n.folds = 5,
    max.trees = 10000
  )

  # model based on full data
  full_brt <- gbm.step(
    data = data.comp,
    gbm.x = 2:ncol(data.comp),
    gbm.y = 1,
    family = "bernoulli",
    site.weights = casewts.full,
    tree.complexity = 3,
    learning.rate = 0.001,
    n.trees = 50,
    n.folds = 5,
    max.trees = 10000
  )

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # testing model fitted on training data####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # make predictions for testdata
  test_mod_brt <- predict(mod_brt,
    data.test,
    type = "response",
    n.trees = mod_brt$gbm.call$best.trees
  )

  # calculating ROC and AUC
  r_brt <- roc(
    response = data.test$occ,
    predictor = as.numeric(test_mod_brt)
  )

  auc_brt <- auc(r_brt)

  # calculating COR
  cor_brt <- cor(
    as.numeric(data.test$occ),
    as.numeric(test_mod_brt)
  )

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # predict to raster ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # calculating prediction to raster based on training data
  rast_pred_brt <- raster::mask(raster::predict(aoi, mod_brt,
    n.trees = mod_brt$gbm.call$best.trees,
    type = "response"
  ),
  rast_ext > 0 & rast_ext < 1,
  maskvalue = F,
  updatevalue = 0
  )

  # calculating prediction to raster based on full data
  rast_pred_full_brt <- raster::mask(raster::predict(aoi, full_brt,
    n.trees = full_brt$gbm.call$best.trees,
    type = "response"
  ),
  rast_ext > 0 & rast_ext < 1,
  maskvalue = F,
  updatevalue = 0
  )

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # create output ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # use template structure for lasso output
  brt <- l

  # fill output
  for (i in 1:length(brt)) {
    if (i == 1) {
      for (k in 1:length(brt$eval)) {
        brt$eval[[k]] <- get(paste0(eval[k], "_brt"))
      }
    } else if (i == 2) {
      for (k in 1:length(brt$full)) {
        brt$full[[k]] <- get(paste0(full[k], "_brt"))
      }
    }
  }

  # return list
  return(brt = brt)
}
