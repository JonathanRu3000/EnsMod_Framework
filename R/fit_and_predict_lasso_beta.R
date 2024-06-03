#' Fit Lasso Model on Ecological Data
#'
#' This function fits a lasso regression model to ecological data. It prepares the data by creating quadratic terms,
#' fits the model on training data, and evaluates it on test data. It also performs predictions on a raster object.
#' The code sources external scripts from valavi/myspatial and relies on code provided by Valavi et al. 2022.
#'
#' @param training_test A list containing elements for training and testing sets.
#' @param numbs A list containing numerical parameters such as weights for training and full data sets.
#' @param rast_ext A raster layer used for masking output predictions.
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
#'taining model evaluations, predictions, and raster outputs.
#'
#' @examples
#' \dontrun{
#' # Assuming training_test, numbs, and rast_ext are predefined:
#' results <- fit_lasso(training_test, numbs, rast_ext)
#' print(results)
#' }
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom raster mask
#' @importFrom Matrix sparse.model.matrix
#' @importFrom pROC roc auc
#' @importFrom dplyr select
#' @export
fit_and_predict_lasso <- function(training_test, numbs, rast_ext) {
  # Ensure the required packages are loaded
  require(glmnet)
  require(raster)
  require(pROC)
  require(dplyr)
  require(Matrix)

  # Correct variable names based on the provided parameters
  data.train <- training_test$training$training_set
  data.test <- training_test$testing$testing_set
  data.comp <- training_test$comp
  casewts.train <- numbs$train$casewts
  casewts.full <- numbs$full$casewts_full
  aoi <- training_test$aoi

  # Source prediction helper from an external URL (mentioned for clarity and reproducibility)
  devtools::source_url("https://raw.githubusercontent.com/rvalavi/myspatial/master/R/prediction_helper.R")

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
  # prepare data ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # function to create quadratic terms for lasso and ridge
  quad_obj <- make_quadratic(data.train, cols = c(2:length(data.train)))
  quad_full <- make_quadratic(data.comp, cols = c(2:length(data.comp)))

  # now we can predict this quadratic object on the training and testing data
  # this make two columns for each covariates used in the transformation
  training_quad <- predict.make_quadratic(quad_obj, newdata = data.train)
  testing_quad <- predict.make_quadratic(quad_obj, newdata = data.test)
  full_quad <- predict.make_quadratic(quad_full, newdata = data.comp)

  # make raster stack with linear and quadratic predictions
  rast_prediction_quad <- predict.make_quadratic(quad_obj, newdata = aoi)

  # convert the data.frames to sparse matrices
  # select all quadratic (and non-quadratic) columns, except the y (occ)
  new_vars <- names(training_quad)[names(training_quad) != "occ"]
  training_sparse <- sparse.model.matrix(~ . - 1, training_quad[, new_vars])
  testing_sparse <- sparse.model.matrix(~ . - 1, testing_quad[, new_vars])
  full_sparse <- sparse.model.matrix(~ . - 1, full_quad[, new_vars])


  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # run models ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # model based on training data
  mod_lasso <- glmnet::cv.glmnet(
    x = training_sparse,
    y = data.train$occ,
    family = "binomial",
    alpha = 1,
    weights = casewts.train,
    standardize = TRUE,
    parallel = F
  )

  # model based on full data
  full_lasso <- glmnet::cv.glmnet(
    x = full_sparse,
    y = data.comp$occ,
    family = "binomial",
    alpha = 1,
    weights = casewts.full,
    standardize = TRUE,
    parallel = F
  )

  #stopCluster(my.cluster)
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # testing model fitted on training data####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # make predictions for testdata
  test_mod_lasso <- predict(mod_lasso, testing_sparse, type = "response", s = "lambda.min")

  # calculating ROC and AUC
  r_lasso <- roc(
    response = testing_quad$occ,
    predictor = as.numeric(test_mod_lasso)
  )

  auc_lasso <- auc(r_lasso)

  # calculating COR
  cor_lasso <- cor(
    as.numeric(testing_quad$occ),
    as.numeric(test_mod_lasso)
  )

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # predict to raster ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # calculating prediction to raster based on training data
  rast_pred_lasso <- raster::mask(predict_glmnet_raster(
    rast_prediction_quad, mod_lasso,
    type = "response",
    slambda = "lambda.min"
  ),
  rast_ext > 0 & rast_ext < 1,
  maskvalue = F,
  updatevalue = 0
  )

  # calculating prediction to raster based on full data
  rast_pred_full_lasso <- raster::mask(predict_glmnet_raster(
    rast_prediction_quad,
    full_lasso,
    type = "response",
    slambda = "lambda.min"
  ),
  rast_ext > 0 & rast_ext < 1,
  maskvalue = F,
  updatevalue = 0
  )

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # create output ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # use template structure for lasso output
  lasso <- l

  # fill output
  for (i in 1:length(lasso)) {
    if (i == 1) {
      for (k in 1:length(lasso$eval)) {
        lasso$eval[[k]] <- get(paste0(eval[k], "_lasso"))
      }
    } else if (i == 2) {
      for (k in 1:length(lasso$full)) {
        lasso$full[[k]] <- get(paste0(full[k], "_lasso"))
      }
    }
  }

  # return list
  return(lasso = lasso)

}
