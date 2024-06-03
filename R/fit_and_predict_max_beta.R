#' Fit Maxent Models for Habitat Suitability Predictions
#'
#' This function fits species distribution models using the Maxent algorithm,
#' and evaluates the models using the \code{\link[ENMeval]{ENMevaluate}} function from the \pkg{ENMeval} package v.2.0.4 (Kass et al. 2021).
#' It incorporates various partitioning methods such as 'random', 'checkerboard1', and 'bootstrap'
#' to assess model performance.
#'
#' @param training_test An object returned by \code{\link{test_train_aoi}} containing training and test datasets.
#' @param backgr.pt Background points as a dataframe with coordinates.
#' @param rast_ext Raster layer used to assess extrapolation risk per grid cell.
#' @param method Type of cross-validation approach used during model evaluation. By default, it uses 'checkerboard1' or 'jackknife' if occurrences are < 3.
#' @seealso \code{\link[ENMeval]{ENMevaluate}} for more details on the evaluation methods used.
#'
#' @references
#' Kass JM, Muscarella R, Galante PJ, Bohl CL, Pinilla‐Buitrago GE, Boria RA, Soley‐Guardia M, Anderson RP. (2021).
#' ENMeval 2.0: Redesigned for customizable and reproducible modeling of species' niches and distributions. Methods in Ecology and Evolution, 12(9), 1602-1608. doi: 10.1111/2041-210X.13628.
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
#' @import dplyr
#' @import raster
#' @import ENMeval
#' @import pROC
#'
#' @examples
#'
#' # fit_and_predict_max(
#' training_test = tta,
#' backgr.pt = backgr$bg,
#' rast_ext = rast_ext,
#' method = NULL
#' )
fit_and_predict_max <- function(training_test, backgr.pt, rast_ext, method=NULL) {
  require(dplyr)
  require(raster)
  require(ENMeval)
  require(pROC)


  data.train.xy = training_test$training$training_xy
  data.test.xy = training_test$testing$testing_xy
  data.test = training_test$testing$testing_set
  backgr.pt = backgr$bg
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
  full <- c("eval_full", "rast_pred_full")

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # prepare data ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # get only  coordinates of presences from training data
  occ_train <- data.train.xy %>%
    filter(occ == "1")

  occ_test <- data.test.xy %>%
    filter(occ == "1")

  occ_full <- data.frame(rbind(occ_test, occ_train))

  method.train <- if(is.null(method)) if(sum(occ_train$occ ==1)>=3) "checkerboard1" else "jackknife" else method
  method.comp <- if(is.null(method)) if(sum(occ_full$occ ==1)>=3) "checkerboard1" else "jackknife" else method
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # run models ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # model based on training data
  eval_max <- ENMeval::ENMevaluate(
    occ = occ_train[, c(1:2)],
    env = aoi,
    bg.coords = backgr.pt[, c(1:2)],
    RMvalues = seq(0.5, 4, 0.5),
    method = method.train,
    fc = c("L", "LQ", "LQH", "LQHP", "LQHPT"),
    parallel = F,
    numCores = parallel::detectCores() - 2,
    parallelType = "doParallel",
    algorithm = "maxnet"
  )

  # filter for optimal model
  opt.seq <- (eval_max@results %>%
                filter(or.10p.avg == min(or.10p.avg)) %>%
                filter(auc.val.avg == max(auc.val.avg)))[1,]

  # get optimal model
  mod_max <- eval.models(eval_max)[[opt.seq$tune.args]]


  # model based on full data
  eval_full_max <- ENMeval::ENMevaluate(
    occ = occ_full[, c(1:2)],
    env = aoi,
    bg.coords = backgr.pt[, c(1:2)],
    RMvalues = seq(0.5, 4, 0.5),
    method = method.comp,
    fc = c("L", "LQ", "LQH", "LQHP", "LQHPT"),
    parallel = F,
    numCores = parallel::detectCores() - 2,
    parallelType = "doParallel",
    algorithm = "maxnet"
  )

  # filter for optimal model
  opt.seq.full <- (eval_full_max@results %>%
                     filter(or.10p.avg == min(or.10p.avg)) %>%
                     filter(auc.val.avg == max(auc.val.avg)))[1,]

  # get optimal model
  mod_full_max <- eval.models(eval_full_max)[[opt.seq.full$tune.args]]

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # testing model fitted on training data####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # make predictions for testdata
  test_mod_max <- raster::predict(mod_max, data.test, type = "cloglog")

  # calculating ROC and AUC
  r_max <- roc(
    response = data.test$occ,
    predictor = as.numeric(test_mod_max)
  )

  auc_max <- auc(r_max)

  # calculating COR
  cor_max <- cor(
    as.numeric(data.test$occ),
    as.numeric(test_mod_max))

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # predict to raster ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # calculating prediction to raster based on training data
  rast_pred_max <- raster::mask(eval_max@predictions[[opt.seq$tune.args]],
                                rast_ext > 0 & rast_ext < 1,
                                maskvalue = F,
                                updatevalue = 0
  )

  # calculating prediction to raster based on full data
  rast_pred_full_max <- raster::mask(eval_full_max@predictions[[opt.seq.full$tune.args]][[1]],
                                     rast_ext > 0 & rast_ext < 1,
                                     maskvalue = F,
                                     updatevalue = 0
  )

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # create output ####
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # use template structure for max output
  max <- l

  # fill output
  for (i in 1:length(max)) {
    if (i == 1) {
      for (k in 1:length(max$eval)) {
        max$eval[[k]] <- get(paste0(eval[k], "_max"))
      }
    } else if (i == 2) {
      for (k in 1:length(max$full)) {
        max$full[[k]] <- get(paste0(full[k], "_max"))
      }
    }
  }

  dir.create(paste0(base_path,
                    "/models/maxent/", names(combine)[i]), recursive = T)

  save(eval_max, file = paste0(
    base_path,
    "/models/maxent/",
    names(combine)[i],
    "/max_model_eval_",
    names(select)[[i]],
    "_",sp_sanitized_2, ".RData"
  ))

  save(eval_full_max, file = paste0(
    base_path,
    "/models/maxent/",
    names(combine)[i],
    "/max_model_fulleval_",
    names(select)[[i]],
    "_",sp_sanitized_2, ".RData"
  ))

  return(max = max)

  rm(list = ls(pattern = "max"))
  gc()

}
