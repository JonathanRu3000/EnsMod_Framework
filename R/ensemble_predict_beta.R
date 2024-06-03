#' Generate Ensemble Predictions
#'
#' This function combines predictions from multiple models to generate an ensemble prediction,
#' considering only models with significant predictive value (AUC > 0.5). It calculates combined
#' ROC, AUC, COR, and creates an ensemble raster prediction.
#'
#' @param combine A list containing model outputs and predictions.
#' @param i Index for the current group of predictors within the combine list.
#' @param tta Testing set object containing observed occurrences.

#'
#' @importFrom dplyr filter select
#' @importFrom raster stack calc
#' @importFrom scales rescale
#' @importFrom pROC roc auc
#' @importFrom glmnet cv.glmnet
#' @export
create_ensemble_model <- function(combine, i=i, tta) {
  ensemble <- list()

  # Extract AUC values from models and identify significant models
  auc_ens <- data.frame(
    model = c("rf", "gam", "lasso", "brt", "max"),
    auc = c(combine[[i]]$rf$eval$AUC, combine[[i]]$gam$eval$AUC, combine[[i]]$lasso$eval$AUC, combine[[i]]$brt$eval$AUC, combine[[i]]$max$eval$AUC)
  )

  filtered <- auc_ens %>% filter(auc > 0.5)
  sorted_list <- combine[[i]][order(match(names(combine[[i]]), auc_ens$model))]

  # Combine predictions from all significant models
  pred_testdata <- data.frame(
    scales::rescale(sorted_list[[1]]$predictions, to = c(0, 1)) * 1,
    scales::rescale(sorted_list[[2]]$predictions, to = c(0, 1)) * 1.2,
    scales::rescale(sorted_list[[3]]$predictions, to = c(0, 1)) * 1.4,
    scales::rescale(sorted_list[[4]]$predictions, to = c(0, 1)) * 1.6,
    scales::rescale(sorted_list[[5]]$predictions, to = c(0, 1)) * 2
  )
  names(pred_testdata) <- auc_ens$model
  pred_testdata <- pred_testdata %>% dplyr::select(all_of(filtered$model))
  mean_prob <- rowSums(pred_testdata) / 7.2

  # Calculate ROC, AUC, and COR for the ensemble model
  r_con <- roc(response = tta$testing$testing_set$occ, predictor = mean_prob)
  auc_con <- auc(r_con)
  cor_con <- cor(as.numeric(tta$testing$testing_set$occ), mean_prob)

  # Calculate ensemble raster prediction
  # calculate ensemble raster prediction
  all <- raster::stack(
    sdmvspecies::rescale(sorted_list[[1]]$raster$raster_full) * 1,
    sdmvspecies::rescale(sorted_list[[2]]$raster$raster_full) * 1.2,
    sdmvspecies::rescale(sorted_list[[3]]$raster$raster_full) * 1.4,
    sdmvspecies::rescale(sorted_list[[4]]$raster$raster_full) * 1.6,
    sdmvspecies::rescale(sorted_list[[5]]$raster$raster_full) * 2
  )

  names(all) <- auc_ens$model
  all <- all[[filtered$model]]
  ensemble$pred <- mean_prob
  ensemble$r <- r_con
  ensemble$auc <- auc_con
  ensemble$cor <- cor_con
  ensemble$rast_pred <- raster::calc(all, sum) / 7.2


  return(ensemble)
}
