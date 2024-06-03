#' Calculate Presence, Absence Counts, and Weights
#'
#' This function calculates the number of presences and absences, as well as weights
#' to handle imbalanced data in species distribution modeling or similar contexts.
#' Weights are adjusted to balance between presence and background (absence) cases.
#'
#' The approach and parts of the code were adapted from Valavi et al., 2022.
#'
#' @param train_test A list object returned by \code{\link{test_train_aoi}}, containing training
#' and complete dataset components.
#'
#' @return A list with two elements, each a list for the training set and the complete set.
#' Each sub-list includes the number of presences (`prNum`), the number of absences (`bgNum`),
#' and a vector of weights (`casewts`) for each observation. Weights help in adjusting for
#' the imbalance between presences and absences.
#'
#' @examples
#' \dontrun{
#' # Assume `test_train_aoi` has been run previously to generate `train_test`
#' results <- calc_numbs(train_test=tta)
#' print(results)
#' }
#'
#' @references
#' Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022).
#' Predictive performance of presence-only species distribution models: A benchmark study with reproducible code.
#' Ecological Monographs, 92(1), e1486. https://doi.org/10.1002/ecm.1486
#'
#'
calc_numbs <- function(train_test){

  training_set <- train_test$training$training_set
  complete_set <- train_test$comp

  # Calculate number of presences and absences
  prNum <- sum(training_set$occ == 1) # number of presences
  bgNum <- sum(training_set$occ == 0) # number of backgrounds

  # Avoid division by zero by adding a small number to denominator if zero
  bgNum_safe <- ifelse(bgNum == 0, .Machine$double.eps, bgNum)

  # Calculate weights
  casewts <- ifelse(training_set$occ == 1, 1, prNum / bgNum_safe)

  # Calculate number of presences and absences for the complete set
  prNum_full <- sum(complete_set$occ == 1) # number of presences
  bgNum_full <- sum(complete_set$occ == 0) # number of backgrounds

  # Avoid division by zero by adding a small number to denominator if zero
  bgNum_full_safe <- ifelse(bgNum_full == 0, .Machine$double.eps, bgNum_full)

  # Calculate weights for the complete set
  casewts_full <- ifelse(complete_set$occ == 1, 1, prNum_full / bgNum_full_safe)

  # Return results in a structured list
  return(list(train = list(prNum = prNum, bgNum = bgNum, casewts = casewts),
              full = list(prNum_full = prNum_full, bgNum_full = bgNum_full, casewts_full = casewts_full)))
}

