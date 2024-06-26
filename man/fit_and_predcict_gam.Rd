% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_and_predict_gam_beta.R
\name{fit_and_predcict_gam}
\alias{fit_and_predcict_gam}
\title{Generalized Additive Models Using mgcv for Habitat suitability Predictions}
\usage{
fit_and_predcict_gam(train_test = NULL, numbs = NULL, rast_ext = NULL)
}
\arguments{
\item{train_test}{An object returned by \code{\link{test_train_aoi}}. A list containing training data, full data, and test data:
$training$training_set: data frame of training data,
$comp: data frame of full data,
$testing$testing_set: data frame of test data,
$aoi: Area of Interest (AOI) as a spatial object for raster prediction.}

\item{numbs}{An object returned by \code{\link{calc_numbs}}A list specifying numeric proportions and case weights for training and full data models:
$train: list with \code{prNum} and \code{bgNum} for training model,
$full: list with \code{prNum_full} and \code{bgNum_full} for full model.}

\item{rast_ext}{Raster object containing information of extrapolation risk returned by \code{\link{extrapolation_risk_assessment}}}
}
\value{
A structured list comprising two major components:
\itemize{
\item \verb{$eval}: Contains the outcomes and evaluations based on the training dataset. It includes:
\itemize{
\item \code{test_model}: The model object returned from training.
\item \code{pred}: Predictions made by the model on the test dataset.
\item \code{roc}: The ROC curve object indicating the performance of the model at various threshold settings.
\item \code{AUC}: The area under the ROC curve, a scalar measure of model performance.
\item \code{COR}: The correlation coefficient between observed outcomes and model predictions, indicating prediction accuracy.
\item \code{raster_test}: Raster predictions based on the training model, masked by \code{rast_ext}.
}
\item \verb{$full}: Contains the outcomes based on the complete dataset (training and testing combined). It includes:
\itemize{
\item \code{full_model}: The model object returned from training on the complete dataset.
\item \code{raster_full}: Raster predictions based on the full model, also masked by \code{rast_ext}.
}
}
}
\description{
This function applies generalized additive models (GAM) using the \code{mgcv} package to training and full datasets,
and computes predictive models for species distribution. It supports parallel processing to improve performance.
The function calculates Area Under the Curve (AUC) and Pearson correlation (COR) for model evaluation,
and performs predictions on raster data for both training and full models. The code is based on the work provided by Valavi et al. (2022)
}
\examples{
train_test <- list(training = list(training_set = train_data),
                   testing = list(testing_set = test_data),
                   comp = full_data, aoi = aoi_raster)
numbs <- list(train = list(casewts = weights_train), full = list(casewts_full = weights_full))
rast_ext <- extent_raster
results <- fit_gam(train_test, numbs, rast_ext)
}
\references{
Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022).
Predictive performance of presence-only species distribution models: A benchmark study with reproducible code.
Ecological Monographs, 92(1), e1486. https://doi.org/10.1002/ecm.1486
}
