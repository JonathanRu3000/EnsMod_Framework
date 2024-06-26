% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_and_predict_brt_beta.R
\name{fit_and_predict_brt}
\alias{fit_and_predict_brt}
\title{Fit a Boosted Regression Trees (BRT) Model for Habitat Suitability Modeling}
\usage{
fit_and_predict_brt(training_test = NULL, numbs = numbs, rast_ext = NULL)
}
\arguments{
\item{training_test}{A list containing separate sub-lists for training and testing datasets along with
any additional parameters needed for the model.}

\item{numbs}{A list containing weights for model fitting on training and complete datasets.}

\item{rast_ext}{A raster object to mask the areas where model predictions are applied, enhancing the
relevance of predictions to specified geographic extents.}
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
This function fits a Boosted Regression Trees (BRT) model to ecological data for habitat suitability analysis.
It processes both training and complete datasets, performs predictions, and evaluates the model using metrics
such as ROC, AUC, and correlation (COR). Additionally, it allows for the application of model predictions to
raster layers. The methodology is based on the approach described by Valavi et al. (2022).
}
\references{
Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022).
Predictive performance of presence-only species distribution models: A benchmark study with reproducible code.
Ecological Monographs, 92(1), e1486. https://doi.org/10.1002/ecm.1486
}
