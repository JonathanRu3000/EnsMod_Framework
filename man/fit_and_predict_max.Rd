% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_and_predict_max_beta.R
\name{fit_and_predict_max}
\alias{fit_and_predict_max}
\title{Fit Maxent Models for Habitat Suitability Predictions}
\usage{
fit_and_predict_max(training_test, backgr.pt, rast_ext, method = NULL)
}
\arguments{
\item{training_test}{An object returned by \code{\link{test_train_aoi}} containing training and test datasets.}

\item{backgr.pt}{Background points as a dataframe with coordinates.}

\item{rast_ext}{Raster layer used to assess extrapolation risk per grid cell.}

\item{method}{Type of cross-validation approach used during model evaluation. By default, it uses 'checkerboard1' or 'jackknife' if occurrences are < 3.}
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
This function fits species distribution models using the Maxent algorithm,
and evaluates the models using the \code{\link[ENMeval]{ENMevaluate}} function from the \pkg{ENMeval} package v.2.0.4 (Kass et al. 2021).
It incorporates various partitioning methods such as 'random', 'checkerboard1', and 'bootstrap'
to assess model performance.
}
\examples{

# fit_and_predict_max(
training_test = tta,
backgr.pt = backgr$bg,
rast_ext = rast_ext,
method = NULL
)
}
\references{
Kass JM, Muscarella R, Galante PJ, Bohl CL, Pinilla‐Buitrago GE, Boria RA, Soley‐Guardia M, Anderson RP. (2021).
ENMeval 2.0: Redesigned for customizable and reproducible modeling of species' niches and distributions. Methods in Ecology and Evolution, 12(9), 1602-1608. doi: 10.1111/2041-210X.13628.
}
\seealso{
\code{\link[ENMeval]{ENMevaluate}} for more details on the evaluation methods used.
}
