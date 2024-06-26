% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/master.R
\name{model_species_distributions}
\alias{model_species_distributions}
\title{Model Species Distributions}
\usage{
model_species_distributions(
  species_data,
  path_climatic,
  path_soil,
  path_topographic,
  output_base_path
)
}
\arguments{
\item{species_data}{Data frame containing species occurrences with columns for decimal latitude and longitude.}

\item{path_climatic}{String path to the directory containing climatic data.}

\item{path_soil}{String path to the directory containing soil data.}

\item{path_topographic}{String path to the directory containing topographic data.}

\item{output_base_path}{Base path where output will be saved.}
}
\value{
Invisible NULL; this function is used for its side effects of saving output.

#' @references
Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022).
Predictive performance of presence-only species distribution models: A benchmark study with reproducible code.
Ecological Monographs, 92(1), e1486. https://doi.org/10.1002/ecm.1486
}
\description{
This function prepares and models species distribution using environmental layers and species occurrence data.
It processes data for each species and applies multiple modeling methods.
}
\examples{
model_species_distributions(occurrences, "path/to/climatic", "path/to/soil", "path/to/topographic", "output/path")
}
