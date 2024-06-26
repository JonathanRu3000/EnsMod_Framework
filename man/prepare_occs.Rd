% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_occs_beta.R
\name{prepare_occs}
\alias{prepare_occs}
\title{Prepare Occurrence Data of a Given Species for Habitat Suitability Modeling}
\usage{
prepare_occs(sp, valid_occurrences)
}
\arguments{
\item{sp}{A string representing the species name(s) to filter the occurrences.}

\item{valid_occurrences}{A data frame containing occurrence data, expected to include columns for species names, decimal longitude, and decimal latitude.
The columns must be named as "species", "decimallongitude", and "decimallatitude".
By default, this is set to \code{valid_occurrences} which should be defined in your environment.}
}
\value{
A list with components for statistics, cleaned occurrences data frame, and spatial points.
}
\description{
This function processes occurrence data for a given species, cleans the coordinates, and summarizes various statistics related to the occurrence data. It also generates a spatial points data frame for further spatial analysis.
Coordinates are cleaned using the \link[CoordinateCleaner]{clean_coordinates} function from the \pkg{CoordinateCleaner} package v.2 (Zizka et al., 2019).
}
\examples{
\dontrun{
species_data <- read.csv("data/processed/species/all_species.csv", dec = ".", sep = ",")
valid_occurrences <- species_data[!is.na(species_data$decimallatitude), ]
result <- \link{prepare_occs}("Panthera leo", valid_occurrences)
print(result)
}

}
\references{
Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., Farooq, H., Herdean, A., Ariza, M., Scharn, R., Svanteson, S., Wengström, N., Zizka, V., & Antonelli, A. (2019). \emph{CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases}. Methods in Ecology and Evolution, 10(5), 744-751. https://doi.org/10.1111/2041-210X.13152
}
