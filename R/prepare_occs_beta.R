#' Prepare Occurrence Data of a Given Species for Habitat Suitability Modeling
#'
#' This function processes occurrence data for a given species, cleans the coordinates, and summarizes various statistics related to the occurrence data. It also generates a spatial points data frame for further spatial analysis.
#' Coordinates are cleaned using the \link[CoordinateCleaner]{clean_coordinates} function from the \pkg{CoordinateCleaner} package v.2 (Zizka et al., 2019).
#'
#' @param sp A string representing the species name(s) to filter the occurrences.
#' @param valid_occurrences A data frame containing occurrence data, expected to include columns for species names, decimal longitude, and decimal latitude.
#' The columns must be named as "species", "decimallongitude", and "decimallatitude".
#' By default, this is set to `valid_occurrences` which should be defined in your environment.
#'
#' @return A list with components for statistics, cleaned occurrences data frame, and spatial points.
#'
#' @references Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., Farooq, H., Herdean, A., Ariza, M., Scharn, R., Svanteson, S., Wengström, N., Zizka, V., & Antonelli, A. (2019). \emph{CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases}. Methods in Ecology and Evolution, 10(5), 744-751. https://doi.org/10.1111/2041-210X.13152
#' @examples
#' \dontrun{
#' species_data <- read.csv("data/processed/species/all_species.csv", dec = ".", sep = ",")
#' valid_occurrences <- species_data[!is.na(species_data$decimallatitude), ]
#' result <- \link{prepare_occs}("Panthera leo", valid_occurrences)
#' print(result)
#' }
#'
#' @import dplyr
#' @import CoordinateCleaner
#' @import raster
#' @import sp
#' @export


prepare_occs <- function(sp, valid_occurrences) {
  # Ensure the valid_occurrences data frame contains the necessary columns
  required_columns <- c("species", "decimallongitude", "decimallatitude")
  if (!all(required_columns %in% names(valid_occurrences))) {
    stop("The data frame valid_occurrences must contain the columns: species, decimallongitude, and decimallatitude")
  }

  # Filter and select relevant columns
  temp_spec <- valid_occurrences %>%
    filter(species == sp) %>%
    dplyr::select(decimallongitude, decimallatitude, species)

  # Clean coordinates and add occurrence column
  temp_spec <- CoordinateCleaner::clean_coordinates(temp_spec, lon = "decimallongitude", lat = "decimallatitude", value = "clean") %>%
    mutate(occ = 1)

  # Keep only one presence point per grid cell and count occurrences per cell
  temp_spec <- temp_spec %>%
    mutate(cells = raster::extract(clim[[1]], .[, c("decimallongitude", "decimallatitude")], cellnumbers = TRUE)[, "cells"]) %>%
    add_count(cells, name = "occ_per_cell") %>%
    distinct(cells, .keep_all = TRUE)

  # Summarize statistics
  stats <- list(
    occs_total = nrow(temp_spec),
    occs_na = sum(is.na(temp_spec$decimallongitude) | is.na(temp_spec$decimallatitude)),
    occs_double = sum(temp_spec$occ_per_cell > 1),
    occs_final = nrow(temp_spec)
  )

  # Create spatial points
  spdf_clean <- SpatialPoints(
    coords = temp_spec[, c("decimallongitude", "decimallatitude")],
    proj4string = CRS(proj4string(clim[[1]]))
  )

  # Return list of results
  return(list(stats = stats, occs = temp_spec[, c("decimallongitude", "decimallatitude", "occ")], points = spdf_clean))
}
