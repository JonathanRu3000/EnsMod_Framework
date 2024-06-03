#' Model Species Distributions
#'
#' This function prepares and models species distribution using environmental layers and species occurrence data.
#' It processes data for each species and applies multiple modeling methods.
#'
#' @param species_data Data frame containing species occurrences with columns for decimal latitude and longitude.
#' @param path_climatic String path to the directory containing climatic data.
#' @param path_soil String path to the directory containing soil data.
#' @param path_topographic String path to the directory containing topographic data.
#' @param output_base_path Base path where output will be saved.
#' @return Invisible NULL; this function is used for its side effects of saving output.
#'
#' #' @references
#' Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022).
#' Predictive performance of presence-only species distribution models: A benchmark study with reproducible code.
#' Ecological Monographs, 92(1), e1486. https://doi.org/10.1002/ecm.1486
#'
#' @export
#' @examples
#' model_species_distributions(occurrences, "path/to/climatic", "path/to/soil", "path/to/topographic", "output/path")
model_species_distributions <- function(species_data, path_climatic, path_soil, path_topographic, output_base_path) {
  # Set seed for reproducibility
  set.seed(123)

  # Load required libraries dynamically
  packages <- c(
    "dplyr", "sp", "rgeos", "geodata", "terra", "caret", "dismo", "randomForest",
    "ranger", "sf", "stringr", "reshape2", "tidyr", "ggpubr", "ggplot2",
    "Matrix", "myspatial", "pROC", "ENMeval", "rgdal", "SDMworkshop", "scales"
  )
  new_packages <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(new_packages) > 0) install.packages(new_packages)
  lapply(packages, library, character.only = TRUE, quietly = TRUE)

  # Load and process raster data
  load_raster_data <- function(path, pattern) {
    files <- list.files(path, pattern = pattern, full.names = TRUE)
    raster::stack(files)
  }

  clim <- load_raster_data(path_climatic, "bio")
  soil <- load_raster_data(path_soil, ".tif$")
  topo <- load_raster_data(path_topographic, ".tif$")
  combined_stack <- raster::stack(clim, soil, topo)


  # Validate and prepare occurrence data
  valid_occurrences <- species_data[!is.na(species_data$decimallatitude), ]
  sp_points <- SpatialPoints(valid_occurrences[, c("decimallongitude", "decimallatitude")])

  # Summary of occurrences per species
  species_summary <- valid_occurrences %>%
    group_by(speciesname) %>%
    summarise(count = n(), nas = sum(is.na(decimallatitude)), .groups = "drop")

  # Filter species with at least 3 occurrences
  filtered_species <- species_summary %>%
    filter(count >= 3) %>%
    pull(speciesname)

  # Import modeling functions and setup environment
  source_files <- c("prepare_occs_2024-05-28", "background_pt_2024-05-28", "select_variables_pergroup_bcselect_2023-03-27", "test_train_aoi_2024-05-28", "calc_numbs_2024-05-28","tukey_soil_2024-05-28", "exdet_extrapolation", "extrapolation_risk_assessment_2024-05-28", "prediction_helper_valavi_etal._2022")
  lapply(paste0("scr/modelling/helpers/", source_files, ".R"), source)
  lapply(paste0("scr/modelling/models/fit_", c("rf", "gam", "lasso", "brt", "max"), ".R"), source)

  # Helper functions
  create_dirs <- function(base_path) {
    dirs <- c("", "/data", "/models", "/stats", "/1981-2010")
    for (d in dirs) {
      dir_path <- paste0(base_path, d)
      if (!dir.exists(dir_path)) {
        dir.create(dir_path)
      }
    }
  }

  handle_error <- function(expr) {
    tryCatch(expr, error = function(e) message("Error: ", e))
  }
  names(valid_occurrences)[2] <- "species"
  sp <- filtered_species[1]
  # Main loop for each species
  for (sp in filtered_species) {
    sp_sanitized_1 <- gsub(" ", "_", sp)
    sp_sanitized_2 <- sub("\\.$", "", sp_sanitized_1)
    base_path <- paste0(output_base_path, sp_sanitized_2)

    create_dirs(base_path)

    # Prepare occurrence data
    s1 <- system.time({
      occs <- handle_error(prepare_occs(sp = sp, valid_occurrences = valid_occurrences))
    })
    if (!is.null(occs)) {
      occs$time <- s1
      handle_error(save(occs, file = paste0(base_path, "/data/occ_details_", sp_sanitized_2, ".RData")))
    }



    # Create background points
    s2 <- system.time({
      backgr <- handle_error(background_pt(data = occs, rast = clim[[1]], n_bg = 0.005, width_cal_area = 100000))
    })

    if (!is.null(backgr)) {
      backgr$time <- s2
      handle_error(save(backgr, file = paste0(base_path, "/data/bg_", sp_sanitized_2, ".RData")))
    }

    # Soil selection
    s3 <- system.time({
      soil_new <- handle_error(tukey_soil(pattern = NULL, points = backgr, soil = soil, growthform = "small_annual"))
    })

    # Variable selection
    groups <- handle_error(list(clim = clim, notclim = stack(topo, soil_new)))

    s4 <- system.time({
      select <- handle_error(select_variables_pergroup(groups = groups, data = backgr))
    })

    if (!is.null(select)) {
      select$time <- s4
      handle_error(save(select, file = paste0(base_path, "/data/var_selection_", sp_sanitized_2, ".RData")))
    }

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Building Models for each group of predictors ####
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    combine <- vector(mode = "list", length = length(groups))
    names(combine) <- names(groups)

    i=1

    for (i in 1:length(combine)) {
      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # create testing/training data and aoi (area of interest) ####
      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # first create a raster stack with the selected predictor variables
      pred <- combined_stack[[select[[i]]$vifnames]]


      # take computation time
      s5 <- system.time({
        # run script
        tta <- handle_error(test_train_aoi(
          data = backgr,
          predictors = pred,
          partition = 0.8,
          aoi_width = 100000
        ))
      })
      if (!is.null(tta)) {
        tta$time <- s5
        handle_error(save(tta, file = paste0(
          base_path, "/data/train_test_aoi_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
      }

      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # Extrapolation risk assessment ####
      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      rast_ext <- handle_error(extrapolation_risk_assessment(pred = pred, backgr = backgr, st = combined_stack, test_train = tta)
)


      handle_error(writeRaster(rast_ext, file = paste0(
        base_path, "/data/extrapolationrisk_1981-2011_", names(select)[[i]], "_", sp_sanitized_2, ".tif"
      ), overwrite = TRUE))

      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # calculate presences and absences in training and complete set ####
      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      s6 <- system.time({
        numbs <- handle_error(calc_numbs(
          tta
        ))
      })


      if (!is.null(numbs)) {
        numbs$time <- s6
        handle_error(save(numbs, file = paste0(
          base_path, "/data/numbs_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
      }

      # Fit Maxent model
      s11 <- system.time({
        max <- handle_error(fit_and_predict_max(
          training_test = tta,
          backgr.pt = backgr$bg,
          rast_ext = rast_ext,
          method = NULL
        ))
      })

      if (!is.null(max)) {
        max$time <- s11
        handle_error(save(max, file = paste0(
          base_path, "/models/max_model_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
        combine[[i]]$max$raster$raster_test <- max$eval$raster_test
        combine[[i]]$max$raster$raster_full <- max$full$raster_full
        combine[[i]]$max$predictions <- max$eval$pred
        combine[[i]]$max$eval$AUC <- max$eval$AUC
      }


      # Fit Random Forest model
      s7 <- system.time({
        rf <- handle_error(fit_and_predict_rf(
          train_test = tta,
          ntrees = 1000,
          numbs = numbs,
          rast_ext = rast_ext
        ))
      })
      if (!is.null(rf)) {
        rf$time <- s7
        handle_error(save(rf, file = paste0(
          base_path, "/models/rf_model_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
        combine[[i]]$rf$raster$raster_test <- rf$eval$raster_test
        combine[[i]]$rf$raster$raster_full <- rf$full$raster_full
        combine[[i]]$rf$predictions <- rf$eval$pred$predictions[, "1"]
        combine[[i]]$rf$eval$AUC <- rf$eval$AUC
      }

      # Fit GAM model
      s8 <- system.time({
        gam <- handle_error(fit_and_predict_gam(
          train_test = tta,
          numbs = numbs,
          rast_ext = rast_ext
        ))
      })
      if (!is.null(gam)) {
        gam$time <- s8
        handle_error(save(gam, file = paste0(
          base_path, "/models/gam_model_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
        combine[[i]]$gam$raster$raster_test <- gam$eval$raster_test
        combine[[i]]$gam$raster$raster_full <- gam$full$raster_full
        combine[[i]]$gam$predictions <- gam$eval$pred
        combine[[i]]$gam$eval$AUC <- gam$eval$AUC
      }

      # Fit Lasso model
      s9 <- system.time({
        lasso <- handle_error(fit_and_predict_lasso(
          training_test = tta,
          numbs = numbs,
          rast_ext = rast_ext
        ))
      })

      if (!is.null(lasso)) {
        lasso$time <- s9
        handle_error(save(lasso, file = paste0(
          base_path, "/models/lasso_model_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
        combine[[i]]$lasso$raster$raster_test <- lasso$eval$raster_test
        combine[[i]]$lasso$raster$raster_full <- lasso$full$raster_full
        combine[[i]]$lasso$predictions <- lasso$eval$pred
        combine[[i]]$lasso$eval$AUC <- lasso$eval$AUC
      }

      # Fit BRT model
      s10 <- system.time({
        brt <- handle_error(fit_and_predict_brt(
          training_test=tta,
          numbs = numbs,
          rast_ext = rast_ext
        ))
      })

      if (!is.null(brt)) {
        brt$time <- s10
        handle_error(save(brt, file = paste0(
          base_path, "/models/brt_model_", names(select)[[i]], "_", sp_sanitized_2, ".RData"
        )))
        combine[[i]]$brt$raster$raster_test <- brt$eval$raster_test
        combine[[i]]$brt$raster$raster_full <- brt$full$raster_full
        combine[[i]]$brt$predictions <- brt$eval$pred
        combine[[i]]$brt$eval$AUC <- brt$eval$AUC
      }

      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # ensemble model for each group of predictors ####
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ensemble <- create_ensemble_model(combine=combine, tta=tta, i=i)


      handle_error(save(ensemble, file = paste0(
        base_path, "/models/ensemble_model_", names(select)[[i]], "_", sp_sanitized_2, ".RData"

        )))

      combine[[i]]$ensemble$raster <- ensemble$rast_pred
      combine[[i]]$ensemble$predictions <- ensemble$pred


      # combine all raster predictions in a stack
      s <- stack(
        brt$full$raster_full, rf$full$raster_full, gam$full$raster_full,
        lasso$full$raster_full, max$full$raster_full, ensemble$rast_pred
      )
      names(s) <- c("brt", "rf", "gam", "lasso", "max", "ensemble")

      for (p in 1:length(s@layers)) {
        terra::writeRaster(s[[p]], paste0(
          base_path, "/1981-2010/", names(s)[p],
          names(select)[[i]], "_", sp_sanitized_2, ".tif"
        ))
      }
    }}
}

