#' Evaluate Sample Size for Environmental Representativeness
#'
#' This function assesses how representative occurrences are for the environmental conditions of a habitat.
#' It calculates the overlap between the full environmental space defined by all occurrence points and the environmental
#' space based on a proportional share of occurrences. Niche overlap is calculated using the Schoener D index.
#' A sample might be considered representative if, for example, less than 90% of the occurrences retrieve more than 90%
#' of the full environmental space based on all occurrences. This function utilizes the `ecospat` package to calculate
#' overlap between niches.
#'
#' @param data_path Character string indicating the path to the climatic, soil, and topographic data files.
#' @param base_path Character string indicating the path to save the output files.
#' @param cores Integer specifying the number of cores to leave aside for parallel processing, with a default of 2.
#' @param iter Integer specifying the number of iterations for stochastic processes, with a default of 500.
#' @return Saves the environmental completeness analysis results and plots to the specified output path.
#' @import ecospat
#' @import dplyr
#' @import ade4
#' @import ggplot2
#' @import stringr
#' @import purrr
#' @import doParallel
#' @import foreach
#' @import raster
#' @export
#' @examples
#' perform_env_completeness_analysis("output/predictions", "data/processed", "output/tables/habitat_completness", 2)

perform_env_completeness_analysis <- function(data_path, base_path, cores = 2, iter=500) {


  # Load climatic data
  clim <- raster::stack(list.files(file.path(data_path, "climatic_na/chelsa/1981-2010/"), pattern = "tif", full.names = TRUE)[-6])

  # Load soil data
  soil <- raster::stack(list.files(file.path(data_path, "physical_soil_na/"), pattern = ".tif$", full.names = TRUE))

  # Load topographic data
  topo <- raster::stack(list.files(file.path(data_path, "topographic_selected_na/"), pattern = ".tif$", full.names = TRUE))

  # Combine soil, climatic, and topographic data
  st <- raster::stack(soil, clim, topo)

  name = list.files(base_path)
  # Function to perform environmental completeness analysis
  foreach::foreach(
    name = name,
    .packages= pack
  ) %do% {
    # List files for data and models
    ldata <- list.files(file.path(base_path, name, "data"), full.names = TRUE, pattern = "RData$") %>%
      keep(~ !str_detect(., "notclim"))

    # load data
    for (p in 1:length(ldata)) {
      load(ldata[p])
    }

    # Select predictors and points
    pred <- st[[c(select$clim$vifnames, select$notclim$vifnames)]]
    points <- backgr

    testniche <- na.omit(data.frame(
      occ = points$pr_bg[, 3],
      raster::extract(pred, points$pr_bg[, 1:2])
    ))

    # Run PCA on the entire dataset
    pca.env <- dudi.pca(testniche[, -1], scannf = FALSE, nf = 2)
    scores.globclim <- pca.env$li
    scores.sp.nat <- suprow(pca.env, testniche[testniche[, 1] == 1, -1])$li
    scores.clim.nat <- suprow(pca.env, testniche[, -1])$li

    # Gridding the native niche
    grid.clim.nat <- ecospat.grid.clim.dyn(scores.globclim, scores.clim.nat, scores.sp.nat, R = 100, th.sp = 0)

    # Test number of occurrence points
    ss <- seq(5, length(scores.sp.nat$Axis1), length.out = 10)
    p <- seq_along(scores.sp.nat$Axis1)
    so <- quantile(p, probs = seq(0.1, 1, 0.02))
    ss <- so[so > 5]


    no_cores <- cores

    # Measure time using system.time
    time_overall1 <- system.time({
      # Parallelize the outer loop
      results<- mclapply(1:length(iter), function(o) {
        sapply(ss, function(i) {
          # Sampling the specific scores
          scores.sp.samp <- suprow(pca.env, testniche[sample(which(testniche[, 1] == 1), i), 2:length(testniche)])$li

          # Gridding the invasive niche
          grid.clim.samp <- ecospat.grid.clim.dyn(
            glob = scores.globclim,
            glob1 = scores.clim.nat,
            sp = scores.sp.samp,
            R = 100,
            th.sp = 0
          )

          # Calculate niche overlap
          D.overlap <- ecospat.niche.overlap(grid.clim.nat, grid.clim.samp, cor = TRUE)$D
          return(D.overlap)
        })
      }, mc.cores = no_cores)

      # Combine the results as in the original script
      final_result <- do.call(cbind, results)
    })

    # Save results
    m <- rowMeans(data.frame(final_result))
    df <- data.frame(occ = ss, D = m) %>% mutate(perc = occ / max(occ) * 100)
    write.csv(df, file.path(base_path, paste0(name, "/stats/sampsize_eval_", name, ".csv")))



# Helper function to plot environmental completeness
  p <- ggplot(df, aes(x = perc, y = D)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 0.9, color = "red", linewidth = 3) +
    scale_x_continuous(name = "Occurrence points (%)", limits = c(0, 100), breaks = seq(10, 100, 10), expand = c(0, 1)) +
    scale_y_continuous(name = "Schoener D index", breaks = seq(0.1, 1, 0.1), limits = c(0, 1)) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      plot.background = element_blank(),
      panel.background = element_blank()
    )
  ggsave(file.path(base_path, name,  paste0("/stats/sampzise_eval_", name, ".png")), plot = p, units = "mm", width = 210, height = 290)
  pdf(file.path(base_path, name,  paste0("/stats/sampzise_eval_", name, ".pdf")), width = 15, height = 10)
  print(p)
  dev.off()

  }
}
