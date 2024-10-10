# Main script for analyzing NFkB image data                +++++++++++++++++
# Authors: Kristina Manzhula, Kai Budde-Sagert, Henrike Rebl
# Created: 2024/05/01
# Last changed: 2024/10/10

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 0. User input (please adapt if necessary) ################################

# TODO: Change to "data" only before uploading to Zenodo
input_directory_images <- file.path("E:/", "PhD", "Data", "Paper",
                                    "2024 NFkB", "imagedata")
output_directory_data <- "results"
output_directory_plots <- "plots"

reanalyze_images <- TRUE

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# General information about microscopy images:

# Channels:
# Red:   cell membrane (pkh26 dye)
# Green: NFkB-, Keap1-, IkBa-, or Nrf2
# Blue: Nucleus (DAPI)


# # 0. Rename files
# 
# czi_files <- list.files(path = input_directory_images,
#                         pattern = "czi", full.names = TRUE)
# 
# czi_file_names <- basename(czi_files)
# czi_files_dir <- unique(dirname(czi_files))
# 
# 
# # new_names <- gsub(pattern = "2023-04-26_(.+)", replacement = "\\1", x = czi_file_names)
# # new_names <- gsub(pattern = "IKBa", replacement = "iiii", x = czi_file_names)
# # new_names <- gsub(pattern = "iiii", replacement = "IkBa", x = czi_file_names)
# # new_names <- gsub(pattern = "Argon180", replacement = "Arg180", x = czi_file_names)
# 
# # add seconds
# new_names <- gsub(pattern = "(.+_NFkB_[0-9]{1,3})_([0-9]{1,2}\\.czi)", replacement = "\\1s_\\2", x = czi_file_names)
# new_names <- gsub(pattern = "(.+_Keap1_[0-9]{1,3})_([0-9]{1,2}\\.czi)", replacement = "\\1s_\\2", x = new_names)
# new_names <- gsub(pattern = "(.+_IkBa_[0-9]{1,3})_([0-9]{1,2}\\.czi)", replacement = "\\1s_\\2", x = new_names)
# new_names <- gsub(pattern = "(.+_Nrf2_[0-9]{1,3})_([0-9]{1,2}\\.czi)", replacement = "\\1s_\\2", x = new_names)
# 
# #add date
# new_names <- gsub(pattern = "(.+)_([0-9]{1,2}\\.czi)", replacement = "\\1_20230724_\\2", x = new_names)
# 
# file.rename(from = czi_files, to = file.path(czi_files_dir, new_names))

# 1. Install and load packages #############################################

# Set groundhog day for reproducibility (see https://groundhogr.com)
groundhog.day <- "2023-01-01"

if(!any(grepl(pattern = "groundhog", x = installed.packages(), ignore.case = TRUE))){
  install.packages("groundhog")
}

# Load packages
library(groundhog)
pkgs <- c("BiocManager", "devtools", "reticulate", "tidyverse")
groundhog.library(pkgs, groundhog.day)

# # Cite packages
# options(citation.bibtex.max=999)
# for(i in 1:length(pkgs)){
#   print(citation(pkgs[i]))
# }
# rm(i)

if(!("EBImage" %in% utils::installed.packages())){
  print("Installing EBImage.")
  BiocManager::install("EBImage")
}
require(EBImage)

# Read in Python package for reading czi files
# (Users will be asked to install miniconda
# when starting for the first time)
if(! "czifile" %in% reticulate::py_list_packages()$package){
  reticulate::py_install("czifile")
}

# Install miniconda if necessary
miniconda_path <- tryCatch(
  {
    reticulate::miniconda_path()
  }, error=function(e) {
    message('Cannot acces the path.')
    return(NA)
  }
)

# if(is.na(miniconda_path())){
# devtools::install_github("hafen/rminiconda")
# rminiconda::install_miniconda(name = "my_python")
# py <- rminiconda::find_miniconda_python("my_python")
# reticulate::use_python(py, required = TRUE)
# }

# Install the R package for reading czi images
devtools::install_github("SFB-ELAINE/readCzi", ref = "v0.4.1")
require(readCzi)

# Install the R package for analysing microscopy data
devtools::install_github("https://github.com/SFB-ELAINE/cellPixels", ref = "v0.2.11")
library(cellPixels)

rm(list = c("groundhog.day", "miniconda_path", "pkgs"))

# 2. Analyze all microscopy images #########################################

if(reanalyze_images){
  
  # Read metadata ans save it in a csv file
  czi_files <- list.files(path = input_directory_images,
                          pattern = "czi", full.names = TRUE)
  
  for(i in 1:length(czi_files)){
    if(i == 1){
      df_metadata <- readCzi::readCziMetadata(input_file <- czi_files[i],
                                              save_metadata = FALSE)
    }else{
      df_dummy    <- readCzi::readCziMetadata(input_file <- czi_files[i],
                                              save_metadata = FALSE)
      df_metadata <- dplyr::bind_rows(df_metadata, df_dummy)
    }
  }
  
  if(i > 1){
    rm(df_dummy)
  }
  rm(i)
  
  split_file_name <- strsplit(x = df_metadata$fileName, split = "_")
  
  df_metadata$cell_line    <- unlist(lapply(split_file_name, `[[`, 1))
  df_metadata$protein      <- unlist(lapply(split_file_name, `[[`, 2))
  df_metadata$exp_group    <- unlist(lapply(split_file_name, `[[`, 3))
  df_metadata$exp_date     <- lubridate::as_date(unlist(lapply(split_file_name, `[[`, 4)))
  df_metadata$image_number <- as.numeric(gsub(pattern = "\\..+", replacement = "", x = unlist(lapply(split_file_name, `[[`, 5))))
  
  rm(split_file_name)
  
  df_metadata <- df_metadata %>% 
    dplyr::relocate(cell_line, protein, exp_group, exp_date, image_number, .after = fileName)
  
  dir.create(output_directory_data, showWarnings = FALSE)
  
  readr::write_csv(x = df_metadata,
                   file = file.path(output_directory_data, "df_metadata_en.csv"))
  readr::write_csv2(x = df_metadata,
                    file = file.path(output_directory_data, "df_metadata_de.csv"))
  
  
  # Analyse czi images
  df_results <- cellPixels::cellPixels(input_dir = input_directory_images,
                                       nucleus_color = "blue",
                                       protein_in_nucleus_color = "green",
                                       protein_in_cytosol_color = "red",
                                       min_nucleus_size = 400,
                                       min_cytosol_size = 600,
                                       number_size_factor = 0.2,
                                       thresh_offset_protein_in_cytosol = 0.01,
                                       add_scale_bar = FALSE)
  
  split_file_name <- strsplit(x = df_results$fileName, split = "_")
  
  df_results$cell_line    <- unlist(lapply(split_file_name, `[[`, 1))
  df_results$protein      <- unlist(lapply(split_file_name, `[[`, 2))
  df_results$exp_group    <- unlist(lapply(split_file_name, `[[`, 3))
  df_results$exp_date     <- lubridate::as_date(unlist(lapply(split_file_name, `[[`, 4)))
  df_results$image_number <- as.numeric(gsub(pattern = "\\..+", replacement = "", x = unlist(lapply(split_file_name, `[[`, 5))))
  
  rm(split_file_name)
  
  df_results <- df_results %>% 
    dplyr::relocate(cell_line, protein, exp_group, exp_date, image_number, .after = fileName)
  
  readr::write_csv(df_results,
                   file = file.path(output_directory_data, "image_analysis_summary_en.csv"))
  
  readr::write_csv2(df_results,
                    file = file.path(output_directory_data, "image_analysis_summary_de.csv"))
  
}else{
  # Import results if all images should not be read again
  df_results <- readr::read_csv(file = file.path(output_directory_data, "image_analysis_summary_en.csv"),
                                name_repair = "universal")
}

# 3. Analyze results #######################################################

# Calculate background-normalized intensities

df_results <- df_results %>% 
  dplyr::mutate(intensity_sum_nucleus_region_red_background_normalized   = intensity_sum_nucleus_region_red/intensity_mean_background_red,
                intensity_sum_nucleus_region_green_background_normalized = intensity_sum_nucleus_region_green/intensity_mean_background_green,
                intensity_sum_nucleus_region_blue_background_normalized  = intensity_sum_nucleus_region_blue/intensity_mean_background_blue)

df_results <- df_results %>% 
  dplyr::mutate(intensity_sum_cytosol_region_red_background_normalized   = intensity_sum_cytosol_region_red/intensity_mean_background_red,
                intensity_sum_cytosol_region_green_background_normalized = intensity_sum_cytosol_region_green/intensity_mean_background_green,
                intensity_sum_cytosol_region_blue_background_normalized  = intensity_sum_cytosol_region_blue/intensity_mean_background_blue)

# Calculate ratio of intensity sum in nucleus / intensity sum in cytosol
df_results <- df_results %>% 
  dplyr::mutate(intensity_sum_ratio_nuc_cyt = intensity_sum_nucleus_region_green_background_normalized/intensity_sum_cytosol_region_green_background_normalized)


# 4. Aggregate data ########################################################

# Nucleus cytosol ratios
df_aggregate_nuc_cyt <- df_results %>% 
  dplyr::group_by(cell_line, protein, exp_group, exp_date) %>% 
  dplyr::summarise(mean_nuc_cyt = mean(intensity_sum_ratio_nuc_cyt))

df_aggregate_nuc_cyt <- df_aggregate_nuc_cyt %>% 
  dplyr::group_by(cell_line, protein, exp_group) %>% 
  dplyr::summarise(mean_of_means_nuc_cyt = mean(mean_nuc_cyt),
                   sd_nuc_cyt = sd(mean_nuc_cyt),
                   n = n())


# Normalized nucleus cytosol ratios
df_aggregate_nuc_cyt_normalized <- df_results %>% 
  dplyr::group_by(cell_line, protein, exp_group, exp_date) %>% 
  dplyr::summarise(mean_nuc_cyt = mean(intensity_sum_ratio_nuc_cyt))

df_aggregate_control <- df_results %>% 
  dplyr::group_by(cell_line, protein, exp_group, exp_date) %>% 
  dplyr::filter(exp_group == "contr") %>% 
  dplyr::summarise(mean_control = mean(intensity_sum_ratio_nuc_cyt)) %>% 
  dplyr::summarise(mean_of_means_control = mean(mean_control),
                   sd_mean_control = sd(mean_control),
                   n = n())

df_aggregate_nuc_cyt_normalized <- left_join(df_aggregate_nuc_cyt_normalized,
                                             df_aggregate_control[,-3],
                                             by = c("cell_line", "protein"))

df_aggregate_nuc_cyt_normalized <- df_aggregate_nuc_cyt_normalized %>% 
  dplyr::mutate(mean_nuc_cyt_normalized = mean_nuc_cyt / mean_of_means_control)

# 4. Plot results ##########################################################

# Nrf2
df_nrf2 <- df_aggregate_nuc_cyt[df_aggregate_nuc_cyt$protein == "Nrf2" &
                                  (df_aggregate_nuc_cyt$exp_group == "contr" |
                                  df_aggregate_nuc_cyt$exp_group == "posC"),]

plot_nrf2 <- ggplot(df_nrf2, aes(x=exp_group, y=mean_of_means_nuc_cyt, color = cell_line, width = 0.5)) + 
  geom_bar(stat="identity", position=position_dodge(width = 0.6), aes(color = cell_line, fill = exp_group), linewidth=2) +
  scale_fill_manual(name = "Exp. group", values = c("white", "lightgray")) +
  scale_color_manual(name = "Cell line", values = c("blue", "orange")) +
  geom_errorbar(aes(ymin=mean_of_means_nuc_cyt , ymax=mean_of_means_nuc_cyt+sd_nuc_cyt), linewidth = 1.5, width=.2,
                position=position_dodge(0.6)) +
  labs(x="", y = "Ratio of Nrf2 (nucleus/cytosol)") +
  theme_bw(base_size = 16)

ggsave(filename = file.path(output_directory_plots, "nrf2_nuc_cyt.pdf"),
       width = 297, height = 210, units = "mm")
ggsave(filename = file.path(output_directory_plots, "nrf2_nuc_cyt.png"),
       width = 297, height = 210, units = "mm")


# Nrf2 (normalized)
df_nrf2_normalized <- df_aggregate_nuc_cyt_normalized[df_aggregate_nuc_cyt_normalized$protein == "Nrf2",]
df_nrf2_normalized$exp_group <- factor(df_nrf2_normalized$exp_group, levels = c("contr", "Arg180", "posC", "60s", "90s", "120s", "150s", "180s"))

plot_nrf2_normalized <- ggplot(df_nrf2_normalized, aes(x=exp_group , y=mean_nuc_cyt_normalized, color = cell_line, fill = exp_group)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, aes(color = cell_line), linewidth=1) + 
  geom_boxplot(linewidth=1) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_fill_manual(name = "Exp. group", values = c("white", "darkgray", "lightgray","white","white","white","white","white")) +
  scale_color_manual(name = "Cell line", values = c("blue", "orange")) +
  facet_grid(. ~ cell_line) +
  coord_cartesian(ylim = c(0.5, 2)) +
  labs(x="", y = "Ratio of Nrf2 (nucleus/cytosol)") +
  theme_bw(base_size = 16)

ggsave(filename = file.path(output_directory_plots, "nrf2_nuc_cyt_normalized.pdf"),
       width = 297, height = 210, units = "mm")
ggsave(filename = file.path(output_directory_plots, "nrf2_nuc_cyt_normalized.png"),
       width = 297, height = 210, units = "mm")


# Keap1 (normalized)
df_keap1_normalized <- df_aggregate_nuc_cyt_normalized[df_aggregate_nuc_cyt_normalized$protein == "Keap1",]
df_keap1_normalized$exp_group <- factor(df_keap1_normalized$exp_group, levels = c("contr", "Arg180", "posC", "60s", "90s", "120s", "150s", "180s"))

plot_keap1_normalized <- ggplot(df_keap1_normalized, aes(x=exp_group , y=mean_nuc_cyt_normalized, color = cell_line, fill = exp_group)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, aes(color = cell_line), linewidth=1) + 
  geom_boxplot(linewidth=1) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_fill_manual(name = "Exp. group", values = c("white", "darkgray", "lightgray","white","white","white","white","white")) +
  scale_color_manual(name = "Cell line", values = c("blue", "orange")) +
  facet_grid(. ~ cell_line) +
  coord_cartesian(ylim = c(0.5, 2)) +
  labs(x="", y = "Ratio of Keap1 (nucleus/cytosol)") +
  theme_bw(base_size = 16)

ggsave(filename = file.path(output_directory_plots, "keap1_nuc_cyt_normalized.pdf"),
       width = 297, height = 210, units = "mm")
ggsave(filename = file.path(output_directory_plots, "keap1_nuc_cyt_normalized.png"),
       width = 297, height = 210, units = "mm")


# NFkB (normalized)
df_nfkb_normalized <- df_aggregate_nuc_cyt_normalized[df_aggregate_nuc_cyt_normalized$protein == "NFkB",]
df_nfkb_normalized$exp_group <- factor(df_nfkb_normalized$exp_group, levels = c("contr", "Arg180", "posC", "60s", "90s", "120s", "150s", "180s"))

plot_nfkb_normalized <- ggplot(df_nfkb_normalized, aes(x=exp_group , y=mean_nuc_cyt_normalized, color = cell_line, fill = exp_group)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, aes(color = cell_line), linewidth=1) + 
  geom_boxplot(linewidth=1) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_fill_manual(name = "Exp. group", values = c("white", "darkgray", "lightgray","white","white","white","white","white")) +
  scale_color_manual(name = "Cell line", values = c("blue", "orange")) +
  facet_grid(. ~ cell_line) +
  coord_cartesian(ylim = c(0.5, 2)) +
  labs(x="", y = "Ratio of NFkB (nucleus/cytosol)") +
  theme_bw(base_size = 16)

ggsave(filename = file.path(output_directory_plots, "nfkb_nuc_cyt_normalized.pdf"),
       width = 297, height = 210, units = "mm")
ggsave(filename = file.path(output_directory_plots, "nfkb_nuc_cyt_normalized.png"),
       width = 297, height = 210, units = "mm")


# IkBa (normalized)
df_ikba_normalized <- df_aggregate_nuc_cyt_normalized[df_aggregate_nuc_cyt_normalized$protein == "IkBa",]
df_ikba_normalized$exp_group <- factor(df_ikba_normalized$exp_group, levels = c("contr", "Arg180", "posC", "60s", "90s", "120s", "150s", "180s"))

plot_ikba_normalized <- ggplot(df_ikba_normalized, aes(x=exp_group , y=mean_nuc_cyt_normalized, color = cell_line, fill = exp_group)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, aes(color = cell_line), linewidth=1) + 
  geom_boxplot(linewidth=1) +
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_fill_manual(name = "Exp. group", values = c("white", "darkgray", "lightgray","white","white","white","white","white")) +
  scale_color_manual(name = "Cell line", values = c("blue", "orange")) +
  facet_grid(. ~ cell_line) +
  coord_cartesian(ylim = c(0.5, 2)) +
  labs(x="", y = "Ratio of IkBa (nucleus/cytosol)") +
  theme_bw(base_size = 16)

ggsave(filename = file.path(output_directory_plots, "ikba_nuc_cyt_normalized.pdf"),
       width = 297, height = 210, units = "mm")
ggsave(filename = file.path(output_directory_plots, "ikba_nuc_cyt_normalized.png"),
       width = 297, height = 210, units = "mm")

