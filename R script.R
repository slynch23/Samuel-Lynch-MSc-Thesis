# Define working directory
input_folder <- "C:/Users/lynch/Documents/Supplementary Info"
output_folder <- input_folder  # Use same folder for both


# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)

# Define output folder
output_folder <- "C:/Users/lynch/Documents/Supplementary Info"
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Load combined results
input_data <- read.csv(file.path(output_folder, "All_Catchments_Combined_Results1.csv"))
colnames(input_data) <- trimws(colnames(input_data))  # Clean column names

# Define relevant indicators and grouping variables
indicators <- c("Q95", "BFI")
grouping_vars <- c("Model", "Scenario", "Catchment", "Period")

# Extract baseline values
baseline_data <- input_data %>%
  filter(Period == "Baseline") %>%
  select(all_of(grouping_vars[-4]), all_of(indicators))

# Calculate percent changes for 2080s only
percent_changes_2080s <- input_data %>%
  filter(Period == "2080s") %>%
  left_join(baseline_data, by = grouping_vars[-4], suffix = c("", "_Baseline")) %>%
  mutate(
    Pct_Q95 = ((Q95 - Q95_Baseline) / Q95_Baseline) * 100,
    Pct_BFI = ((BFI - BFI_Baseline) / BFI_Baseline) * 100
  ) %>%
  select(Model, Scenario, Catchment, Period, Pct_Q95, Pct_BFI)

# Save percent change results
write.csv(percent_changes_2080s, file.path(output_folder, "percent_changes_Q95_BFI_2080s.csv"), row.names = FALSE)
cat("Saved percent change results to CSV.\n")

#########################################################
# 1. Combined Side-by-Side Boxplot (Q95 + BFI)
#########################################################
combined_long <- percent_changes_2080s %>%
  pivot_longer(cols = starts_with("Pct_"), names_to = "Metric", values_to = "Change") %>%
  mutate(
    Metric = recode(Metric,
                    "Pct_Q95" = "Q95 (Low Flow)",
                    "Pct_BFI" = "BFI (Base Flow Index)")
  )

gg_combined_side_by_side <- ggplot(combined_long, aes(x = Scenario, y = Change, fill = Scenario)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(
    title = "2080s Percent Change in Q95 and BFI by SSP Scenario",
    x = "SSP Scenario",
    y = "Percent Change",
    fill = "Scenario"
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 1) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

ggsave(
  file.path(output_folder, "Boxplot_2080s_Combined_Q95_BFI.pdf"),
  plot = gg_combined_side_by_side,
  width = 12, height = 6
)
cat("Saved side-by-side combined boxplot.\n")

#########################################################
# 2. Combined Horizontal Boxplot (Q95 + BFI, stacked)
#########################################################
horizontal_combined <- percent_changes_2080s %>%
  pivot_longer(cols = starts_with("Pct_"), names_to = "Metric", values_to = "Change") %>%
  mutate(
    Metric = recode(Metric,
                    "Pct_Q95" = "Q95 (Low Flow)",
                    "Pct_BFI" = "BFI (Base Flow Index)"),
    Catchment = factor(Catchment),
    Scenario = factor(Scenario)
  )

gg_horizontal_combined <- ggplot(horizontal_combined, aes(x = Change, y = Catchment, fill = Scenario)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Metric, scales = "free_x") +
  labs(
    title = "Percent Change in Low Flow Indicators by Catchment (2080s)",
    x = "Percent Change (%)",
    y = "Catchment",
    fill = "SSP"
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 1) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

ggsave(
  file.path(output_folder, "Combined_Horizontal_Boxplot_Q95_BFI_2080s.pdf"),
  plot = gg_horizontal_combined,
  width = 16, height = 12
)

cat("Saved combined horizontal boxplot for Q95 and BFI.\n")

##Cliffs Delta Analysis

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(effsize)
library(readr)
library(purrr)
library(boot)

setwd("C:/Users/lynch/Documents/Supplementary Info")

# Load filtered dataset
df_cleaned <- read.csv("percent_changes_Q95_BFI_2080s.csv") %>%
  mutate(
    Scenario = tolower(trimws(Scenario)),
    Catchment = as.character(Catchment)
  )

# Define indicators and scenario pairs
indicators <- c("Pct_Q95", "Pct_BFI")
scenario_pairs <- list(
  "SSP5 vs SSP1" = c("ssp5", "ssp1"),
  "SSP3 vs SSP1" = c("ssp3", "ssp1")
)

# Cliff’s delta function
cliff_stat <- function(data, indices) {
  d <- data[indices, ]
  tryCatch(
    cliff.delta(Pct ~ Scenario, data = d)$estimate,
    error = function(e) NA_real_
  )
}

# Main analysis function (no individual plots)
run_comparison <- function(indicator, scenarios, label, data) {
  df_indicator <- data %>%
    filter(Scenario %in% scenarios) %>%
    mutate(
      Scenario = factor(Scenario, levels = scenarios),
      Pct = as.numeric(.data[[indicator]])
    )
  
  # Wilcoxon test
  wilcox_df <- df_indicator %>%
    group_by(Catchment) %>%
    group_map(~ {
      if (n_distinct(.x$Scenario) == 2) {
        tibble(
          Catchment = .y$Catchment,
          p_value = wilcox.test(Pct ~ Scenario, data = .x)$p.value
        )
      } else {
        tibble(Catchment = .y$Catchment, p_value = NA_real_)
      }
    }) %>%
    bind_rows() %>%
    mutate(Comparison = label, Indicator = indicator)
  
  # Cliff’s Delta + bootstrapped CI
  valid_groups <- df_indicator %>%
    group_by(Catchment) %>%
    filter(n_distinct(Scenario) == 2) %>%
    group_split()
  
  cliff_bootstrap <- map_dfr(valid_groups, function(group_data) {
    boot_result <- boot(data = group_data, statistic = cliff_stat, R = 1000)
    ci <- tryCatch(boot.ci(boot_result, type = "perc"), error = function(e) NULL)
    
    tibble(
      Catchment = unique(group_data$Catchment),
      Period = "2080s",
      cliffs_delta = boot_result$t0,
      ci_low = if (!is.null(ci)) ci$percent[4] else NA_real_,
      ci_high = if (!is.null(ci)) ci$percent[5] else NA_real_,
      Indicator = indicator,
      Comparison = label
    )
  }) %>% filter(!is.na(cliffs_delta))
  
  return(list(wilcox = wilcox_df, cliff = cliff_bootstrap))
}

# Run comparisons
all_results <- map(indicators, function(ind) {
  map2(scenario_pairs, names(scenario_pairs), function(pair, label) {
    run_comparison(ind, pair, label, df_cleaned)
  })
})

# Combine results
wilcox_results <- map_dfr(all_results, ~ map_dfr(.x, "wilcox"))
cliff_results <- map_dfr(all_results, ~ map_dfr(.x, "cliff"))

# Save CSV outputs
write_csv(wilcox_results, "wilcoxon_results_2080s.csv")
write_csv(cliff_results, "cliffs_delta_results_2080s.csv")
cat("✅ Cliff’s Delta + Wilcoxon results saved.\n")

# Recode and factor for plotting
cliff_results <- cliff_results %>%
  mutate(
    Indicator = recode(Indicator,
                       "Pct_Q95" = "Q95 (Low Flow)",
                       "Pct_BFI" = "BFI (Base Flow Index)"),
    Comparison = factor(Comparison, levels = c("SSP5 vs SSP1", "SSP3 vs SSP1"))
  )

# Plot: One per Indicator, colored by Comparison
for (metric in unique(cliff_results$Indicator)) {
  plot_data <- cliff_results %>% filter(Indicator == metric)
  
  p <- ggplot(plot_data, aes(x = Catchment, y = cliffs_delta, color = Comparison)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  position = position_dodge(width = 0.5),
                  width = 0.25) +
    labs(
      title = paste("Cliff’s Delta with 95% CI (2080s) -", metric),
      x = "Catchment",
      y = "Cliff's Delta",
      color = "Comparison"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  ggsave(
    filename = paste0("Combined_Cliffs_Delta_2080s_", gsub(" ", "_", metric), ".png"),
    plot = p,
    width = 12,
    height = 6,
    bg = "white"
  )
}


###RP10 Low Flow Analysis

input_folder <- "C:/Users/lynch/Documents/Supplementary Info"


# Get all catchment files
catchment_files <- list.files(input_folder, pattern = "S\\d+_flow_gr4j\\.12gcms\\.csv", full.names = TRUE)

all_results <- list()

for (file in catchment_files) {
  catchment_name <- sub("^S(\\d+)_.*", "S\\1", basename(file))
  print(paste("Processing catchment:", catchment_name))
  
  flow_data <- read.csv(file)
  scenario_columns <- grep("^CM", colnames(flow_data), value = TRUE)
  print(paste("Found", length(scenario_columns), "scenario columns."))
  
  if (length(scenario_columns) == 0) {
    warning(paste("No scenario columns found for", catchment_name))
    next
  }
  
  try({
    # Only perform low flow analysis
    low_flow_results <- map_dfr(scenario_columns, ~ perform_moving_window(flow_data, .x, flow_type = "low"))
    print(paste("Rows in low_flow_results:", nrow(low_flow_results)))
    
    if (nrow(low_flow_results) > 0) {
      all_results[[catchment_name]] <- low_flow_results
      
      # Save individual result
      write_csv(low_flow_results, file.path(output_folder, paste0(catchment_name, "_low_flow_moving_window_results.csv")))
      
      # Combine all low flow results
      if (length(all_results) > 0) {
        all_catchments_results <- bind_rows(all_results, .id = "Catchment")
        write_csv(all_catchments_results, file.path(output_folder, "all_catchments_low_flow_moving_window_results.csv"))
        print("✅ All catchments combined and saved.")
      } else {
        warning("❌ No results found to save.")
      }

library(tidyverse)

# Load data
data <- read.csv("C:/Users/lynch/Documents/Supplementary Info/all_catchments_low_flow_moving_window_results.csv")

# --- STEP 1: Prepare RP10 low flow data ---
long_data <- data %>%
  pivot_longer(cols = starts_with("RP"), names_to = "ReturnPeriod", values_to = "Flow") %>%
  filter(FlowType == "low", ReturnPeriod == "RP10") %>%
  mutate(
    SSP = str_extract(Scenario, "ssp\\d"),
    GCM = str_extract(Scenario, "CM\\d+")
  )

# --- STEP 2: GCM medians per SSP ---
gcm_ssp_medians <- long_data %>%
  group_by(StartYear, SSP, GCM) %>%
  summarise(GCM_Median = median(Flow, na.rm = TRUE), .groups = "drop")

# --- STEP 3: SSP-wide medians ---
ssp_medians <- long_data %>%
  group_by(StartYear, SSP) %>%
  summarise(SSP_Median = median(Flow, na.rm = TRUE), .groups = "drop")

# --- STEP 4: Plot ---
rp10_plot <- ggplot() +
  geom_line(data = gcm_ssp_medians,
            aes(x = StartYear, y = GCM_Median, group = interaction(GCM, SSP), color = SSP),
            linetype = "dotted", size = 0.8, alpha = 0.6) +
  geom_line(data = ssp_medians,
            aes(x = StartYear, y = SSP_Median, color = SSP),
            size = 1.5) +
  scale_color_manual(values = c(
    "ssp1" = "#1b9e77",
    "ssp3" = "#d95f02",
    "ssp5" = "#7570b3"
  )) +
  labs(
    title = "RP10 Low Flows – Solid = SSP Median • Dotted = GCM Medians (colored by SSP)",
    x = "Start Year", y = "Flow (m³/s)", color = "SSP"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

ggsave(
  filename = "RP10_Low_Flows_SSP_Comparison.png",
  plot = rp10_plot,
  path = "C:/Users/lynch/Documents/Supplementary Info",
  width = 12,
  height = 8,
  dpi = 600,
  units = "in",
  bg = "white"  # <- Ensures white background
)


##############percent change maps
# Load libraries
library(ggplot2)
library(sf)
library(dplyr)
library(readr)
library(viridis)

# Step 1: Load the data
df <- read.csv("median_summary_2080s.csv")

# Step 2: Convert to spatial object
points_sf <- st_as_sf(df, coords = c("Easting", "Northing"), crs = 29902)

library(rnaturalearth)
library(rnaturalearthdata)

ireland <- ne_countries(scale = "medium", country = "Ireland", returnclass = "sf")
ireland <- st_transform(ireland, crs = 29902)  # match with your CRS


# Step 4: Filter NAs
points_sf_bfi <- points_sf %>% filter(!is.na(Median_Pct_BFI))

# Step 5: Plot BFI changes
ggplot() +
  geom_sf(data = ireland, fill = "grey95", color = "black") +
  geom_sf(data = points_sf_bfi, aes(color = Median_Pct_BFI), size = 3) +
  scale_color_viridis_c(
    option = "A",  # Perceptual color palette
    direction = -1,
    name = "% Change BFI"
  ) +
  facet_wrap(~Scenario, nrow = 1) +
  labs(
    title = "Median % Change in BFI (2080s) by SSP",
    x = "", y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold")
  )




# Filter NAs
points_sf_q95 <- points_sf %>% filter(!is.na(Median_Pct_Q95))

# Plot Q95 changes
ggplot() +
  geom_sf(data = ireland, fill = "grey95", color = "black") +
  geom_sf(data = points_sf_q95, aes(color = Median_Pct_Q95), size = 3) +
  scale_color_viridis_c(
    option = "A",
    direction = -1,
    name = "% Change Q95"
  ) +
  facet_wrap(~Scenario, nrow = 1) +
  labs(
    title = "Median % Change in Q95 (2080s) by SSP",
    x = "", y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold")
  )


############
#cliffs delta scores map
# Load libraries
library(readxl)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# Load Cliff's Delta data
data <- read_excel("combined_all_sheets.xlsx")

# Convert to sf: Irish Grid (EPSG:29903) -> WGS84 (EPSG:4326)
data_sf <- st_as_sf(data, coords = c("Easting", "Northing"), crs = 29903)
data_sf <- st_transform(data_sf, crs = 4326)

# Get Ireland outline
ireland <- ne_countries(scale = "medium", country = "Ireland", returnclass = "sf")
ireland <- st_transform(ireland, crs = 4326)

# Subset datasets for plotting
q95_ssp5 <- data_sf %>% filter(Indicator == "Pct_Q95", Comparison == "SSP5 vs SSP1")
q95_ssp3 <- data_sf %>% filter(Indicator == "Pct_Q95", Comparison == "SSP3 vs SSP1")
bfi_ssp5 <- data_sf %>% filter(Indicator == "Pct_BFI", Comparison == "SSP5 vs SSP1")
bfi_ssp3 <- data_sf %>% filter(Indicator == "Pct_BFI", Comparison == "SSP3 vs SSP1")

# Function to create warm-colored Cliff's Delta maps (-0.5 to -1)
make_map <- function(data, title) {
  ggplot() +
    geom_sf(data = ireland, fill = NA, color = "grey40") +
    geom_sf(data = data, aes(color = cliffs_delta), size = 3) +
    scale_color_gradientn(
      colours = c("#7f3b08", "#b35806", "#e08214", "#fdb863", "#fee0b6", "#ffffbf"),
      values = scales::rescale(seq(-1, -0.5, length.out = 6)),
      limits = c(-1, -0.5),
      oob = scales::squish,
      name = "Cliff's Δ"
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

# Create the four maps
map_q95_ssp5 <- make_map(q95_ssp5, "Pct_Q95 – SSP5 vs SSP1")
map_q95_ssp3 <- make_map(q95_ssp3, "Pct_Q95 – SSP3 vs SSP1")
map_bfi_ssp5 <- make_map(bfi_ssp5, "Pct_BFI – SSP5 vs SSP1")
map_bfi_ssp3 <- make_map(bfi_ssp3, "Pct_BFI – SSP3 vs SSP1")

# Combine into 2x2 layout
combined_map <- (map_q95_ssp5 | map_q95_ssp3) / (map_bfi_ssp5 | map_bfi_ssp3)

# Display
combined_map

# Optionally save
ggsave("Combined_Cliffs_Delta_Maps_Orange_Yellow.png", plot = combined_map,
       width = 14, height = 10, dpi = 300, bg = "white")


##########7 outlier catchment Wilcoxon Test and Analysis
# Load required library
library(dplyr)

# Load your full dataset
df <- read.csv("stations_full.csv")

# Define the 7 catchments of interest
focus_ids <- c(12001, 16011, 16008, 25001, 27002, 36019, 25006)

# Corrected assignment of focus group
df$FocusGroup <- ifelse(df$Station.No. %in% focus_ids, "Focus", "Other")

# Exclude non-descriptor columns
excluded_cols <- c("Station.No.", "Period", "Comparison", "Indicator", "cliffs_delta", "FocusGroup")
descriptor_vars <- df %>%
  select(where(is.numeric)) %>%
  select(-any_of(excluded_cols))

# Run Wilcoxon test for each descriptor
wilcox_results <- sapply(names(descriptor_vars), function(var) {
  test <- wilcox.test(df[[var]] ~ df$FocusGroup)
  c(W = test$statistic, p = test$p.value)
}) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  arrange(p)

# Print sorted results
print(wilcox_results)

library(dplyr)

# Load your dataset
df <- read.csv("stations_full.csv")

# Define focus catchments
focus_ids <- c(12001, 16011, 16008, 25001, 27002, 36019, 25006)

# Create FocusGroup label
df$FocusGroup <- ifelse(df$`Station.No.` %in% focus_ids, "Focus", "Other")

# List of descriptor columns (excluding metadata and cliffs_delta)
excluded_cols <- c("Station.No.", "Period", "Comparison", "Indicator", "cliffs_delta", "FocusGroup")
descriptor_vars <- df %>%
  select(where(is.numeric)) %>%
  select(-any_of(excluded_cols)) %>%
  names()

# Compute medians for each variable by group
medians <- df %>%
  group_by(FocusGroup) %>%
  summarise(across(all_of(descriptor_vars), median, na.rm = TRUE)) %>%
  t() %>%
  as.data.frame()

# Clean up row and column names
colnames(medians) <- medians[1, ]
medians <- medians[-1, ]
medians <- tibble::rownames_to_column(medians, var = "Variable")

print(medians)


#######Spearman Correlation Analysis for Cliffs Delta Scores
cols <- read.csv("stations_full.csv")
str(cols)

# Load required libraries
library(dplyr)
library(tidyr)

# Define exclusions
excluded_cols <- c("Station.No.", "Period", "Comparison", "Indicator", "cliffs_delta")

# Identify numeric predictor columns
numeric_cols <- cols %>%
  select(where(is.numeric)) %>%
  select(-any_of(excluded_cols)) %>%
  colnames()

# Create all Comparison + Indicator combinations of interest
combos <- cols %>%
  filter(Comparison %in% c("SSP5 vs SSP1", "SSP3 vs SSP1"),
         Indicator %in% c("Pct_Q95", "Pct_BFI")) %>%
  select(Comparison, Indicator) %>%
  distinct()

# Initialize list to store correlation results
all_correlations <- list()

# Loop through each combo and compute correlations
for(i in 1:nrow(combos)) {
  comp <- combos$Comparison[i]
  ind <- combos$Indicator[i]
  
  subset_data <- cols %>%
    filter(Comparison == comp, Indicator == ind)
  
  cor_values <- sapply(numeric_cols, function(col) {
    cor(subset_data[[col]], subset_data$cliffs_delta, method = "spearman")
  })
  
  cor_df <- data.frame(
    Variable = names(cor_values),
    Spearman_Correlation = as.numeric(cor_values),
    Comparison = comp,
    Indicator = ind
  )
  
  all_correlations[[length(all_correlations) + 1]] <- cor_df
}

# Combine and sort all results
final_correlations <- bind_rows(all_correlations) %>%
  arrange(Indicator, Comparison, desc(abs(Spearman_Correlation)))

# View the results
print(final_correlations)

####Spearman Correlation Analysis for median percentage changes
mediancol <- read.csv("median_summary_2080s.csv")
str(mediancol)

library(dplyr)

# Define target variables for correlation
indicator_vars <- c("Median_Pct_Q95", "Median_Pct_BFI")

# Split data by Scenario
split_data <- group_split(mediancol, Scenario)

# Calculate correlations per Scenario
cor_by_scenario <- lapply(split_data, function(df) {
  numeric_data <- df %>%
    select(where(is.numeric)) %>%
    select(-Catchment)
  
  result <- sapply(indicator_vars, function(ind) {
    sapply(numeric_data %>% select(-all_of(indicator_vars)), function(x) {
      cor(x, numeric_data[[ind]], method = "spearman", use = "complete.obs")
    })
  })
  
  # Add Scenario name
  scenario_name <- unique(df$Scenario)
  result_df <- as.data.frame(result)
  result_df$Variable <- rownames(result)
  result_df$Scenario <- scenario_name
  result_df
})

# Combine all results into one data frame
cor_all <- bind_rows(cor_by_scenario)

# View result
print(cor_all)


