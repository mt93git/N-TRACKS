# =========================================================================
# N-TRACKS: Cohort-Wide Markov Behavior Transition Analysis
# =========================================================================
# This script reads the cohort-wide frame-by-frame tracks, maps each file
# to its biological cohort group, computes category-specific transition matrices,
# and generates directed transition network diagrams for comparative analysis.

library(dplyr)
library(tidyr)
library(readr)
library(igraph)
library(ggplot2)

# Resolve paths relative to working directory or environment variable
project_root <- Sys.getenv("NTRACKS_ROOT")
if (project_root == "") {
  current_wd <- getwd()
  if (file.exists(file.path(current_wd, "configs", "batch_manifest_golden.csv"))) {
    project_root <- current_wd
  } else if (file.exists(file.path(dirname(current_wd), "configs", "batch_manifest_golden.csv"))) {
    project_root <- dirname(current_wd)
  } else if (file.exists(file.path(dirname(dirname(current_wd)), "configs", "batch_manifest_golden.csv"))) {
    project_root <- dirname(dirname(current_wd))
  } else {
    project_root <- current_wd
  }
}

frame_candidates <- c(
  file.path(project_root, "output", "MASTER_RESULTS_FRAME_BY_FRAME.csv"),
  file.path(project_root, "outputs_reference", "pipeline2_mac_intel_outputs", "MASTER_RESULTS_FRAME_BY_FRAME.csv")
)
frame_data_path <- ""
for (fc in frame_candidates) {
  if (file.exists(fc)) {
    frame_data_path <- fc
    break
  }
}

atlas_candidates <- c(
  file.path(project_root, "src", "ntracks", "atlas", "Reference_Atlas2.rds"),
  file.path(project_root, "src", "r_integration", "data_morphomics", "Reference_Atlas2.rds"),
  file.path(project_root, "atlas", "Reference_Atlas2.rds")
)
atlas_path <- ""
for (ac in atlas_candidates) {
  if (file.exists(ac)) {
    atlas_path <- ac
    break
  }
}

out_dir <- file.path(project_root, "output")
if (!dir.exists(out_dir)) {
  out_dir <- file.path(project_root, "outputs_reference", "pipeline2_mac_intel_outputs")
}

if (frame_data_path == "" || !file.exists(frame_data_path)) {
  stop("✖ ERROR: Cohort frame-by-frame data not found in candidate paths.")
}
if (atlas_path == "" || !file.exists(atlas_path)) {
  stop("✖ ERROR: Reference Atlas RDS not found in candidate paths.")
}

message("Loading Reference Atlas and cohort frame-by-frame tracks...")
atlas <- readRDS(atlas_path)
frame_df <- read_csv(frame_data_path, show_col_types = FALSE)

# Rename ID column if necessary
if ("ID" %in% colnames(frame_df) && !("CellID" %in% colnames(frame_df))) {
  frame_df <- frame_df %>% rename(CellID = ID)
}

# --- 1. Demographic Group Mapping ---
GROUP_MAP <- c(
  "D16_20250930_Control_LPS" = "Control_LPS", "D5_20251001_Control_LPS" = "Control_LPS", "D5_20250924_Control_LPS" = "Control_LPS", "D0_20251006_Control_LPS" = "Control_LPS", "D12_250124_Control_LPS" = "Control_LPS",
  "D16_20250930_Control_LTB4" = "Control_LTB4", "D3_241211_Control_LTB4" = "Control_LTB4", "D7_241213_Control_LTB4" = "Control_LTB4", "D2_241217_Control_LTB4" = "Control_LTB4",
  "D18_250219_Control_LTB4_RvD4" = "Control_LTB4_RvD4", "D17_250214_Control_LTB4_RvD4" = "Control_LTB4_RvD4", "D0_241216_Control_LTB4_RvD4" = "Control_LTB4_RvD4", "D6_250115_Control_LTB4_RvD4" = "Control_LTB4_RvD4", "D1_241112_Control_LTB4_RvD4" = "Control_LTB4_RvD4",
  "D0_20251014_Control_PGE2" = "Control_PGE2", "D0_20251006_Control_PGE2" = "Control_PGE2",
  "D7_250211_Control_PMA" = "Control_PMA", "D5_250122_Control_PMA" = "Control_PMA",
  "D5_250122_Control_RvD4" = "Control_LTB4",
  "D0_20251008_Control_TL0259" = "Control_TL0259", "D5_20251001_Control_TL0259" = "Control_TL0259",
  "D0_20251002_Control_TL0259_LPS" = "Control_TL0259_LPS",
  "P4_250503_FluA" = "FluA",
  "P2_20250303_Hospitalized_Control_Control" = "Hospitalized_Control_Control", "P6_251103_Hospitalized_Control_Control" = "Hospitalized_Control_Control", "P12_250408_Hospitalized_Control_Control" = "Hospitalized_Control_Control", "P5_251103_Hospitalized_Control_Control" = "Hospitalized_Control_Control", "P1_250304_Hospitalized_Control_Control" = "Hospitalized_Control_Control", "P3_250403_Hospitalized_Control_Control" = "Hospitalized_Control_Control",
  "P13_20251021_Infection_Wo_Sepsis_Only" = "Infection_Wo_Sepsis_Only",
  "P9_20251015_Infection_Wo_Sepsis_Only_Control" = "Infection_Wo_Sepsis_Only_Control", "P14_250416_Infection_Wo_Sepsis_Only_Control" = "Infection_Wo_Sepsis_Only_Control",
  "P9_20251015_Infection_Wo_Sepsis_Only_TL0259" = "Infection_Wo_Sepsis_Only_TL0259",
  "D16_20250930_LPS" = "LPS", "D5_20250924_LPS" = "LPS", "D14_250127_LPS" = "LPS", "D17_250214_LPS" = "LPS", "D0_20251006_LPS" = "LPS", "D12_250124_LPS" = "LPS",
  "D3_241211_LTB4" = "LTB4", "D1_241112_LTB4" = "LTB4", "D7_241213_LTB4" = "LTB4", "D9_250110_LTB4" = "LTB4", "D2_241217_LTB4" = "LTB4",
  "D16_250128_LTB4_RvD4" = "LTB4_RvD4", "D1_241112_LTB4_RvD4" = "LTB4_RvD4", "D9_250110_LTB4_RvD4" = "LTB4_RvD4", "D2_250114_LTB4_RvD4" = "LTB4_RvD4", "D10_250123_LTB4_RvD4" = "LTB4_RvD4", "D1_241219_LTB4_RvD4" = "LTB4_RvD4", "D2_250116_LTB4_RvD4" = "LTB4_RvD4", "D3_241114_LTB4_RvD4" = "LTB4_RvD4", "D5_250122_LTB4_RvD4" = "LTB4_RvD4", "D12_250124_LTB4_RvD4" = "LTB4_RvD4", "D6_250115_LTB4_RvD4" = "LTB4_RvD4",
  "D2_20251010_PGE2" = "PGE2", "D0_20251014_PGE2" = "PGE2", "D0_20251006_PGE2" = "PGE2",
  "D16_250128_PMA" = "PMA", "D7_250211_PMA" = "PMA", "D1_250116_PMA" = "PMA", "D2_250114_PMA" = "PMA", "D5_250122_PMA" = "PMA",
  "P8_250307_Septic_Shock" = "Septic_Shock", "P11_20251017_Septic_Shock" = "Septic_Shock",
  "P22_20251013_Septic_Shock_Control" = "Septic_Shock_Control",
  "P22_20251013_Septic_Shock_TL0259" = "Septic_Shock_TL0259", "P22_20251017_Septic_Shock_TL0259" = "Septic_Shock_TL0259",
  "P18_20250929_Severe_Sepsis" = "Severe_Sepsis",
  "D2_20251010_TL0259" = "TL0259",
  "D0_20251002_TL0259_LPS" = "TL0259_LPS",
  "D0_20251014_TL0259_PGE2" = "TL0259_PGE2"
)

# Apply group mappings
frame_df <- frame_df %>%
  mutate(
    group_clean = GROUP_MAP[SampleID],
    group_clean = ifelse(is.na(group_clean), "Unassigned_Group", group_clean)
  )

# Define high-level biological categories
frame_df <- frame_df %>%
  mutate(
    category = case_when(
      # 1. Healthy untreated / control baselines
      grepl("^Control_", group_clean) ~ "Healthy_Untreated",
      
      # 2. Healthy treated (excluding recovery)
      group_clean %in% c("LPS", "LTB4", "PGE2", "PMA", "TL0259", "TL0259_LPS", "TL0259_PGE2") ~ "Healthy_Treated",
      
      # 3. Healthy treated + recovery (Resolvin D4)
      group_clean == "LTB4_RvD4" ~ "Healthy_Treated_Recovery",
      
      # 4. Clinical septic / infection patients (baselines/untreated)
      group_clean %in% c("Septic_Shock", "Severe_Sepsis", "Infection_Wo_Sepsis_Only", "Hospitalized_Control_Control", "FluA", "Septic_Shock_Control", "Infection_Wo_Sepsis_Only_Control") ~ "Clinical_Patient",
      
      # 5. Clinical septic / infection patients + treatment (recovery)
      group_clean %in% c("Septic_Shock_TL0259", "Infection_Wo_Sepsis_Only_TL0259") ~ "Clinical_Patient_Treatment_Recovery",
      
      TRUE ~ "Unassigned_Category"
    )
  )

message("Demographic distribution by category:")
print(table(frame_df$category))

# --- 2. Sliding Window Feature Calculation ---
feature_list <- atlas$feature_list
K_VEL <- 1.288988
K_SPH <- 1.004726
WINDOW_SIZE <- 5

message("Processing rolling window features cohort-wide (optimized)...")

# Split frame_df by cell to avoid extremely slow repeated filtering in R
cell_splits <- split(frame_df, paste(frame_df$SourceFile, frame_df$CellID, sep = "_"))
message(paste("   > Splitted into", length(cell_splits), "cells. Processing rolling windows..."))

processed_cells <- lapply(cell_splits, function(cell_data) {
  cell_data <- cell_data %>% arrange(Frame)
  n_frames <- nrow(cell_data)
  if (n_frames < 2) return(NULL)
  
  centroids <- as.matrix(cell_data[, c("Centroid_X", "Centroid_Y", "Centroid_Z")])
  displacements <- sqrt(rowSums(diff(centroids)^2))
  velocities <- c(0, displacements * K_VEL)
  velocities[1] <- velocities[2]
  cell_data$Velocity <- velocities
  
  cell_features <- data.frame(matrix(ncol = length(feature_list), nrow = n_frames))
  colnames(cell_features) <- feature_list
  
  for (t in 1:n_frames) {
    win_start <- max(1, t - floor(WINDOW_SIZE/2))
    win_end   <- min(n_frames, t + floor(WINDOW_SIZE/2))
    win_idx   <- win_start:win_end
    n_win     <- length(win_idx)
    win_data  <- cell_data[win_idx, ]
    win_vels  <- win_data$Velocity
    
    n_steps_win <- n_win - 1
    if (n_steps_win < 1) n_steps_win <- 1
    factor_linear <- 30 / n_steps_win
    factor_sqrt   <- sqrt(30 / n_steps_win)
    
    cell_features$TraveledDistance[t] <- sum(win_vels) * factor_linear
    cell_features$MeanVelocity[t]     <- mean(win_vels)
    cell_features$StdVelocity[t]      <- if (n_win >= 2) sd(win_vels) else 0
    cell_features$MaxVelocity[t]      <- max(win_vels)
    cell_features$MinVelocity[t]      <- min(win_vels)
    cell_features$CVvelocity[t]       <- if (cell_features$MeanVelocity[t] > 1e-6) cell_features$StdVelocity[t] / cell_features$MeanVelocity[t] else 0
    
    p_start <- as.numeric(win_data[1, c("Centroid_X", "Centroid_Y", "Centroid_Z")])
    p_end   <- as.numeric(win_data[n_win, c("Centroid_X", "Centroid_Y", "Centroid_Z")])
    lin_dist <- sqrt(sum((p_end - p_start)^2)) * K_VEL
    cell_features$LinearDistance[t]   <- lin_dist * factor_sqrt
    tort <- if (lin_dist > 1e-6) (sum(win_vels) * factor_linear) / (lin_dist * factor_sqrt) else 1.0
    cell_features$Tortuosity[t]       <- max(1.0, tort)
    
    cell_features$MeanVolume[t]       <- mean(win_data$Volume)
    cell_features$StdVolume[t]        <- if (n_win >= 2) sd(win_data$Volume) else 0
    cell_features$MeanSurfaceArea[t]  <- mean(win_data$SurfaceArea)
    cell_features$StdSurfaceArea[t]   <- if (n_win >= 2) sd(win_data$SurfaceArea) else 0
    
    sphericities <- win_data$Sphericity * K_SPH
    cell_features$MeanSphericity[t]   <- mean(sphericities)
    cell_features$StdSphericity[t]    <- if (n_win >= 2) sd(sphericities) else 0
    
    maj <- win_data$AxisLength_1
    mid <- win_data$AxisLength_2
    mn  <- win_data$AxisLength_3
    
    cell_features$MeanMajorAxis[t]    <- mean(maj)
    cell_features$StdMajorAxis[t]     <- if (n_win >= 2) sd(maj) else 0
    cell_features$MeanMidAxis[t]      <- mean(mid)
    cell_features$StdMidAxis[t]       <- if (n_win >= 2) sd(mid) else 0
    cell_features$MeanMinAxis[t]      <- mean(mn)
    cell_features$StdMinAxis[t]       <- if (n_win >= 2) sd(mn) else 0
    
    ar1 <- maj / (mid + 1e-6)
    ar2 <- maj / (mn + 1e-6)
    ar3 <- mid / (mn + 1e-6)
    prol   <- (maj - mid) / (maj + 1e-6)
    oblate <- (mid - mn) / (mid + 1e-6)
    
    cell_features$MeanAR1[t]           <- mean(ar1)
    cell_features$StdAR1[t]            <- 1e-5
    cell_features$MeanAR2[t]           <- mean(ar2)
    cell_features$StdAR2[t]            <- 1e-5
    cell_features$MeanAR3[t]           <- mean(ar3)
    cell_features$StdAR3[t]            <- 1e-5
    cell_features$MeanProlateEllipt[t] <- mean(prol)
    cell_features$StdProlateEllipt[t]  <- 1e-5
    cell_features$MeanOblateEllipt[t]  <- mean(oblate)
    cell_features$StdOblateEllipt[t]   <- 1e-5
    
    if (n_win >= 2) {
      dx_win <- diff(win_data$Centroid_X)
      dy_win <- diff(win_data$Centroid_Y)
      angles <- atan2(dy_win, dx_win) * 180 / pi
      cell_features$AngleOfDirection_Sum[t] <- sum(abs(angles)) * factor_linear
      cell_features$AngleOfDirection_Std[t] <- sd(angles)
    } else {
      cell_features$AngleOfDirection_Sum[t] <- 0
      cell_features$AngleOfDirection_Std[t] <- 0
    }
  }
  
  cell_features$SourceFile <- cell_data$SourceFile[1]
  cell_features$CellID     <- cell_data$CellID[1]
  cell_features$Frame      <- cell_data$Frame
  cell_features$category   <- cell_data$category[1]
  return(cell_features)
})

# Filter out NULLs
processed_cells <- processed_cells[!sapply(processed_cells, is.null)]
all_feats <- bind_rows(processed_cells)

scaled_feats <- scale(all_feats[, feature_list], center = atlas$scaling_parameters$center, scale = atlas$scaling_parameters$scale)
scaled_feats[is.na(scaled_feats) | is.nan(scaled_feats) | is.infinite(scaled_feats)] <- 0

# --- 3. Projection Mapping ---
message("Projecting cells into behavior space...")
centroids <- atlas$kmeans_model$centers
predict_cluster <- function(row) {
  distances <- apply(centroids, 1, function(c) sqrt(sum((row - c)^2)))
  return(which.min(distances))
}
all_feats$Cluster <- apply(scaled_feats, 1, predict_cluster)
behavior_labels <- c("B2", "B4", "B1", "B3", "B5")
all_feats$Behavior <- behavior_labels[all_feats$Cluster]

# --- 4. Transition Hysteresis Filter ---
message("Smoothing state sequences (optimized)...")
cell_list <- split(all_feats, paste(all_feats$SourceFile, all_feats$CellID, sep = "_"))
smoothed_list <- lapply(cell_list, function(cell_data) {
  cell_data <- cell_data %>% arrange(Frame)
  states <- cell_data$Behavior
  n <- length(states)
  if (n >= 4) {
    for (i in 2:(n - 2)) {
      curr <- states[i]
      prev <- states[i - 1]
      if (curr != prev) {
        if (states[i + 1] != curr || states[i + 2] != curr) {
          states[i] <- prev
        }
      }
    }
  }
  cell_data$SmoothedBehavior <- states
  return(cell_data)
})
all_feats <- bind_rows(smoothed_list)

# Distinct list of successfully processed cells for demographic iteration
unique_cells <- all_feats %>% distinct(SourceFile, CellID, category)

# --- 5. Category-Specific Markov Transition Analyses ---
categories <- c(
  "Healthy_Untreated",
  "Healthy_Treated",
  "Healthy_Treated_Recovery",
  "Clinical_Patient",
  "Clinical_Patient_Treatment_Recovery"
)

states_ordered <- c("B1", "B2", "B3", "B4", "B5")
official_colors <- c(B1 = "#E8E654", B2 = "#5EC468", B3 = "#C780CF", B4 = "#DB5390", B5 = "#5DACF5")

darken_color <- function(color_hex, factor = 0.6) {
  rgb_vals <- col2rgb(color_hex)
  darker_rgb <- rgb_vals * factor
  darker_rgb[darker_rgb > 255] <- 255
  return(rgb(darker_rgb[1,], darker_rgb[2,], darker_rgb[3,], maxColorValue = 255))
}

save_heatmap <- function(prob_matrix, cat_name, out_dir) {
  # Convert matrix to long format dataframe for ggplot2
  df <- as.data.frame(as.table(prob_matrix))
  colnames(df) <- c("From", "To", "Probability")
  
  # Convert Probability to percentage labels
  df$Label <- sapply(df$Probability, function(w) {
    pct <- w * 100
    if (pct %% 1 == 0) sprintf("%.0f%%", pct) else sprintf("%.1f%%", pct)
  })
  
  p <- ggplot(df, aes(x = To, y = factor(From, levels = rev(states_ordered)), fill = Probability)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = Label), color = ifelse(df$Probability > 0.4, "white", "#1e293b"), fontface = "bold", size = 4) +
    scale_fill_gradientn(colors = c("#f8fafc", "#e2e8f0", "#3b82f6", "#1d4ed8"), limits = c(0, 1)) +
    labs(title = paste("Transition Heatmap:", cat_name), x = "To Behavior", y = "From Behavior") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(face = "bold", size = 9, color = "#475569"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  # Save PNG
  ggsave(file.path(out_dir, paste0("transition_heatmap_", cat_name, ".png")), plot = p, width = 6, height = 5.2, dpi = 150)
  # Save PDF
  ggsave(file.path(out_dir, paste0("transition_heatmap_", cat_name, ".pdf")), plot = p, width = 6, height = 5.2)
}

for (cat_name in categories) {
  cat_feats <- all_feats %>% filter(category == cat_name)
  cat_cells <- unique_cells %>% filter(category == cat_name)
  
  if (nrow(cat_feats) == 0) {
    message("   > Category '", cat_name, "' is empty in the current output. Skipping.")
    next
  }
  
  message("\nProcessing Category: ", cat_name)
  
  # Calculate transition counts
  trans_matrix <- matrix(0, nrow = 5, ncol = 5, dimnames = list(states_ordered, states_ordered))
  for (idx in 1:nrow(cat_cells)) {
    sf <- cat_cells$SourceFile[idx]
    cid <- cat_cells$CellID[idx]
    
    cell_data <- cat_feats %>%
      filter(SourceFile == sf, CellID == cid) %>%
      arrange(Frame)
    
    states <- cell_data$SmoothedBehavior
    n <- length(states)
    if (n < 2) next
    for (i in 1:(n - 1)) {
      s_from <- states[i]
      s_to   <- states[i + 1]
      if (s_from %in% states_ordered && s_to %in% states_ordered) {
        trans_matrix[s_from, s_to] <- trans_matrix[s_from, s_to] + 1
      }
    }
  }
  
  row_sums <- rowSums(trans_matrix)
  prob_matrix <- trans_matrix / ifelse(row_sums == 0, 1, row_sums)
  
  write.csv(prob_matrix, file.path(out_dir, paste0("transition_matrix_", cat_name, ".csv")))
  print(prob_matrix)
  save_heatmap(prob_matrix, cat_name, out_dir)
  
  # Plot directed network of Markov state transitions
  edge_list <- data.frame(from = character(), to = character(), weight = numeric(), is_loop = logical(), stringsAsFactors = FALSE)
  for (r in states_ordered) {
    for (c in states_ordered) {
      weight <- prob_matrix[r, c]
      if (r == c) {
        if (weight > 0) {
          edge_list <- rbind(edge_list, data.frame(from = r, to = c, weight = weight, is_loop = TRUE, stringsAsFactors = FALSE))
        }
      } else if (weight > 0.005) { # threshold for transitions
          edge_list <- rbind(edge_list, data.frame(from = r, to = c, weight = weight, is_loop = FALSE, stringsAsFactors = FALSE))
      }
    }
  }
  
  if (nrow(edge_list) > 0) {
    g <- graph_from_data_frame(edge_list, directed = TRUE)
    v_names <- V(g)$name
    vertex_colors <- official_colors[v_names]
    
    # Get edge source and target node names
    edge_ends <- ends(g, E(g), names = TRUE)
    from_nodes <- edge_ends[, 1]
    to_nodes   <- edge_ends[, 2]
    
    # Set edge styles:
    # 1. Colors match target node's behavior color
    E(g)$color <- official_colors[to_nodes]
    E(g)$width <- ifelse(E(g)$is_loop, 1.5, E(g)$weight * 12 + 1)
    E(g)$arrow.size <- ifelse(E(g)$is_loop, 0.4, 0.5) # keep loop arrows small, regular arrows standard
    
    # 2. Text color of edge labels matches the edge/target node color (darkened for readability)
    E(g)$label.color <- sapply(official_colors[to_nodes], darken_color, factor = 0.7)
    
    # Format labels: remove decimal if integer (e.g. 88% instead of 88.0%)
    E(g)$label <- sapply(E(g)$weight, function(w) {
      pct <- w * 100
      if (pct %% 1 == 0) sprintf("%.0f%%", pct) else sprintf("%.1f%%", pct)
    })
    
    # Shift bidirectional labels to prevent overlap (position closer to target node)
    E(g)$label.pos <- 0.5
    for (i in 1:ecount(g)) {
      u <- from_nodes[i]
      v <- to_nodes[i]
      if (u != v) {
        opposite_idx <- which(from_nodes == v & to_nodes == u)
        if (length(opposite_idx) > 0) {
          # Shift label towards target node (0.65) to put it next to arrowhead
          E(g)$label.pos[i] <- 0.65
        }
      }
    }
    
    # Position vertices in a clean pentagon layout starting from the top and going clockwise:
    angles_v <- seq(pi/2, pi/2 - 2*pi, length.out = 6)[1:5]
    names(angles_v) <- states_ordered
    
    # Align self-loop angles radially outwards so they don't overlay node text
    E(g)$loop.angle <- angles_v[from_nodes]
    
    lay <- cbind(cos(angles_v), sin(angles_v))
    rownames(lay) <- states_ordered
    lay_filtered <- lay[V(g)$name, , drop = FALSE]
    
    # Draw PNG Plot (using 2-pass overlay to draw arrows on the first plan)
    plot_path <- file.path(out_dir, paste0("transition_plot_", cat_name, ".png"))
    png(plot_path, width = 800, height = 800, res = 150)
    
    par(mar = c(2, 2, 2, 2))
    # Pass 1: Draw vertices and labels only (edges invisible)
    plot(g,
         layout = lay_filtered,
         vertex.color = vertex_colors,
         vertex.frame.color = "white",
         vertex.label.color = "white",
         vertex.label.font = 2,
         vertex.label.cex = 0.8,
         vertex.size = 30,
         edge.color = "transparent",
         edge.label = NA,
         xlim = c(-1.5, 1.5), # expanded margins to prevent loop and number clipping
         ylim = c(-1.5, 1.5),
         main = paste("Transitions:", cat_name))
         
    # Pass 2: Draw edges, arrowheads, and percentage labels on top
    plot(g,
         layout = lay_filtered,
         vertex.color = "transparent",
         vertex.frame.color = "transparent",
         vertex.label = NA,
         vertex.size = 30, # ensure edges stop exactly at the node boundary
         edge.color = E(g)$color,
         edge.width = E(g)$width,
         edge.arrow.size = E(g)$arrow.size,
         edge.arrow.width = 1.0,
         edge.label = E(g)$label,
         edge.label.color = E(g)$label.color,
         edge.label.font = 2,
         edge.label.cex = 0.7,
         edge.curved = 0.08, # curve slightly to separate bidirectional arrows
         edge.loop.radius = 0.25,
         edge.label.position = E(g)$label.pos,
         xlim = c(-1.5, 1.5),
         ylim = c(-1.5, 1.5),
         add = TRUE)
    dev.off()
    message("   > Saved PNG plot to: ", plot_path)
    
    # Draw PDF Plot (using same 2-pass overlay)
    plot_path_pdf <- file.path(out_dir, paste0("transition_plot_", cat_name, ".pdf"))
    pdf(plot_path_pdf, width = 6, height = 6)
    
    par(mar = c(2, 2, 2, 2))
    plot(g,
         layout = lay_filtered,
         vertex.color = vertex_colors,
         vertex.frame.color = "white",
         vertex.label.color = "white",
         vertex.label.font = 2,
         vertex.label.cex = 0.8,
         vertex.size = 30,
         edge.color = "transparent",
         edge.label = NA,
         xlim = c(-1.5, 1.5),
         ylim = c(-1.5, 1.5),
         main = paste("Transitions:", cat_name))
         
    plot(g,
         layout = lay_filtered,
         vertex.color = "transparent",
         vertex.frame.color = "transparent",
         vertex.label = NA,
         vertex.size = 30, # ensure edges stop exactly at the node boundary
         edge.color = E(g)$color,
         edge.width = E(g)$width,
         edge.arrow.size = E(g)$arrow.size,
         edge.arrow.width = 1.0,
         edge.label = E(g)$label,
         edge.label.color = E(g)$label.color,
         edge.label.font = 2,
         edge.label.cex = 0.7,
         edge.curved = 0.08,
         edge.loop.radius = 0.25,
         edge.label.position = E(g)$label.pos,
         xlim = c(-1.5, 1.5),
         ylim = c(-1.5, 1.5),
         add = TRUE)
    dev.off()
    message("   > Saved PDF plot to: ", plot_path_pdf)
  }
}

message("\n✔ SUCCESS: Cohort-wide category transition analysis complete.")
# =========================================================================
