# =========================================================================
# N-TRACKS: Markov Behavior Transition Analysis on Real Data
# =========================================================================

library(dplyr)
library(tidyr)
library(readr)
library(igraph)

# Define file paths (Dynamic & Cluster-Portable)
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
output_matrix_path <- file.path(out_dir, "markov_transition_matrix.csv")
output_plot_path <- file.path(out_dir, "markov_transition_plot.png")

cat("Loading reference atlas and frame-by-frame data...\n")
if (frame_data_path == "" || !file.exists(frame_data_path)) {
  stop("Error: Frame-by-frame tracking data not found in candidate paths.")
}
if (atlas_path == "" || !file.exists(atlas_path)) {
  stop("Error: Reference atlas not found in candidate paths.")
}

atlas <- readRDS(atlas_path)
frame_df <- read_csv(frame_data_path, show_col_types = FALSE)

# Rename ID column to CellID if necessary
if ("ID" %in% colnames(frame_df) && !("CellID" %in% colnames(frame_df))) {
  frame_df <- frame_df %>% rename(CellID = ID)
}

cat("Loaded", nrow(frame_df), "frame entries across", length(unique(frame_df$CellID)), "unique cell trajectories.\n")

# Calibration coefficients
K_VEL <- 1.288988
K_SPH <- 1.004726
WINDOW_SIZE <- 5  # Sliding window size (2 frames before, 2 frames after)

# Feature list used for K-means in the atlas
feature_list <- atlas$feature_list

cat("Computing windowed features for each cell frame...\n")

# Process each cell trajectory
processed_cells <- list()
unique_cells <- frame_df %>% distinct(SourceFile, CellID)

for (idx in 1:nrow(unique_cells)) {
  sf <- unique_cells$SourceFile[idx]
  cid <- unique_cells$CellID[idx]
  
  cell_data <- frame_df %>%
    filter(SourceFile == sf, CellID == cid) %>%
    arrange(Frame)
  
  n_frames <- nrow(cell_data)
  if (n_frames < 2) next
  
  # Calculate instantaneous velocities
  centroids <- as.matrix(cell_data[, c("Centroid_X", "Centroid_Y", "Centroid_Z")])
  displacements <- sqrt(rowSums(diff(centroids)^2))
  velocities <- c(0, displacements * K_VEL)
  # Smooth the initial velocity (set frame 1 velocity equal to frame 2)
  if (n_frames >= 2) {
    velocities[1] <- velocities[2]
  }
  
  cell_data$Velocity <- velocities
  
  # Compute rolling window features
  cell_features <- data.frame(matrix(ncol = length(feature_list), nrow = n_frames))
  colnames(cell_features) <- feature_list
  
  for (t in 1:n_frames) {
    # Symmetric sliding window centered at t, bounded by [1, n_frames]
    win_start <- max(1, t - floor(WINDOW_SIZE/2))
    win_end <- min(n_frames, t + floor(WINDOW_SIZE/2))
    win_idx <- win_start:win_end
    n_win <- length(win_idx)
    
    win_data <- cell_data[win_idx, ]
    win_vels <- win_data$Velocity
    
    # Extrapolation factor for 30 steps (typical full trajectory has 31 frames/30 steps)
    n_steps_win <- n_win - 1
    if (n_steps_win < 1) n_steps_win <- 1
    factor_linear <- 30 / n_steps_win
    factor_sqrt <- sqrt(30 / n_steps_win)
    
    # 1. Velocities
    cell_features$TraveledDistance[t] <- sum(win_vels) * factor_linear
    cell_features$MeanVelocity[t] <- mean(win_vels)
    cell_features$StdVelocity[t] <- if (n_win >= 2) sd(win_vels) else 0
    cell_features$MaxVelocity[t] <- max(win_vels)
    cell_features$MinVelocity[t] <- min(win_vels)
    cell_features$CVvelocity[t] <- if (cell_features$MeanVelocity[t] > 1e-6) cell_features$StdVelocity[t] / cell_features$MeanVelocity[t] else 0
    
    # 2. Linear displacement and Tortuosity
    p_start <- as.numeric(win_data[1, c("Centroid_X", "Centroid_Y", "Centroid_Z")])
    p_end <- as.numeric(win_data[n_win, c("Centroid_X", "Centroid_Y", "Centroid_Z")])
    lin_dist <- sqrt(sum((p_end - p_start)^2)) * K_VEL
    cell_features$LinearDistance[t] <- lin_dist * factor_sqrt
    tort <- if (lin_dist > 1e-6) (sum(win_vels) * factor_linear) / (lin_dist * factor_sqrt) else 1.0
    cell_features$Tortuosity[t] <- max(1.0, tort)
    
    # 3. Shape features
    cell_features$MeanVolume[t] <- mean(win_data$Volume)
    cell_features$StdVolume[t] <- if (n_win >= 2) sd(win_data$Volume) else 0
    
    cell_features$MeanSurfaceArea[t] <- mean(win_data$SurfaceArea)
    cell_features$StdSurfaceArea[t] <- if (n_win >= 2) sd(win_data$SurfaceArea) else 0
    
    sphericities <- win_data$Sphericity * K_SPH
    cell_features$MeanSphericity[t] <- mean(sphericities)
    cell_features$StdSphericity[t] <- if (n_win >= 2) sd(sphericities) else 0
    
    maj <- win_data$AxisLength_1
    mid <- win_data$AxisLength_2
    mn <- win_data$AxisLength_3
    
    cell_features$MeanMajorAxis[t] <- mean(maj)
    cell_features$StdMajorAxis[t] <- if (n_win >= 2) sd(maj) else 0
    
    cell_features$MeanMidAxis[t] <- mean(mid)
    cell_features$StdMidAxis[t] <- if (n_win >= 2) sd(mid) else 0
    
    cell_features$MeanMinAxis[t] <- mean(mn)
    cell_features$StdMinAxis[t] <- if (n_win >= 2) sd(mn) else 0
    
    # Aspect ratios and Ellipticities
    ar1 <- maj / (mid + 1e-6)
    ar2 <- maj / (mn + 1e-6)
    ar3 <- mid / (mn + 1e-6)
    prol <- (maj - mid) / (maj + 1e-6)
    oblate <- (mid - mn) / (mid + 1e-6)
    
    cell_features$MeanAR1[t] <- mean(ar1)
    cell_features$StdAR1[t] <- 1e-5  # Dummy to match Reference Atlas training data
    
    cell_features$MeanAR2[t] <- mean(ar2)
    cell_features$StdAR2[t] <- 1e-5  # Dummy to match Reference Atlas training data
    
    cell_features$MeanAR3[t] <- mean(ar3)
    cell_features$StdAR3[t] <- 1e-5  # Dummy to match Reference Atlas training data
    
    cell_features$MeanProlateEllipt[t] <- mean(prol)
    cell_features$StdProlateEllipt[t] <- 1e-5  # Dummy to match Reference Atlas training data
    
    cell_features$MeanOblateEllipt[t] <- mean(oblate)
    cell_features$StdOblateEllipt[t] <- 1e-5  # Dummy to match Reference Atlas training data
    
    # 4. Directional angles
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
  
  # Append metadata to features
  cell_features$SourceFile <- sf
  cell_features$CellID <- cid
  cell_features$Frame <- cell_data$Frame
  
  processed_cells[[length(processed_cells) + 1]] <- cell_features
}

all_features_df <- bind_rows(processed_cells)

cat("Normalizing features against reference atlas parameters...\n")
scaled_features <- all_features_df %>%
  dplyr::select(all_of(feature_list)) %>%
  as.matrix()

# Perform standard scaling
scaled_features <- scale(scaled_features, 
                         center = atlas$scaling_parameters$center, 
                         scale = atlas$scaling_parameters$scale)

# Handle NA/NaN/Inf values in scaled features (fill with 0, i.e., the mean of scaled features)
scaled_features[is.na(scaled_features) | is.nan(scaled_features) | is.infinite(scaled_features)] <- 0

cat("Predicting frame-level behaviors using K-Means centroids...\n")
centroids <- atlas$kmeans_model$centers

# Predict function (Euclidean distance to centers)
predict_cluster <- function(row) {
  if (any(is.na(row))) return(3) # Fallback to Cluster 3 (B1)
  distances <- apply(centroids, 1, function(c) sqrt(sum((row - c)^2)))
  val <- which.min(distances)
  if (length(val) == 0) return(3)
  return(val)
}

predicted_idx <- apply(scaled_features, 1, predict_cluster)
predicted_idx <- unlist(predicted_idx)

# Map Cluster Index to Biological labels:
# Cluster 1 -> B2 (Confined Erratic)
# Cluster 2 -> B4 (Flattened/Dynamic)
# Cluster 3 -> B1 (Sessile/Spherical)
# Cluster 4 -> B3 (Directional/Migratory)
# Cluster 5 -> B5 (Adherent/Oscillatory)
mapping <- c("B2", "B4", "B1", "B3", "B5")
biological_states <- mapping[predicted_idx]

all_features_df$RawState <- biological_states

cat("Smoothing behavior sequences with a 3-frame transition filter...\n")

# Transition Smoothing Algorithm:
# Only change state if cell stays in another behavior for >3 frames (at least 4 frames)
smooth_states <- function(states) {
  N <- length(states)
  if (N == 0) return(states)
  smoothed <- states
  curr_state <- states[1]
  t <- 2
  while (t <= N) {
    new_state <- states[t]
    if (new_state != curr_state) {
      if (t + 3 <= N && all(states[t:(t+3)] == new_state)) {
        curr_state <- new_state
        smoothed[t:(t+3)] <- new_state
        t <- t + 4
      } else {
        smoothed[t] <- curr_state
        t <- t + 1
      }
    } else {
      smoothed[t] <- curr_state
      t <- t + 1
    }
  }
  return(smoothed)
}

smoothed_df <- all_features_df %>%
  group_by(SourceFile, CellID) %>%
  arrange(Frame) %>%
  mutate(SmoothedState = smooth_states(RawState)) %>%
  ungroup()

# Calculate transition counts
states_unique <- c("B1", "B2", "B3", "B4", "B5")
transition_counts <- matrix(0, nrow = 5, ncol = 5, dimnames = list(states_unique, states_unique))

cat("Computing Markov transition probability matrix...\n")
cell_tracks <- smoothed_df %>% distinct(SourceFile, CellID)

for (idx in 1:nrow(cell_tracks)) {
  sf <- cell_tracks$SourceFile[idx]
  cid <- cell_tracks$CellID[idx]
  
  track_states <- smoothed_df %>%
    filter(SourceFile == sf, CellID == cid) %>%
    arrange(Frame) %>%
    pull(SmoothedState)
  
  N_states <- length(track_states)
  if (N_states < 2) next
  
  for (t in 1:(N_states - 1)) {
    from_state <- track_states[t]
    to_state <- track_states[t+1]
    transition_counts[from_state, to_state] <- transition_counts[from_state, to_state] + 1
  }
}

# Normalize to probabilities
transition_probs <- transition_counts
for (i in 1:5) {
  row_sum <- sum(transition_probs[i, ])
  if (row_sum > 0) {
    transition_probs[i, ] <- transition_probs[i, ] / row_sum
  }
}

# Save probability matrix
write.csv(transition_probs, output_matrix_path, row.names = TRUE)
cat("✔ SUCCESS: Markov transition probability matrix saved to: ", output_matrix_path, "\n\n")

print(transition_probs)

# Plot directed graph using igraph
cat("\nGenerating Markov state transition diagram...\n")
behaviors_label_full <- c(
  "B1: Sessile/Spherical",
  "B2: Confined Erratic",
  "B3: Directional/Migratory",
  "B4: Flattened/Dynamic",
  "B5: Adherent/Oscillatory"
)

# Colors matching the reference atlas palette exactly
behavior_colors <- c(
  "B1" = "#E8E654", # Yellow
  "B2" = "#5EC468", # Green
  "B3" = "#C780CF", # Purple
  "B4" = "#DB5390", # Pink
  "B5" = "#5DACF5"  # Blue
)

# Create an edge list for transitions (excluding self-transitions for the network plot)
g_edges <- data.frame(from = character(), to = character(), prob = numeric(), stringsAsFactors = FALSE)
for (i in states_unique) {
  for (j in states_unique) {
    if (i != j && transition_probs[i, j] > 0.005) { # Only include transitions > 0.5%
      g_edges <- rbind(g_edges, data.frame(from = i, to = j, prob = transition_probs[i, j]))
    }
  }
}

# Build graph object
g <- graph_from_data_frame(g_edges, directed = TRUE, vertices = data.frame(name = states_unique, label = behaviors_label_full))

# Edge properties
E(g)$width <- E(g)$prob * 15
E(g)$label <- paste0(round(E(g)$prob * 100, 1), "%")
E(g)$color <- adjustcolor("black", alpha.f = 0.5)

# Vertex properties
V(g)$color <- behavior_colors[V(g)$name]
V(g)$size <- 50
V(g)$label.color <- "black"
V(g)$label.font <- 2

# Save transition diagram
png(output_plot_path, width = 800, height = 800, res = 100)
par(mar = c(1, 1, 1, 1))

# Use circular layout
plot(g, 
     layout = layout_in_circle,
     vertex.shape = "circle",
     vertex.label.dist = 0,
     edge.arrow.size = 0.8,
     edge.curved = 0.2,
     edge.label.color = "red",
     edge.label.font = 2,
     main = "Markov Chain State Transition Diagram (Filter: > 3 Frames)")

dev.off()
cat("✔ SUCCESS: Directed transition network plot saved to: ", output_plot_path, "\n")
