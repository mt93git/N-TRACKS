# ---
# ---
# ================================================================================
# SHINY APP: Script G (Scrollable UI)
#
# SCRIPT UPDATE:
#
# --- USER REQUEST [2025-10-08] (4) ---
# 1. SCROLLABLE SIDEBAR: Wrapped the contents of the sidebarPanel in a div
#    with a fixed height (90% of viewport height) and vertical scroll enabled.
#    This keeps the main plot panel visible while scrolling through the
#    extensive sidebar controls.
#
# ... (Previous changes remain) ...
# ================================================================================

# --- 1. Load Required Libraries ---
cat("Loading required libraries...\n")
packages_to_load <- c(
  "shiny", "shinythemes", "tidyverse", "pheatmap", "viridis", "grid", "data.table", "scales",
  "shinyjs", "colourpicker", "RColorBrewer", "MASS", "shinyjqui", "ggrepel", "DT", "readxl", "uwot",
  "plotly", "ggpubr" # ggpubr added for arranging plots
)
for (pkg in packages_to_load) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# --- 2. Constants, Colors, and Lookups ---
HEATMAP_FEATURE_ORDER <- c(
  'TraveledDistance', 'MeanVelocity', 'StdVelocity', 'MaxVelocity', 'MinVelocity',
  'CVvelocity', 'LinearDistance', 'Tortuosity', 'MeanVolume', 'MeanSurfaceArea',
  'MeanSphericity', 'MeanProlateEllipt', 'MeanOblateEllipt', 'MeanMajorAxis',
  'MeanMidAxis', 'MeanMinAxis', 'MeanAR1', 'MeanAR2', 'MeanAR3',
  'DistanceToSurface_Mean', 'AngleOfDirection_Sum', 'StdVolume', 'StdSurfaceArea',
  'StdSphericity', 'StdProlateEllipt', 'StdOblateEllipt', 'StdMajorAxis',
  'StdMidAxis', 'StdMinAxis', 'StdAR1', 'StdAR2', 'StdAR3',
  'DistanceToSurface_Std', 'AngleOfDirection_Std'
)
FEATURE_NAME_LOOKUP <- setNames(
  c(
    "Traveled distance", "Mean velocity", "Std velocity", "Max velocity", "Min velocity",
    "CV velocity", "Linear distance", "Tortuosity (log)", "Volume (mean)",
    "Superficial area (mean)", "Sphericity (mean)", "Prolate ellipticity (mean)",
    "Oblate ellipticity (mean)", "Principal axis length (mean)",
    "Second principal axis length (mean)", "Third principal axis length (mean)",
    "AR1 (mean)", "AR2 (mean)", "AR3 (mean)", "Distance to surface (mean)",
    "Angle of direction (sum)", "Volume (std)", "Superficial area (std)",
    "Sphericity (std)", "Prolate ellipticity (std)", "Oblate ellipticity (std)",
    "Principal axis length (std)", "Second principal axis length (std)",
    "Third principal axis length (std)", "AR1 (std)", "AR2 (std)", "AR3 (std)",
    "Distance to surface (std)", "Angle of direction (std)"
  ),
  HEATMAP_FEATURE_ORDER
)
BEHAVIOR_COLORS <- setNames(
  c('#0173b2', '#de8f05', '#029e73', '#d90000', '#cc78bc', '#949494'),
  c('B1: Sessile/Spherical', 'B2: Confined/Erratic', 'B3: Directional/Migratory', 'B4: Flattened/Dynamic', 'B5: Adherent/Oscillatory', 'Unassigned')
)
DISTINCT_COLORS <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#AA4499", "#DDCC77", "#117733", "#882255", "#44AA99", "#88CCEE", "#332288", "#CC6677",
  "#999933", "#661100", "#6699CC", "#888888"
)

# --- 3. Helper Functions & UI Modules ---
darken_color <- function(color_hex, factor=0.5) {
  tryCatch({
    rgb_vals <- col2rgb(color_hex)
    darker_rgb <- rgb_vals * factor
    rgb(darker_rgb[1,], darker_rgb[2,], darker_rgb[3,], maxColorValue=255)
  }, error = function(e) color_hex)
}

calculate_density_polygons <- function(data, x_col, y_col, group_col, smoothing, limits, normalize = FALSE, nlevels = 20) {
  all_groups <- na.omit(unique(as.character(data[[group_col]])))
  polygon_data <- lapply(all_groups, function(current_group) {
    df_subset <- data %>% filter(.data[[group_col]] == current_group)
    if (nrow(df_subset) < 5) return(NULL)
    bw_x <- smoothing * IQR(df_subset[[x_col]], na.rm=TRUE)
    bw_y <- smoothing * IQR(df_subset[[y_col]], na.rm=TRUE)
    if (is.na(bw_x) || bw_x <= 0.01) bw_x <- 0.1
    if (is.na(bw_y) || bw_y <= 0.01) bw_y <- 0.1
    kde_result <- tryCatch({
      MASS::kde2d(df_subset[[x_col]], df_subset[[y_col]], n=100, h=c(bw_x, bw_y), lims=limits)
    }, error = function(e) NULL)
    if (is.null(kde_result)) return(NULL)
    
    if (isTRUE(normalize)) {
      sorted_z <- sort(kde_result$z, decreasing = TRUE)
      cumulative_z <- cumsum(sorted_z)
      total_sum <- sum(kde_result$z)
      normalized_cumulative_z <- cumulative_z / total_sum
      prob_levels <- seq(0.01, 0.99, length.out = nlevels)
      level_indices <- sapply(prob_levels, function(p) which.min(abs(normalized_cumulative_z - p)))
      contour_levels <- sorted_z[level_indices]
      contour_levels <- unique(sort(contour_levels))
      contourLines(kde_result$x, kde_result$y, kde_result$z, levels = contour_levels)
    } else {
      contourLines(kde_result$x, kde_result$y, kde_result$z, nlevels = nlevels)
    }
  })
  names(polygon_data) <- all_groups
  return(polygon_data)
}

densityControlsUI <- function(id_prefix, s_density_val, s_norm_val, pad_val, smooth_val, s_glow_val, glow_op_val, glow_ext_val, glow_ow_val, s_core_val, core_op_val, core_size_val, core_dark_val, core_ow_val, s_dense_val, dense_op_val, dense_size_val, dense_dark_val, dense_ow_val) {
  ns <- function(id) paste0(id_prefix, "_", id)
  tagList(
    checkboxInput(ns("show_density"), "Enable Density Layers", value = s_density_val),
    conditionalPanel(
      condition = paste0("input['", ns("show_density"), "'] == true"),
      checkboxInput(ns("normalize_density"), "Normalize Density by Group (%)", value = s_norm_val),
      hr(),
      h4("Frame & Shape Controls"),
      sliderInput(ns("padding_factor"), "Plot Padding (%):", min = 5, max = 100, value = pad_val, step = 1),
      sliderInput(ns("smoothing"),"Smoothing:",min=0.1,max=3.0,value = smooth_val,step=0.05), hr(),
      h4("Layer Visibility & Appearance"),
      checkboxInput(ns("show_glow"), "Show Peripheral Glow", value = s_glow_val),
      sliderInput(ns("glow_opacity"), "Glow Opacity:", min=0, max=1, value = glow_op_val, step=0.01),
      sliderInput(ns("glow_extent"), "Glow Extent (# Outer Rings):", min=1, max=20, value = glow_ext_val, step=1),
      sliderInput(ns("glow_outline_width"), "Glow Outline Width:", min=0, max=2, value = glow_ow_val, step=0.1),
      checkboxInput(ns("show_core"), "Show Main Core", value = s_core_val),
      sliderInput(ns("core_opacity"),"Main Core Opacity:",min=0,max=1,value = core_op_val, step=0.01),
      sliderInput(ns("core_size"), "Main Core Size (# Inner Rings):", min=1, max=20, value = core_size_val, step=1),
      sliderInput(ns("darken_factor"), "Main Core Darkening:", min=0, max=1, value = core_dark_val, step=0.01),
      sliderInput(ns("core_outline_width"), "Core Outline Width:", min=0, max=2, value = core_ow_val, step=0.1),
      checkboxInput(ns("show_dense_core"), "Show Dense Core (Top)", value = s_dense_val),
      sliderInput(ns("dense_core_opacity"),"Dense Core Opacity:",min=0,max=1,value = dense_op_val, step=0.01),
      sliderInput(ns("dense_core_size"), "Dense Core Size (# Innermost):", min = 1, max = 10, value = dense_size_val, step = 1),
      sliderInput(ns("dense_core_darkening"), "Dense Core Darkening:", min = 0, max = 1, value = dense_dark_val, step = 0.01),
      sliderInput(ns("dense_core_outline_width"), "Dense Core Outline Width:", min=0, max=2, value = dense_ow_val, step=0.1)
    )
  )
}

pointControlsUI <- function(id_prefix, alpha_val, size_val) {
  ns <- function(id) paste0(id_prefix, "_", id)
  tagList(
    sliderInput(ns("point_alpha"), "Point Opacity:", min = 0, max = 1, value = alpha_val, step = 0.05),
    sliderInput(ns("dot_size"), "Point Size:", min = 0.5, max = 10, value = size_val, step = 0.1)
  )
}

syncedColorInput <- function(id_prefix, label, value) {
  text_input_id <- paste0(id_prefix, "_hex")
  picker_input_id <- paste0(id_prefix, "_picker")
  tagList(
    tags$label(label, `for` = text_input_id),
    div(style = "display: flex; align-items: center;",
        div(style = "flex-grow: 1;",
            textInput(text_input_id, label = NULL, value = value)
        ),
        div(style = "margin-left: 5px;",
            colourInput(picker_input_id, label = NULL, value = value, showColour = "background")
        )
    )
  )
}

# --- 4. UI Definition ---
ui <- fluidPage(theme = shinytheme("yeti"),
                useShinyjs(),
                titlePanel("ACME Definitive Explorer (v49.0 - Unified Curation)"),
                sidebarLayout(
                  sidebarPanel(width = 3,
                               # --- UI MODIFICATION: Add a scrollable div wrapper ---
                               tags$div(
                                 style = "height: 90vh; overflow-y: auto;", # Set fixed height and enable vertical scrolling
                                 h4("Global Controls"),
                                 p("Start by creating a new analysis or loading a saved Atlas."),
                                 fileInput("file_upload", "1a. Upload Base Data (.csv/.xlsx)", accept = c(".csv", ".xlsx")),
                                 p("OR", style="text-align: center;"),
                                 fileInput("atlas_file_upload", "1b. Load Reference Atlas (.rds)", accept = ".rds"),
                                 hr(),
                                 h4("UMAP Projection"),
                                 p("Requires a loaded or generated Atlas.", style="font-size: 0.9em; color: grey;"),
                                 fileInput("projection_file_upload", "2. Upload Data to Project (.csv/.xlsx)", accept = c(".csv", ".xlsx")),
                                 actionButton("project_data_btn", "Project New Data", icon = icon("bullseye"), class = "btn-primary", width="100%"),
                                 hr(),
                                 tabsetPanel(
                                   id = "curation_tabs",
                                   tabPanel("Sample Selection",
                                            wellPanel(
                                              h5("Select Samples to Include:"),
                                              uiOutput("sample_selection_ui")
                                            )
                                   ),
                                   # --- UNIFIED CURATION TAB ---
                                   tabPanel("Data Curation",
                                            wellPanel(
                                              h4("Curation Application Target"),
                                              fluidRow(
                                                column(6, checkboxInput("apply_curation_to_base", "Apply to Base Atlas", value = TRUE)),
                                                column(6, checkboxInput("apply_curation_to_projected", "Apply to Projected Data", value = TRUE))
                                              ),
                                              hr(),
                                              h4("Filtering & Outlier Controls"),
                                              numericInput("filter_min_timepoints", "Filter by Min. Timepoints:", value = 0, min = 0, max = 100),
                                              checkboxInput("filter_na_strict", "Apply Strict NA Filter", value = TRUE),
                                              hr(),
                                              radioButtons("outlier_method", "Outlier Handling (Base Atlas Only):",
                                                           choices = c("None", "Winsorize", "Clip (Filter)"),
                                                           selected = "None", inline = TRUE),
                                              conditionalPanel(
                                                condition = "input.outlier_method == 'Winsorize' || input.outlier_method == 'Clip (Filter)'",
                                                fluidRow(
                                                  column(6, numericInput("winsor_lower", "Lower (%)", value = 5, min = 0, max = 49, step = 1)),
                                                  column(6, numericInput("winsor_upper", "Upper (%)", value = 95, min = 51, max = 100, step = 1))
                                                )
                                              ),
                                              hr(),
                                              checkboxInput("enable_aberrant_filter", "Enable Aberrant Value Filters", value = FALSE),
                                              conditionalPanel(
                                                condition = "input.enable_aberrant_filter == true",
                                                sliderInput("filter_volume", "Mean Volume Range:", min=50, max=6200, value=c(50, 1500)),
                                                sliderInput("filter_surface_area", "Mean Surface Area Range:", min=70, max=3600, value=c(70, 1000)),
                                                sliderInput("filter_distance", "Traveled Distance Range:", min=0, max=400, value=c(0, 250)),
                                                sliderInput("filter_tortuosity", "Tortuosity Range:", min=1, max=220, value=c(1, 40))
                                              ),
                                              hr(),
                                              h4("Clustering & UMAP (Base Atlas Only)"),
                                              numericInput("kmeans_k", "Number of K-Means Clusters (k):", value = 5, min = 2, max = 20, step = 1),
                                              hr(),
                                              h5("UMAP Parameters"),
                                              numericInput("umap_n_neighbors", "Neighbors (n_neighbors):", value = 150, min = 2, max = 200, step = 1),
                                              sliderInput("umap_min_dist", "Min. Distance (min_dist):", value = 0.2, min = 0, max = 1, step = 0.01),
                                              sliderInput("umap_spread", "Spread:", value = 1.0, min = 0.1, max = 5, step = 0.1),
                                              selectInput("umap_metric", "Distance Metric:",
                                                          choices = c("euclidean", "cosine", "manhattan", "hamming"),
                                                          selected = "euclidean")
                                            )
                                   ),
                                   tabPanel("Sub-Clustering",
                                            wellPanel(
                                              uiOutput("sub_cluster_ui")
                                            )
                                   ),
                                   tabPanel("Styling",
                                            wellPanel(
                                              uiOutput("behavior_style_ui"),
                                              hr(),
                                              uiOutput("condition_color_ui"),
                                              hr(),
                                              uiOutput("cumulative_style_ui")
                                            )
                                   )
                                 ),
                                 actionButton("apply_curation", "Apply Curation & Re-Process", icon = icon("cogs"), class = "btn-success", width="100%"),
                                 hr(),
                                 uiOutput("curation_summary_ui"),
                                 hr(),
                                 uiOutput("dynamic_sidebar_ui")
                               ) # End of scrollable div
                  ),
                  mainPanel(width = 9,
                            uiOutput("main_panel_content_ui")
                  )
                )
)

# --- 5. Server Logic ---
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 200*1024^2)
  
  app_state <- reactiveValues(
    raw_data = NULL,
    processed_data = NULL,
    analysis_atlas = NULL,
    heatmap_plot = NULL,
    umap_plot = NULL,
    filter_plot = NULL,
    behavior_pie_plots = list(),
    condition_in_behavior_pie_plots = list(),
    current_k = NULL,
    sample_id_list = NULL,
    behavior_styles = NULL,
    condition_styles = NULL,
    merged_conditions_map = list(),
    cumulative_group_styles = tibble(group = character(), color = character())
  )
  
  sidebar_state <- reactiveValues()
  
  filter_variables <- c("active_group", "display_cluster")
  
  setupColorSyncObserver <- function(input, session, id_prefix) {
    text_input_id <- paste0(id_prefix, "_hex")
    picker_input_id <- paste0(id_prefix, "_picker")
    is_hex_color <- function(s) {
      grepl("^#([A-Fa-f0-9]{8}|[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", s, ignore.case = TRUE)
    }
    observeEvent(input[[text_input_id]], {
      text_val <- input[[text_input_id]]
      if (is_hex_color(text_val) && !is.null(input[[picker_input_id]]) && text_val != input[[picker_input_id]]) {
        updateColourInput(session, picker_input_id, value = text_val)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    observeEvent(input[[picker_input_id]], {
      picker_val <- input[[picker_input_id]]
      if (!is.null(input[[text_input_id]]) && picker_val != input[[text_input_id]]) {
        updateTextInput(session, text_input_id, value = picker_val)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
  }
  
  dynamic_behavior_palette <- reactive({
    req(app_state$behavior_styles)
    setNames(app_state$behavior_styles$color, app_state$behavior_styles$cluster_id)
  })
  
  dynamic_behavior_labels <- reactive({
    req(app_state$behavior_styles)
    setNames(app_state$behavior_styles$label, app_state$behavior_styles$cluster_id)
  })
  
  dynamic_condition_palette <- reactive({
    req(app_state$condition_styles)
    setNames(app_state$condition_styles$color, app_state$condition_styles$condition)
  })
  
  # --- WORKFLOW 1A: Load Base Data from CSV/XLSX ---
  observeEvent(input$file_upload, {
    req(input$file_upload)
    file_path <- input$file_upload$datapath
    df <- if (tools::file_ext(file_path) == "xlsx") {
      readxl::read_excel(file_path)
    } else {
      data.table::fread(file_path)
    }
    if(!("experiment_id" %in% names(df))) {
      showNotification("Error: Uploaded file must contain an 'experiment_id' column.", type="error", duration=15)
      return()
    }
    grouping_col_name <- NULL
    if ('group_clean' %in% names(df)) {
      grouping_col_name <- "group_clean"
      showNotification("Used 'group_clean' for condition labels.", type='message', duration=8)
    } else if ('group' %in% names(df)) {
      grouping_col_name <- "group"
      showNotification("Warning: 'group_clean' not found. Using raw 'group' column.", type='warning', duration=8)
    } else {
      df$primary_group <- 'N/A'
      showNotification("Warning: No 'group_clean' or 'group' column found. All samples assigned to 'N/A' group.", type='warning', duration=10)
    }
    if (!is.null(grouping_col_name)) {
      names(df)[names(df) == grouping_col_name] <- "primary_group"
    }
    app_state$raw_data <- as_tibble(df)
    id_map <- app_state$raw_data %>% distinct(primary_group, experiment_id)
    app_state$sample_id_list <- split(id_map$experiment_id, id_map$primary_group)
    
    # Reset all analysis products
    app_state$processed_data <- NULL
    app_state$analysis_atlas <- NULL
    app_state$current_k <- NULL
    app_state$behavior_styles <- NULL
    app_state$condition_styles <- NULL
    app_state$merged_conditions_map <- list()
    app_state$cumulative_group_styles <- tibble(group = character(), color = character())
    
    showNotification("Base data loaded. Please apply curation settings.", type="message")
  })
  
  # --- WORKFLOW 1B: Load Atlas from RDS ---
  observeEvent(input$atlas_file_upload, {
    req(input$atlas_file_upload)
    withProgress(message = 'Loading Reference Atlas...', value = 0.2, {
      tryCatch({
        loaded_atlas <- readRDS(input$atlas_file_upload$datapath)
        
        # Validation
        incProgress(0.2, detail = "Validating atlas structure...")
        req_names <- c("curated_data", "scaling_parameters", "umap_model", "kmeans_model", "metadata")
        if(!all(req_names %in% names(loaded_atlas))) {
          showNotification("Error: Loaded RDS file is not a valid Reference Atlas.", type="error", duration=15)
          return()
        }
        
        # Reset states
        incProgress(0.2, detail = "Resetting current state...")
        app_state$raw_data <- NULL
        app_state$processed_data <- NULL
        app_state$analysis_atlas <- NULL
        app_state$behavior_styles <- NULL
        app_state$condition_styles <- NULL
        
        # Populate states from Atlas
        incProgress(0.2, detail = "Populating app state...")
        app_state$analysis_atlas <- loaded_atlas
        app_state$processed_data <- loaded_atlas$curated_data
        app_state$raw_data <- loaded_atlas$curated_data # Use curated data as the base for selections
        app_state$current_k <- length(unique(loaded_atlas$curated_data$kmeans_cluster))
        
        # Restore styles and UI elements
        style_meta <- loaded_atlas$metadata$visualization_styles
        app_state$behavior_styles <- tibble(cluster_id = names(style_meta$behavior_labels), label = style_meta$behavior_labels, color = style_meta$behavior_palette[names(style_meta$behavior_labels)])
        app_state$condition_styles <- tibble(condition = names(style_meta$condition_palette), color = style_meta$condition_palette)
        
        id_map <- app_state$processed_data %>% distinct(primary_group, experiment_id)
        app_state$sample_id_list <- split(id_map$experiment_id, id_map$primary_group)
        
        incProgress(0.2, detail = "Done.")
        showNotification("Reference Atlas loaded successfully.", type="message")
        
      }, error = function(e) {
        showNotification(paste("Error loading Atlas:", e$message), type="error", duration=15)
      })
    })
  })
  
  output$sample_selection_ui <- renderUI({
    req(app_state$sample_id_list)
    sample_list <- app_state$sample_id_list
    lapply(names(sample_list), function(group) {
      group_id <- gsub("[^A-Za-z0-9]", "_", group)
      tagList(
        checkboxInput(paste0("group_", group_id), label = tags$strong(group), value = TRUE),
        div(style="margin-left: 20px; max-height: 150px; overflow-y: auto; border: 1px solid #ddd; padding: 5px; border-radius: 4px;",
            checkboxGroupInput(inputId = paste0("ids_", group_id), label = NULL,
                               choices = sample_list[[group]],
                               selected = sample_list[[group]]))
      )
    })
  })
  
  observe({
    req(app_state$sample_id_list)
    lapply(names(app_state$sample_id_list), function(group) {
      group_id <- gsub("[^A-Za-z0-9]", "_", group)
      observeEvent(input[[paste0("group_", group_id)]], {
        if (isTRUE(input[[paste0("group_", group_id)]])) {
          updateCheckboxGroupInput(session, paste0("ids_", group_id), selected = app_state$sample_id_list[[group]])
        } else {
          updateCheckboxGroupInput(session, paste0("ids_", group_id), selected = character(0))
        }
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    })
  })
  
  selected_samples_reactive <- reactive({
    req(app_state$sample_id_list)
    unlist(lapply(names(app_state$sample_id_list), function(group) {
      group_id <- gsub("[^A-Za-z0-9]", "_", group)
      input[[paste0("ids_", group_id)]]
    }))
  })
  
  reprocess_data_clusters_umap <- function(data_to_process) {
    metrics_to_use <- intersect(HEATMAP_FEATURE_ORDER, names(data_to_process))
    metrics_to_use <- setdiff(metrics_to_use, c("DistanceToSurface_Mean", "DistanceToSurface_Std"))
    numeric_metrics_to_use <- metrics_to_use[sapply(data_to_process[metrics_to_use], is.numeric)]
    
    if(nrow(data_to_process) < input$kmeans_k || length(numeric_metrics_to_use) == 0) {
      showNotification("Not enough data or numeric features after filtering for the selected k.", type="error")
      return(NULL)
    }
    
    features_for_scaling <- data_to_process %>% dplyr::select(all_of(numeric_metrics_to_use))
    scaled_features <- scale(features_for_scaling)
    scaling_parameters <- list(
      center = attr(scaled_features, "scaled:center"),
      scale = attr(scaled_features, "scaled:scale")
    )
    
    set.seed(42)
    kmeans_model <- kmeans(scaled_features, centers = input$kmeans_k, nstart = 25)
    data_to_process$kmeans_cluster <- paste0("B", kmeans_model$cluster)
    
    n_neighbors_val <- input$umap_n_neighbors
    if (n_neighbors_val >= nrow(scaled_features)) {
      n_neighbors_val <- nrow(scaled_features) - 1
      showNotification(paste("UMAP n_neighbors was capped at", n_neighbors_val, " (rows-1) to prevent errors."), type="warning")
    }
    
    set.seed(42)
    umap_model <- uwot::umap(
      scaled_features,
      n_neighbors = n_neighbors_val,
      min_dist = input$umap_min_dist,
      spread = input$umap_spread,
      metric = input$umap_metric,
      ret_model = TRUE
    )
    
    data_to_process$feature_umap_x <- umap_model$embedding[, 1]
    data_to_process$feature_umap_y <- umap_model$embedding[, 2]
    data_to_process$cell_id <- 1:nrow(data_to_process)
    data_to_process$data_source <- "Base" # Assign data source
    
    if('Biological_Behavior' %in% names(data_to_process)){
      data_to_process <- data_to_process %>% dplyr::select(-Biological_Behavior)
    }
    
    return(list(
      processed_data = data_to_process,
      scaling_parameters = scaling_parameters,
      kmeans_model = kmeans_model,
      umap_model = umap_model,
      feature_list = numeric_metrics_to_use
    ))
  }
  
  observeEvent(input$apply_curation, {
    req(app_state$raw_data, input$kmeans_k, selected_samples_reactive())
    withProgress(message = 'Curating and Re-processing Data...', value = 0, {
      incProgress(0.1, detail = "Filtering samples...")
      selected_ids <- selected_samples_reactive()
      validate(need(length(selected_ids) > 0, "No samples selected. Please select at least one sample to analyze."))
      data <- app_state$raw_data %>% filter(experiment_id %in% selected_ids)
      validate(need(nrow(data) > 0, "No data remaining after sample selection."))
      
      # --- CONDITIONAL FILTERING FOR BASE DATA ---
      # This block is now gated by the 'apply_curation_to_base' checkbox.
      if(isTRUE(input$apply_curation_to_base)) {
        incProgress(0.1, detail = "Applying base data filters...")
        if (input$filter_min_timepoints > 0) { data <- data %>% filter(NumTimepoints >= input$filter_min_timepoints) }
        if (isTRUE(input$filter_na_strict)) {
          na_check_features <- intersect(HEATMAP_FEATURE_ORDER, names(data));
          na_check_features <- setdiff(na_check_features, c("DistanceToSurface_Mean", "DistanceToSurface_Std"))
          data <- data %>% drop_na(all_of(na_check_features))
        }
        if (isTRUE(input$enable_aberrant_filter)) {
          data <- data %>%
            filter(between(MeanVolume, input$filter_volume[1], input$filter_volume[2])) %>%
            filter(between(MeanSurfaceArea, input$filter_surface_area[1], input$filter_surface_area[2])) %>%
            filter(between(TraveledDistance, input$filter_distance[1], input$filter_distance[2])) %>%
            filter(between(Tortuosity, input$filter_tortuosity[1], input$filter_tortuosity[2]))
        }
        if (input$outlier_method != "None") {
          features_to_process <- intersect(HEATMAP_FEATURE_ORDER, names(data))
          features_to_process <- features_to_process[sapply(data[features_to_process], is.numeric)]
          lower_b <- input$winsor_lower / 100
          upper_b <- input$winsor_upper / 100
          if (input$outlier_method == "Winsorize") {
            incProgress(0.05, detail = "Winsorizing outliers...")
            winsorize <- function(x, l=lower_b, u=upper_b) {
              q <- quantile(x, c(l,u), na.rm=TRUE)
              x[x < q[1]] <- q[1]
              x[x > q[2]] <- q[2]
              return(x)
            }
            data <- data %>% mutate(across(all_of(features_to_process), ~ winsorize(.)))
          } else if (input$outlier_method == "Clip (Filter)") {
            incProgress(0.05, detail = "Clipping outliers...")
            for (feat in features_to_process) {
              if (is.numeric(data[[feat]])) {
                q_vals <- quantile(data[[feat]], c(lower_b, upper_b), na.rm = TRUE)
                if(all(!is.na(q_vals))) {
                  data <- data %>% filter(.data[[feat]] >= q_vals[1] & .data[[feat]] <= q_vals[2])
                }
              }
            }
          }
        }
      } else {
        incProgress(0.2, detail = "Skipping base data filters...")
        showNotification("Curation filters were not applied to the base data as per checkbox setting.", type="message")
      }
      
      incProgress(0.6, detail = "Running Clustering & UMAP...")
      analysis_result <- reprocess_data_clusters_umap(data)
      
      if (!is.null(analysis_result)) {
        app_state$analysis_atlas <- analysis_result
        processed_data <- analysis_result$processed_data
        
        processed_data$display_cluster <- processed_data$kmeans_cluster
        processed_data$active_group <- processed_data$primary_group
        
        app_state$processed_data <- processed_data
        app_state$current_k <- input$kmeans_k
        
        app_state$merged_conditions_map <- list()
        app_state$cumulative_group_styles <- tibble(group = character(), color = character())
        
        k_val <- app_state$current_k
        default_styles_df <- tibble(
          cluster_id = c("B1", "B2", "B3", "B4", "B5"),
          label = c("B1: Sessile/Spherical", "B2: Confined Erratic", "B3: Directional/Migratory", "B4: Flattened/Dynamic", "B5: Adherent/Oscillatory"),
          color = c("#E8E654", "#5EC468", "#C780CF", "#DB5390", "#5DACF5")
        )
        if (k_val == 5) {
          new_b_styles <- default_styles_df
        } else {
          new_b_styles <- tibble(
            cluster_id = paste0("B", 1:k_val),
            label = paste0("Behavior ", 1:k_val),
            color = DISTINCT_COLORS[1:k_val]
          )
        }
        app_state$behavior_styles <- new_b_styles
        
        unique_conditions <- sort(unique(processed_data$primary_group))
        new_c_styles <- tibble(
          condition = unique_conditions,
          color = rep(RColorBrewer::brewer.pal(max(3, min(8, length(unique_conditions))), "Accent"), length.out = length(unique_conditions))
        )
        if (!is.null(app_state$condition_styles)) {
          old_styles <- app_state$condition_styles
          new_c_styles <- new_c_styles %>%
            left_join(old_styles, by = "condition", suffix = c("", ".old")) %>%
            mutate(color = ifelse(!is.na(color.old), color.old, color)) %>%
            dplyr::select(condition, color)
        }
        app_state$condition_styles <- new_c_styles
        
        app_state$heatmap_plot <- NULL
        app_state$umap_plot <- NULL
        app_state$filter_plot <- NULL
        app_state$behavior_pie_plots <- list()
      }
      incProgress(0.1, detail = "Done.")
    })
    if (!is.null(app_state$processed_data)) {
      showNotification("Curation and re-processing complete!", type="message")
    }
  })
  
  # --- WORKFLOW 2: Project New Data onto Atlas ---
  observeEvent(input$project_data_btn, {
    req(input$projection_file_upload, app_state$analysis_atlas)
    
    withProgress(message = 'Projecting new data...', value = 0, {
      incProgress(0.1, detail = "Resetting previous projection...")
      if ("data_source" %in% names(app_state$processed_data)) {
        app_state$processed_data <- app_state$processed_data %>%
          filter(data_source == "Base")
      }
      
      incProgress(0.1, detail = "Loading projection data...")
      proj_file_path <- input$projection_file_upload$datapath
      
      proj_data <- if (tools::file_ext(proj_file_path) == "xlsx") {
        readxl::read_excel(proj_file_path)
      } else {
        data.table::fread(proj_file_path)
      }
      
      if (!("experiment_id" %in% names(proj_data))) {
        showNotification("Error: Projection file must contain an 'experiment_id' column.", type = "error", duration = 15)
        return()
      }
      
      grouping_col_name <- NULL
      if ('group_clean' %in% names(proj_data)) {
        grouping_col_name <- "group_clean"
      } else if ('group' %in% names(proj_data)) {
        grouping_col_name <- "group"
      }
      
      if (!is.null(grouping_col_name)) {
        names(proj_data)[names(proj_data) == grouping_col_name] <- "primary_group"
      } else {
        proj_data$primary_group <- proj_data$experiment_id
      }
      
      proj_data$experiment_id <- as.character(proj_data$experiment_id)
      proj_data$primary_group <- as.character(proj_data$primary_group)
      
      # --- CONDITIONAL FILTERING FOR PROJECTED DATA ---
      # This block is now gated by the 'apply_curation_to_projected' checkbox.
      if(isTRUE(input$apply_curation_to_projected)){
        incProgress(0.1, detail = "Applying projection filters...")
        tracks_before_filter <- nrow(proj_data)
        
        if (input$filter_min_timepoints > 0) {
          proj_data <- proj_data %>% filter(NumTimepoints >= input$filter_min_timepoints)
        }
        if (isTRUE(input$filter_na_strict)) {
          na_check_features <- intersect(HEATMAP_FEATURE_ORDER, names(proj_data))
          na_check_features <- setdiff(na_check_features, c("DistanceToSurface_Mean", "DistanceToSurface_Std"))
          proj_data <- proj_data %>% drop_na(all_of(na_check_features))
        }
        if (isTRUE(input$enable_aberrant_filter)) {
          proj_data <- proj_data %>%
            filter(between(MeanVolume, input$filter_volume[1], input$filter_volume[2])) %>%
            filter(between(MeanSurfaceArea, input$filter_surface_area[1], input$filter_surface_area[2])) %>%
            filter(between(TraveledDistance, input$filter_distance[1], input$filter_distance[2])) %>%
            filter(between(Tortuosity, input$filter_tortuosity[1], input$filter_tortuosity[2]))
        }
        
        tracks_after_filter <- nrow(proj_data)
        showNotification(paste0("Projection Filter: Removed ", tracks_before_filter - tracks_after_filter, " tracks. ", tracks_after_filter, " remaining."), type="message")
        
        if (tracks_after_filter == 0) {
          showNotification("No data remains in the projection dataset after filtering.", type="error", duration=15)
          return()
        }
      } else {
        incProgress(0.1, detail = "Skipping projection filters...")
        showNotification("Curation filters were not applied to the projected data as per checkbox setting.", type="message")
      }
      
      metrics_to_use <- app_state$analysis_atlas$feature_list
      if (!all(metrics_to_use %in% names(proj_data))) {
        missing_cols <- setdiff(metrics_to_use, names(proj_data))
        showNotification(paste("Error: Projection data is missing required columns:", paste(missing_cols, collapse = ", ")), type = "error", duration = 15)
        return()
      }
      
      incProgress(0.2, detail = "Scaling new data...")
      proj_data_numeric <- proj_data %>% dplyr::select(all_of(metrics_to_use))
      proj_scaled <- scale(proj_data_numeric, center = app_state$analysis_atlas$scaling_parameters$center, scale = app_state$analysis_atlas$scaling_parameters$scale)
      
      incProgress(0.2, detail = "Projecting onto UMAP...")
      projected_coords <- uwot::umap_transform(X = proj_scaled, model = app_state$analysis_atlas$umap_model)
      
      incProgress(0.1, detail = "Predicting cluster assignments...")
      centroids <- app_state$analysis_atlas$kmeans_model$centers
      
      predict_cluster <- function(point, centers) {
        distances <- apply(centers, 1, function(c) sqrt(sum((point - c)^2)))
        return(which.min(distances))
      }
      predicted_clusters_idx <- apply(proj_scaled, 1, predict_cluster, centers = centroids)
      
      max_base_cell_id <- max(app_state$processed_data$cell_id, na.rm = TRUE)
      
      projected_data_full <- as_tibble(proj_data) %>%
        mutate(
          feature_umap_x = projected_coords[, 1],
          feature_umap_y = projected_coords[, 2],
          kmeans_cluster = paste0("B", predicted_clusters_idx),
          display_cluster = kmeans_cluster,
          active_group = primary_group,
          data_source = "Projected",
          cell_id = max_base_cell_id + row_number()
        )
      
      incProgress(0.1, detail = "Integrating data...")
      app_state$processed_data <- bind_rows(app_state$processed_data, projected_data_full)
      
      id_map_updated <- app_state$processed_data %>% distinct(primary_group, experiment_id)
      app_state$sample_id_list <- split(id_map_updated$experiment_id, id_map_updated$primary_group)
      
      new_conditions <- setdiff(unique(projected_data_full$primary_group), app_state$condition_styles$condition)
      if(length(new_conditions) > 0) {
        new_colors_needed <- length(new_conditions)
        available_colors <- setdiff(DISTINCT_COLORS, app_state$condition_styles$color)
        if(length(available_colors) < new_colors_needed) {
          available_colors <- rep(DISTINCT_COLORS, length.out = new_colors_needed)
        }
        
        new_styles_df <- tibble(
          condition = new_conditions,
          color = available_colors[1:new_colors_needed]
        )
        app_state$condition_styles <- bind_rows(app_state$condition_styles, new_styles_df)
      }
      
      incProgress(0.1, detail = "Done.")
    })
    
    if(exists("projected_coords")) {
      showNotification(paste("Successfully projected and integrated", nrow(projected_coords), "new data points from", length(unique(proj_data$experiment_id)), "experiment(s)."), type = "message")
    }
  })
  
  output$curation_summary_ui <- renderUI({
    req(app_state$processed_data)
    
    base_data <- app_state$processed_data %>% filter(data_source == "Base")
    projected_data <- app_state$processed_data %>% filter(data_source == "Projected")
    
    tagList(
      p(strong("Base Tracks in Analysis:"), br(), nrow(base_data)),
      if (nrow(projected_data) > 0) {
        p(strong("Projected Tracks in Analysis:"), br(), nrow(projected_data))
      },
      p(strong("Total Tracks Displayed:"), br(), nrow(app_state$processed_data))
    )
  })
  
  output$main_panel_content_ui <- renderUI({
    if (is.null(app_state$processed_data) && is.null(app_state$raw_data)) return(wellPanel(style="text-align:center; padding: 50px;", h3("Please upload base data or load an Atlas to begin.")))
    if (is.null(app_state$processed_data) && !is.null(app_state$raw_data)) return(wellPanel(style="text-align:center; padding: 50px;", h3("Please select samples and apply curation settings.")))
    
    tabsetPanel(id = "main_tabs",
                tabPanel("UMAP Comparator", value = "umap_comparator", jqui_resizable(plotOutput("umap_plot", height="600px")), hr(), DT::dataTableOutput("umap_data_table")),
                tabPanel("Heatmap Explorer", value = "heatmap_explorer",
                         tabsetPanel(
                           id = "heatmap_subtabs",
                           tabPanel("Global View",
                                    jqui_resizable(plotOutput("heatmap_plot", height="800px"))
                           ),
                           tabPanel("Filtered View",
                                    jqui_resizable(plotOutput("filtered_heatmap_plot", height="800px"))
                           )
                         )
                ),
                tabPanel("Interactive Data Filter", value = "interactive_filter",
                         tabsetPanel(
                           type = "tabs",
                           tabPanel("Merged View",
                                    jqui_resizable(plotOutput("filter_plot", height="550px"))
                           ),
                           tabPanel("Individual Sample View",
                                    div(style = "text-align: center; margin-top: 10px; margin-bottom: 10px;",
                                        downloadButton("download_individual_plots_pdf", "Download All as PDF", class = "btn-info")
                                    ),
                                    jqui_resizable(uiOutput("individual_sample_plots_ui"))
                           )
                         ),
                         hr(),
                         tabsetPanel(
                           type = "tabs",
                           tabPanel("Behavior Distribution by Condition",
                                    h4("Behavior Distribution by Condition", style = "text-align: center; font-weight: bold;"),
                                    p("Shows the composition of behaviors within each selected condition.", style="text-align:center; color:grey;"),
                                    downloadButton("download_behavior_dist_pdf", "Download Plots as PDF", class="btn-sm", style="float: right; margin: 10px;"),
                                    jqui_resizable(uiOutput("behavior_pie_charts_ui"))
                           ),
                           tabPanel("Condition Distribution by Behavior",
                                    h4("Condition Distribution by Behavior", style = "text-align: center; font-weight: bold;"),
                                    p("Shows the composition of conditions within each behavior cluster.", style="text-align:center; color:grey;"),
                                    downloadButton("download_condition_dist_pdf", "Download Plots as PDF", class="btn-sm", style="float: right; margin: 10px;"),
                                    uiOutput("behavior_condition_dist_ui")
                           )
                         )
                )
    )
  })
  
  sidebar_input_ids <- c(
    "dist_summary_level", "filter_active_group", "filter_display_cluster",
    "heatmap_cluster_order", "umap_color", "umap_quantile_lower", "umap_quantile_upper",
    "umap_layer_order", "umap_points_point_alpha", "umap_points_dot_size",
    "umap_density_show_density", "umap_density_normalize_density", "umap_density_padding_factor", "umap_density_smoothing",
    "umap_density_show_glow", "umap_density_glow_opacity", "umap_density_glow_extent",
    "umap_density_show_core", "umap_density_core_opacity", "umap_density_core_size",
    "umap_density_darken_factor", "umap_density_show_dense_core", "umap_density_dense_core_opacity",
    "umap_density_dense_core_size", "umap_density_dense_core_darkening",
    "umap_density_glow_outline_width", "umap_density_core_outline_width", "umap_density_dense_core_outline_width",
    "cumulative_filter", "dist_plot_mode", "filter_layer_order", "filter_points_point_alpha",
    "filter_points_dot_size", "filter_density_show_density", "filter_density_normalize_density",
    "filter_density_padding_factor", "filter_density_smoothing", "filter_density_show_glow",
    "filter_density_glow_opacity", "filter_density_glow_extent", "filter_density_show_core",
    "filter_density_core_opacity", "filter_density_core_size", "filter_density_darken_factor",
    "filter_density_show_dense_core", "filter_density_dense_core_opacity",
    "filter_density_dense_core_size", "filter_density_dense_core_darkening",
    "filter_density_glow_outline_width", "filter_density_core_outline_width", "filter_density_dense_core_outline_width",
    "filter_grey_points_point_alpha", "filter_grey_points_dot_size",
    "filter_color_continuous", "filter_quantile_lower", "filter_quantile_upper",
    "filter_density_greyscale", "donut_hole_size"
  )
  lapply(sidebar_input_ids, function(id) {
    observeEvent(input[[id]], {
      sidebar_state[[id]] <- input[[id]]
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
  })
  
  output$dynamic_sidebar_ui <- renderUI({
    req(app_state$processed_data, input$main_tabs, app_state$current_k)
    data <- app_state$processed_data
    `%||%` <- rlang::`%||%`
    if (input$main_tabs == "heatmap_explorer") {
      tagList(
        h4("Heatmap Controls"),
        orderInput(
          inputId = "heatmap_cluster_order",
          label = "Drag to Reorder Heatmap Clusters:",
          items = isolate(sidebar_state$heatmap_cluster_order) %||% app_state$behavior_styles$cluster_id
        ),
        hr(),
        downloadButton("download_heatmap_pdf", "Download Global Heatmap (PDF)", class = "btn-info", width = "100%"),
        br(),br(),
        downloadButton("download_filtered_heatmap_pdf", "Download Filtered Heatmap (PDF)", class = "btn-info", width = "100%")
      )
    } else if (input$main_tabs == "umap_comparator") {
      available_cols <- intersect(HEATMAP_FEATURE_ORDER, names(data))
      umap_feature_choices <- setNames(names(FEATURE_NAME_LOOKUP[available_cols]), FEATURE_NAME_LOOKUP[available_cols])
      categorical_choices <- c('Condition' = 'active_group', 'Behavior'='display_cluster')
      if('Biological_Behavior' %in% names(data)){
        categorical_choices <- c(categorical_choices, 'Biological Behavior' = 'Biological_Behavior')
      }
      tagList(
        h4("UMAP Plot Controls"),
        selectInput('umap_color','Color by:',
                    choices = list(`Categorical` = categorical_choices, `Continuous` = umap_feature_choices),
                    selected = isolate(sidebar_state$umap_color) %||% 'display_cluster'),
        conditionalPanel(
          condition = paste0("input.umap_color && ['", paste(HEATMAP_FEATURE_ORDER, collapse="','"), "'].includes(input.umap_color)"),
          hr(),
          h4("Color Scale Controls"),
          p("Set the quantile range for the gradient.", style="font-size: 0.9em; color: grey;"),
          fluidRow(
            column(6, numericInput("umap_quantile_lower", "Lower (%)", value = isolate(sidebar_state$umap_quantile_lower) %||% 0, min = 0, max = 100, step = 1)),
            column(6, numericInput("umap_quantile_upper", "Upper (%)", value = isolate(sidebar_state$umap_quantile_upper) %||% 100, min = 0, max = 100, step = 1))
          )
        ),
        hr(), h4("Point Appearance"),
        pointControlsUI("umap_points",
                        alpha_val = isolate(sidebar_state$umap_points_point_alpha) %||% 0.7,
                        size_val = isolate(sidebar_state$umap_points_dot_size) %||% 2.5),
        hr(),
        h4("Layer Order"),
        orderInput("umap_layer_order", "Drag to Reorder Layers",
                   items = isolate(sidebar_state$umap_layer_order) %||% c("Points", "Glow", "Core", "Dense Core")),
        hr(),
        densityControlsUI("umap_density",
                          s_density_val = isolate(sidebar_state$umap_density_show_density) %||% FALSE,
                          s_norm_val = isolate(sidebar_state$umap_density_normalize_density) %||% FALSE,
                          pad_val= isolate(sidebar_state$umap_density_padding_factor) %||% 10,
                          smooth_val = isolate(sidebar_state$umap_density_smoothing) %||% 0.7,
                          s_glow_val = isolate(sidebar_state$umap_density_show_glow) %||% TRUE,
                          glow_op_val = isolate(sidebar_state$umap_density_glow_opacity) %||% 0.14,
                          glow_ext_val = isolate(sidebar_state$umap_density_glow_extent) %||% 20,
                          glow_ow_val = isolate(sidebar_state$umap_density_glow_outline_width) %||% 0,
                          s_core_val = isolate(sidebar_state$umap_density_show_core) %||% TRUE,
                          core_op_val = isolate(sidebar_state$umap_density_core_opacity) %||% 0.1,
                          core_size_val = isolate(sidebar_state$umap_density_core_size) %||% 18,
                          core_dark_val = isolate(sidebar_state$umap_density_darken_factor) %||% 0.92,
                          core_ow_val = isolate(sidebar_state$umap_density_core_outline_width) %||% 0,
                          s_dense_val = isolate(sidebar_state$umap_density_show_dense_core) %||% TRUE,
                          dense_op_val = isolate(sidebar_state$umap_density_dense_core_opacity) %||% 0.24,
                          dense_size_val = isolate(sidebar_state$umap_density_dense_core_size) %||% 5,
                          dense_dark_val = isolate(sidebar_state$umap_density_dense_core_darkening) %||% 0.65,
                          dense_ow_val = isolate(sidebar_state$umap_density_dense_core_outline_width) %||% 0),
        hr(),
        h4("Downloads"),
        downloadButton("download_umap_pdf", "Download UMAP (PDF)", class = "btn-info", width = "100%"),
        br(),br(),
        downloadButton("download_culled_data", "Download Plotted Data (.csv)", class="btn-primary", width="100%"),
        br(),br(),
        downloadButton("download_atlas", "Export as Reference Atlas (.rds)", class = "btn-success", icon = icon("save"), width = "100%")
      )
    } else if (input$main_tabs == "interactive_filter") {
      filter_checkbox_ui <- lapply(filter_variables, function(var_name) {
        if(var_name %in% names(data) && !all(is.na(data[[var_name]]))) {
          choices <- sort(unique(data[[var_name]]))
          if (var_name == "display_cluster") {
            dynamic_labels <- dynamic_behavior_labels()
            names(choices) <- sapply(choices, function(c) { lbl <- dynamic_labels[c]; if (is.na(lbl) || is.null(lbl)) c else lbl })
          }
          label_text <- if (var_name == "display_cluster") "Filter by Behavior" else "Filter by Condition"
          checkboxGroupInput(
            inputId = paste0("filter_", var_name),
            label = label_text,
            choices = choices,
            selected = isolate(sidebar_state[[paste0("filter_", var_name)]])
          )
        }
      })
      
      available_cols_filter <- intersect(HEATMAP_FEATURE_ORDER, names(data))
      filter_feature_choices <- setNames(
        names(FEATURE_NAME_LOOKUP[available_cols_filter]),
        FEATURE_NAME_LOOKUP[available_cols_filter]
      )
      
      tagList(
        h4("Create Merged Condition"),
        wellPanel(
          uiOutput("condition_merging_ui"),
          uiOutput("merged_conditions_list_ui")
        ),
        hr(),
        h4("Interactive Filtering Controls"),
        p("Select conditions/behaviors to highlight on the UMAP and generate summary plots.", style="font-size: 0.9em; color: grey;"),
        
        selectInput('filter_color_continuous', 'Color by (Highlighted Points):',
                    choices = c('Categorical Highlight' = 'categorical', `Continuous Feature` = filter_feature_choices),
                    selected = isolate(sidebar_state$filter_color_continuous) %||% 'categorical'),
        conditionalPanel(
          condition = "input.filter_color_continuous != 'categorical'",
          wellPanel(
            style = "background-color: #f8f9fa;",
            h5("Continuous Color Scale"),
            p("Set the quantile range for the gradient.", style="font-size: 0.9em; color: grey;"),
            fluidRow(
              column(6, numericInput("filter_quantile_lower", "Lower (%)", value = isolate(sidebar_state$filter_quantile_lower) %||% 0, min = 0, max = 100, step = 1)),
              column(6, numericInput("filter_quantile_upper", "Upper (%)", value = isolate(sidebar_state$filter_quantile_upper) %||% 100, min = 0, max = 100, step = 1))
            )
          )
        ),
        
        checkboxInput("cumulative_filter", "Enable Cumulative (AND) Filtering", value = isolate(sidebar_state$cumulative_filter) %||% FALSE),
        hr(),
        h4("Distribution Plot Mode"),
        radioButtons("dist_plot_mode", "Display Mode:",
                     choices = c("Relative Percentage (Donut Chart)", "Frequency (Bar Chart)"),
                     selected = isolate(sidebar_state$dist_plot_mode) %||% "Relative Percentage (Donut Chart)"),
        
        conditionalPanel(
          condition = "input.dist_plot_mode == 'Relative Percentage (Donut Chart)'",
          sliderInput("donut_hole_size", "Donut Hole Size:",
                      min = 0, max = 0.8, value = isolate(sidebar_state$donut_hole_size) %||% 0, step = 0.05)
        ),
        
        hr(),
        h4("Distribution Summary Level"),
        radioButtons("dist_summary_level", "Display by:",
                     choices = c("Group by Condition", "Show Individual Samples"),
                     selected = isolate(sidebar_state$dist_summary_level) %||% "Group by Condition"),
        hr(),
        h4("Point Appearance (Highlighted)"),
        pointControlsUI("filter_points",
                        alpha_val = isolate(sidebar_state$filter_points_point_alpha) %||% 0.05,
                        size_val = isolate(sidebar_state$filter_points_dot_size) %||% 8.1),
        hr(),
        h4("Default Grey Dots Appearance"),
        pointControlsUI("filter_grey_points",
                        alpha_val = isolate(sidebar_state$filter_grey_points_point_alpha) %||% 0.35,
                        size_val = isolate(sidebar_state$filter_grey_points_dot_size) %||% 6.5),
        hr(),
        h4("Layer Order"),
        orderInput("filter_layer_order", "Drag to Reorder Layers",
                   items = isolate(sidebar_state$filter_layer_order) %||% c("Points", "Glow", "Core", "Dense Core")),
        hr(),
        densityControlsUI("filter_density",
                          s_density_val = isolate(sidebar_state$filter_density_show_density) %||% FALSE,
                          s_norm_val = isolate(sidebar_state$filter_density_normalize_density) %||% FALSE,
                          pad_val = isolate(sidebar_state$filter_density_padding_factor) %||% 12,
                          smooth_val = isolate(sidebar_state$filter_density_smoothing) %||% 0.5,
                          s_glow_val = isolate(sidebar_state$filter_density_show_glow) %||% TRUE,
                          glow_op_val = isolate(sidebar_state$filter_density_glow_opacity) %||% 0.23,
                          glow_ext_val = isolate(sidebar_state$filter_density_glow_extent) %||% 3,
                          glow_ow_val = isolate(sidebar_state$filter_density_glow_outline_width) %||% 0.2,
                          s_core_val = isolate(sidebar_state$filter_density_show_core) %||% TRUE,
                          core_op_val = isolate(sidebar_state$filter_density_core_opacity) %||% 0.07,
                          core_size_val = isolate(sidebar_state$filter_density_core_size) %||% 18,
                          core_dark_val = isolate(sidebar_state$filter_density_darken_factor) %||% 0.67,
                          core_ow_val = isolate(sidebar_state$filter_density_core_outline_width) %||% 0.1,
                          s_dense_val = isolate(sidebar_state$filter_density_show_dense_core) %||% TRUE,
                          dense_op_val = isolate(sidebar_state$filter_density_dense_core_opacity) %||% 0.3,
                          dense_size_val = isolate(sidebar_state$filter_density_dense_core_size) %||% 5,
                          dense_dark_val = isolate(sidebar_state$filter_density_dense_core_darkening) %||% 0.56,
                          dense_ow_val = isolate(sidebar_state$filter_density_dense_core_outline_width) %||% 0.2),
        conditionalPanel(
          condition = "input.filter_density_show_density == true",
          checkboxInput("filter_density_greyscale", "Use Greyscale for Density",
                        value = isolate(sidebar_state$filter_density_greyscale) %||% TRUE)
        ),
        hr(),
        filter_checkbox_ui,
        hr(),
        h4("Filter by Sample"),
        uiOutput("filter_sample_selection_ui"),
        hr(),
        downloadButton("download_filter_plot_pdf", "Download Filter Plot (PDF)", class = "btn-info", width = "100%")
      )
    }
  })
  
  output$filter_sample_selection_ui <- renderUI({
    req(app_state$processed_data)
    id_map <- app_state$processed_data %>% distinct(primary_group, experiment_id)
    sample_list <- split(id_map$experiment_id, id_map$primary_group)
    
    lapply(names(sample_list), function(group) {
      group_id <- gsub("[^A-Za-z0-9]", "_", group)
      tagList(
        checkboxInput(paste0("filter_group_", group_id), label = tags$strong(group), value = TRUE),
        div(style="margin-left: 20px; max-height: 150px; overflow-y: auto; border: 1px solid #ddd; padding: 5px; border-radius: 4px;",
            checkboxGroupInput(inputId = paste0("filter_ids_", group_id), label = NULL,
                               choices = sample_list[[group]],
                               selected = sample_list[[group]]))
      )
    })
  })
  
  observe({
    req(app_state$processed_data)
    id_map <- app_state$processed_data %>% distinct(primary_group, experiment_id)
    sample_list <- split(id_map$experiment_id, id_map$primary_group)
    
    lapply(names(sample_list), function(group) {
      group_id <- gsub("[^A-Za-z0-9]", "_", group)
      observeEvent(input[[paste0("filter_group_", group_id)]], {
        if (isTRUE(input[[paste0("filter_group_", group_id)]])) {
          updateCheckboxGroupInput(session, paste0("filter_ids_", group_id), selected = sample_list[[group]])
        } else {
          updateCheckboxGroupInput(session, paste0("filter_ids_", group_id), selected = character(0))
        }
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    })
  })
  
  filtered_samples_reactive <- reactive({
    req(app_state$processed_data)
    id_map <- app_state$processed_data %>% distinct(primary_group, experiment_id)
    sample_list <- split(id_map$experiment_id, id_map$primary_group)
    unlist(lapply(names(sample_list), function(group) {
      group_id <- gsub("[^A-Za-z0-9]", "_", group)
      input[[paste0("filter_ids_", group_id)]]
    }))
  })
  
  output$sub_cluster_ui <- renderUI({
    if (is.null(app_state$processed_data)) {
      return(p("Process data first to enable sub-clustering."))
    }
    
    parent_choices <- app_state$behavior_styles$cluster_id
    
    tagList(
      h4("Hierarchical Sub-Clustering"),
      p("Select a parent behavior and split it into a specified number of sub-behaviors.", style="font-size: 0.9em; color: grey;"),
      hr(),
      selectInput("parent_cluster_select", "1. Select Parent Behavior:", choices = parent_choices),
      numericInput("sub_cluster_k", "2. Number of Sub-Behaviors (k):", value = 2, min = 2, max = 10, step = 1),
      actionButton("apply_sub_clustering", "3. Generate Sub-Behaviors", icon=icon("sitemap"), class="btn-primary", width="100%"),
      hr(),
      actionButton("clear_sub_clustering", "Clear All Sub-Behaviors", icon=icon("undo"), class="btn-warning", width="100%")
    )
  })
  
  observeEvent(input$apply_sub_clustering, {
    req(input$parent_cluster_select, input$sub_cluster_k, app_state$processed_data)
    
    parent_cluster <- input$parent_cluster_select
    k_sub <- input$sub_cluster_k
    
    data_subset_indices <- which(app_state$processed_data$display_cluster == parent_cluster)
    data_subset <- app_state$processed_data[data_subset_indices, ]
    
    if (nrow(data_subset) < k_sub) {
      showNotification(paste("Not enough data points in", parent_cluster, "to create", k_sub, "clusters."), type="error")
      return()
    }
    
    withProgress(message = paste("Sub-clustering", parent_cluster, "..."), value = 0.3, {
      metrics_to_use <- app_state$analysis_atlas$feature_list
      numeric_metrics_to_use <- metrics_to_use[sapply(data_subset[metrics_to_use], is.numeric)]
      
      scaled_features_sub <- scale(data_subset %>% dplyr::select(all_of(numeric_metrics_to_use)))
      
      set.seed(42)
      kmeans_sub <- kmeans(scaled_features_sub, centers = k_sub, nstart = 25)
      
      sub_labels <- paste0(parent_cluster, ".", kmeans_sub$cluster)
      
      app_state$processed_data$display_cluster[data_subset_indices] <- sub_labels
      
      parent_style <- app_state$behavior_styles %>% filter(cluster_id == parent_cluster)
      
      new_sub_styles <- tibble(
        cluster_id = sort(unique(sub_labels)),
        label = sort(unique(sub_labels)),
        color = colorRampPalette(c(darken_color(parent_style$color, 1.2), darken_color(parent_style$color, 0.5)))(k_sub)
      )
      
      app_state$behavior_styles <- app_state$behavior_styles %>%
        filter(cluster_id != parent_cluster) %>%
        bind_rows(new_sub_styles)
      
      incProgress(0.7, detail="Done.")
    })
    
    showNotification(paste(parent_cluster, "was successfully subdivided into", k_sub, "sub-behaviors."), type="message")
  })
  
  observeEvent(input$clear_sub_clustering, {
    req(app_state$processed_data)
    
    app_state$processed_data$display_cluster <- app_state$processed_data$kmeans_cluster
    
    k_val <- app_state$current_k
    default_styles_df <- tibble(
      cluster_id = c("B1", "B2", "B3", "B4", "B5"),
      label = c("B1: Sessile/Spherical", "B2: Confined Erratic", "B3: Directional/Migratory", "B4: Flattened/Dynamic", "B5: Adherent/Oscillatory"),
      color = c("#E8E654", "#5EC468", "#C780CF", "#DB5390", "#5DACF5")
    )
    if (k_val == 5) {
      original_styles <- default_styles_df
    } else {
      original_styles <- tibble(
        cluster_id = paste0("B", 1:k_val),
        label = paste0("Behavior ", 1:k_val),
        color = DISTINCT_COLORS[1:k_val]
      )
    }
    app_state$behavior_styles <- original_styles
    
    showNotification("All sub-clustering has been cleared.", type="message")
  })
  
  output$condition_merging_ui <- renderUI({
    req(app_state$processed_data)
    
    available_conditions <- app_state$condition_styles$condition
    unmerged_conditions <- setdiff(available_conditions, c(unlist(app_state$merged_conditions_map), "Projected Data"))
    
    tagList(
      selectizeInput("conditions_to_merge", "1. Select conditions to merge:",
                     choices = unmerged_conditions, multiple = TRUE),
      textInput("new_merged_condition_name", "2. Name for merged condition:"),
      actionButton("apply_condition_merge", "3. Create Merged Condition",
                   icon = icon("object-group"), class = "btn-primary", width="100%")
    )
  })
  
  output$merged_conditions_list_ui <- renderUI({
    req(app_state$merged_conditions_map)
    if (length(app_state$merged_conditions_map) == 0) return(NULL)
    
    tags$div(
      hr(),
      tags$h5("Active Merged Conditions:"),
      lapply(names(app_state$merged_conditions_map), function(merged_name) {
        tags$div(
          style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px; padding: 5px; background-color: #f8f9fa; border-radius: 4px;",
          tags$span(tags$strong(merged_name), " (", paste(app_state$merged_conditions_map[[merged_name]], collapse=", "), ")"),
          actionButton(paste0("unmerge_", make.names(merged_name)), "", icon = icon("times"),
                       class = "btn-danger btn-xs", style="height: 24px; width: 24px; padding: 1px;")
        )
      })
    )
  })
  
  observeEvent(input$apply_condition_merge, {
    req(input$conditions_to_merge, input$new_merged_condition_name)
    
    selected_conds <- input$conditions_to_merge
    new_name <- trimws(input$new_merged_condition_name)
    
    if (length(selected_conds) < 2) {
      showNotification("Please select at least two conditions to merge.", type="warning")
      return()
    }
    if (nchar(new_name) == 0) {
      showNotification("Please provide a name for the new merged condition.", type="warning")
      return()
    }
    if (new_name %in% app_state$condition_styles$condition) {
      showNotification("The new name is already an existing condition name. Please choose a unique name.", type="error")
      return()
    }
    
    app_state$merged_conditions_map[[new_name]] <- selected_conds
    
    app_state$processed_data <- app_state$processed_data %>%
      mutate(active_group = if_else(primary_group %in% selected_conds, new_name, active_group))
    
    new_color <- setdiff(DISTINCT_COLORS, app_state$condition_styles$color)[1]
    if (is.na(new_color)) new_color <- "#000000"
    
    app_state$condition_styles <- app_state$condition_styles %>%
      filter(!condition %in% selected_conds) %>%
      add_row(condition = new_name, color = new_color)
    
    updateSelectizeInput(session, "conditions_to_merge", selected = "")
    updateTextInput(session, "new_merged_condition_name", value = "")
    
    showNotification(paste0("Conditions merged into '", new_name, "'."), type="message")
  })
  
  observe({
    req(app_state$merged_conditions_map)
    merged_names <- names(app_state$merged_conditions_map)
    
    lapply(merged_names, function(m_name) {
      observeEvent(input[[paste0("unmerge_", make.names(m_name))]], {
        
        original_conds <- app_state$merged_conditions_map[[m_name]]
        
        app_state$processed_data <- app_state$processed_data %>%
          mutate(active_group = if_else(active_group == m_name, primary_group, active_group))
        
        original_styles <- app_state$raw_data %>%
          filter(primary_group %in% original_conds) %>%
          distinct(primary_group)
        
        app_state$condition_styles <- app_state$condition_styles %>%
          filter(condition != m_name)
        
        for(orig_cond in original_conds) {
          if (!orig_cond %in% app_state$condition_styles$condition) {
            new_color <- setdiff(DISTINCT_COLORS, app_state$condition_styles$color)[1]
            if(is.na(new_color)) new_color <- "#808080"
            app_state$condition_styles <- app_state$condition_styles %>%
              add_row(condition = orig_cond, color = new_color)
          }
        }
        
        app_state$merged_conditions_map[[m_name]] <- NULL
        
        showNotification(paste0("Unmerged '", m_name, "'."), type="message")
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    })
  })
  
  output$behavior_style_ui <- renderUI({
    req(app_state$behavior_styles)
    style_df <- app_state$behavior_styles
    ui_elements <- lapply(1:nrow(style_df), function(i) {
      cluster_id <- style_df$cluster_id[i]
      default_label <- style_df$label[i]
      default_color <- style_df$color[i]
      color_input_prefix <- paste0("style_color_", make.names(cluster_id))
      setupColorSyncObserver(input, session, color_input_prefix)
      fluidRow(
        column(6, textInput(inputId = paste0("style_label_", make.names(cluster_id)), label = paste0(cluster_id, " Label"), value = default_label)),
        column(6, syncedColorInput(id_prefix = color_input_prefix, label = paste0(cluster_id, " Color"), value = default_color))
      )
    })
    tagList(
      h4("Customize Behavior Display"),
      p("Edit the labels and colors used for clusters in plots.", style="font-size: 0.9em; color: grey;"),
      hr(),
      ui_elements,
      hr(),
      actionButton("apply_styles", "Apply Style Changes", icon=icon("paint-brush"), class="btn-info", width="100%")
    )
  })
  
  observeEvent(input$apply_styles, {
    req(app_state$behavior_styles)
    updated_styles <- app_state$behavior_styles
    for (i in 1:nrow(updated_styles)) {
      cluster_id <- updated_styles$cluster_id[i]
      label_input_id <- paste0("style_label_", make.names(cluster_id))
      color_input_id <- paste0("style_color_", make.names(cluster_id), "_hex")
      if (!is.null(input[[label_input_id]])) {
        updated_styles$label[i] <- input[[label_input_id]]
      }
      if (!is.null(input[[color_input_id]])) {
        updated_styles$color[i] <- input[[color_input_id]]
      }
    }
    app_state$behavior_styles <- updated_styles
    showNotification("Behavior styles updated.", type = "message")
  })
  
  output$condition_color_ui <- renderUI({
    req(app_state$condition_styles)
    style_df <- app_state$condition_styles
    ui_elements <- lapply(1:nrow(style_df), function(i) {
      condition_name <- style_df$condition[i]
      default_color <- style_df$color[i]
      base_id <- make.names(paste0("color_condition_", condition_name))
      setupColorSyncObserver(input, session, base_id)
      syncedColorInput(id_prefix = base_id, label = paste0("Color for: ", condition_name), value = default_color)
    })
    tagList(
      h4("Customize Condition Colors"),
      wellPanel(ui_elements, hr(), actionButton("apply_condition_colors", "Apply Condition Colors", icon=icon("palette"), class="btn-info", width="100%"))
    )
  })
  
  observeEvent(input$apply_condition_colors, {
    req(app_state$condition_styles)
    updated_styles <- app_state$condition_styles
    for (i in 1:nrow(updated_styles)) {
      condition_name <- updated_styles$condition[i]
      input_id <- make.names(paste0("color_condition_", condition_name, "_hex"))
      if (!is.null(input[[input_id]])) {
        updated_styles$color[i] <- input[[input_id]]
      }
    }
    app_state$condition_styles <- updated_styles
    showNotification("Condition colors updated.", type = "message")
  })
  
  output$cumulative_style_ui <- renderUI({
    req(app_state$cumulative_group_styles)
    styles_df <- app_state$cumulative_group_styles
    
    if (nrow(styles_df) == 0) return(NULL)
    
    ui_elements <- lapply(1:nrow(styles_df), function(i) {
      group_name <- styles_df$group[i]
      default_color <- styles_df$color[i]
      base_id <- make.names(paste0("color_cumulative_", group_name))
      setupColorSyncObserver(input, session, base_id)
      syncedColorInput(id_prefix = base_id, label = paste0("Color for: ", str_wrap(group_name, 30)), value = default_color)
    })
    
    tagList(
      h4("Customize Cumulative Group Colors"),
      p("These groups are created when 'Enable Cumulative (AND) Filtering' is active.", style="font-size: 0.9em; color: grey;"),
      wellPanel(
        ui_elements,
        hr(),
        actionButton("apply_cumulative_colors", "Apply Cumulative Colors", icon=icon("palette"), class="btn-info", width="100%")
      )
    )
  })
  
  observeEvent(input$apply_cumulative_colors, {
    req(app_state$cumulative_group_styles)
    updated_styles <- app_state$cumulative_group_styles
    
    for (i in 1:nrow(updated_styles)) {
      group_name <- updated_styles$group[i]
      input_id <- make.names(paste0("color_cumulative_", group_name, "_hex"))
      if (!is.null(input[[input_id]])) {
        updated_styles$color[i] <- input[[input_id]]
      }
    }
    app_state$cumulative_group_styles <- updated_styles
    showNotification("Cumulative group colors updated.", type = "message")
  })
  
  output$download_culled_data <- downloadHandler(
    filename = function() { paste0('ACME_Plotted_Data_', Sys.Date(), '.csv') },
    content = function(file) { req(app_state$processed_data); fwrite(app_state$processed_data, file) }
  )
  
  output$download_atlas <- downloadHandler(
    filename = function() { "Reference_Atlas.rds" },
    content = function(file) {
      req(app_state$analysis_atlas, app_state$behavior_styles, app_state$condition_styles)
      
      atlas_to_save <- app_state$analysis_atlas
      
      atlas_to_save$curated_data <- atlas_to_save$processed_data %>% filter(data_source == "Base")
      atlas_to_save$processed_data <- NULL
      
      visualization_styles <- list(
        behavior_palette = dynamic_behavior_palette(),
        condition_palette = dynamic_condition_palette(),
        behavior_labels = dynamic_behavior_labels()
      )
      
      atlas_to_save$metadata <- list(
        feature_list = app_state$analysis_atlas$feature_list,
        visualization_styles = visualization_styles
      )
      
      saveRDS(atlas_to_save, file)
    }
  )
  
  generate_heatmap_object <- function(data_for_heatmap) {
    clean_data <- data_for_heatmap %>% filter(!is.na(display_cluster) & !is.na(active_group))
    req(nrow(clean_data) > 0)
    
    dynamic_palette <- dynamic_behavior_palette()
    dynamic_labels <- dynamic_behavior_labels()
    features_present <- app_state$analysis_atlas$feature_list
    
    annotation_col_df <- data.frame(
      Behavior = clean_data$display_cluster,
      Condition = clean_data$active_group,
      DataSource = clean_data$data_source,
      row.names = clean_data$cell_id
    )
    
    custom_order <- input$heatmap_cluster_order
    annotation_col_df$Behavior <- factor(annotation_col_df$Behavior, levels = custom_order)
    col_order <- order(annotation_col_df$Behavior, annotation_col_df$DataSource, annotation_col_df$Condition)
    
    matrix_data <- clean_data %>% dplyr::select(all_of(features_present)) %>% as.matrix() %>% t()
    colnames(matrix_data) <- clean_data$cell_id
    
    matrix_ordered <- matrix_data[, col_order, drop = FALSE]
    annotation_ordered <- annotation_col_df[col_order, , drop = FALSE]
    
    cluster_color_map <- dynamic_palette
    condition_color_map <- dynamic_condition_palette()
    names(cluster_color_map) <- dynamic_labels[names(cluster_color_map)]
    annotation_ordered$Behavior <- recode(as.character(annotation_ordered$Behavior), !!!dynamic_labels)
    
    ann_colors <- list(
      Behavior = cluster_color_map,
      Condition = condition_color_map,
      DataSource = c("Base" = "grey90", "Projected" = "#e41a1c")
    )
    
    final_matrix <- t(scale(t(matrix_ordered)))
    final_matrix[final_matrix > 2] <- 2
    final_matrix[final_matrix < -2] <- -2
    final_matrix[is.na(final_matrix) | is.infinite(final_matrix)] <- 0
    
    display_rownames <- FEATURE_NAME_LOOKUP[rownames(final_matrix)]
    
    pheatmap(final_matrix, main = "Cell Feature Heatmap", scale = "none",
             color = viridis(100), breaks = seq(-2, 2, length.out = 101),
             annotation_col = annotation_ordered, annotation_colors = ann_colors,
             cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
             border_color = NA, silent = TRUE, labels_row = display_rownames)
  }
  
  observe({
    req(app_state$processed_data, input$heatmap_cluster_order, app_state$current_k, input$main_tabs == "heatmap_explorer")
    withProgress(message = 'Generating Global Heatmap...', value=0.2, {
      p_obj <- generate_heatmap_object(app_state$processed_data)
      app_state$heatmap_plot <- p_obj
    })
  })
  
  output$heatmap_plot <- renderPlot({
    if(is.null(app_state$heatmap_plot)) {
      return(ggplot() + theme_void() + annotate("text", x=1, y=1, label="Heatmap will generate automatically.", size=5))
    }
    grid.newpage()
    grid.draw(app_state$heatmap_plot$gtable)
  })
  
  filtered_heatmap_data <- reactive({
    req(app_state$processed_data, sidebar_state$filter_active_group)
    validate(
      need(
        length(sidebar_state$filter_active_group) > 0,
        "To view the filtered heatmap, please select one or more conditions in the 'Interactive Data Filter' tab."
      )
    )
    app_state$processed_data %>%
      filter(active_group %in% sidebar_state$filter_active_group)
  })
  
  filtered_heatmap_object <- reactive({
    req(filtered_heatmap_data(), input$heatmap_cluster_order)
    validate(
      need(
        nrow(filtered_heatmap_data()) > 0,
        "No data remains for the selected conditions after curation. Cannot generate heatmap."
      )
    )
    generate_heatmap_object(filtered_heatmap_data())
  })
  
  output$filtered_heatmap_plot <- renderPlot({
    plot_obj <- filtered_heatmap_object()
    req(plot_obj)
    grid.newpage()
    grid.draw(plot_obj$gtable)
  })
  
  umap_plot_object <- reactive({
    req(app_state$processed_data, input$umap_color, app_state$current_k, !is.null(input$umap_points_point_alpha), !is.null(input$umap_points_dot_size), !is.null(input$umap_density_padding_factor), !is.null(input$umap_layer_order))
    plot_data <- app_state$processed_data; color_col <- input$umap_color
    is_categorical <- color_col %in% c('active_group', 'display_cluster', 'Biological_Behavior')
    x_col <- "feature_umap_x"; y_col <- "feature_umap_y"; padding <- input$umap_density_padding_factor / 100
    
    full_x_range <- range(plot_data[[x_col]], na.rm=TRUE)
    full_y_range <- range(plot_data[[y_col]], na.rm=TRUE)
    
    x_pad <- (full_x_range[2] - full_x_range[1]) * padding
    y_pad <- (full_y_range[2] - full_y_range[1]) * padding
    padded_limits <- c(full_x_range[1] - x_pad, full_x_range[2] + x_pad, full_y_range[1] - y_pad, full_y_range[2] + y_pad)
    
    p <- ggplot() + theme_minimal() + theme(panel.grid = element_blank()) + labs(title="UMAP of Curated Data", x="UMAP 1", y="UMAP 2")
    
    base_data_for_density <- plot_data %>% filter(data_source == "Base")
    polygon_data <- NULL; if(is_categorical && isTRUE(input$umap_density_show_density)) { polygon_data <- calculate_density_polygons(base_data_for_density, x_col, y_col, color_col, input$umap_density_smoothing, padded_limits, normalize = isTRUE(input$umap_density_normalize_density)) }
    
    palette <- if(color_col == 'display_cluster') { dynamic_behavior_palette() } else if (color_col == 'Biological_Behavior') { BEHAVIOR_COLORS } else { dynamic_condition_palette() }
    layer_order <- input$umap_layer_order
    
    for (layer in layer_order) {
      p <- switch(layer,
                  "Glow" = {
                    if (isTRUE(input$umap_density_show_glow) && !is.null(polygon_data)) {
                      for(group_name in names(polygon_data)) {
                        polys <- polygon_data[[group_name]]; if(is.null(polys)) next;
                        group_color <- palette[as.character(group_name)]
                        glow_polys <- head(polys, input$umap_density_glow_extent)
                        for(poly in glow_polys) {
                          poly_df <- as.data.frame(poly)
                          outline_width <- input$umap_density_glow_outline_width %||% 0
                          if (outline_width > 0) {
                            p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=group_color, alpha=input$umap_density_glow_opacity, color=darken_color(group_color, 0.4), linewidth=outline_width)
                          } else {
                            p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=group_color, alpha=input$umap_density_glow_opacity, color=NA)
                          }
                        }
                      }
                    }; p
                  },
                  "Core" = {
                    if (isTRUE(input$umap_density_show_core) && !is.null(polygon_data)) {
                      for(group_name in names(polygon_data)) {
                        polys <- polygon_data[[group_name]]; if(is.null(polys)) next;
                        group_color <- palette[as.character(group_name)]
                        core_polys <- tail(polys, input$umap_density_core_size)
                        for(poly in core_polys) {
                          poly_df <- as.data.frame(poly)
                          outline_width <- input$umap_density_core_outline_width %||% 0
                          if (outline_width > 0) {
                            p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=darken_color(group_color, input$umap_density_darken_factor), alpha=input$umap_density_core_opacity, color=darken_color(group_color, 0.4), linewidth=outline_width)
                          } else {
                            p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=darken_color(group_color, input$umap_density_darken_factor), alpha=input$umap_density_core_opacity, color=NA)
                          }
                        }
                      }
                    }; p
                  },
                  "Dense Core" = {
                    if (isTRUE(input$umap_density_show_dense_core) && !is.null(polygon_data)) {
                      for(group_name in names(polygon_data)) {
                        polys <- polygon_data[[group_name]]; if(is.null(polys)) next;
                        group_color <- palette[as.character(group_name)]
                        dense_core_polys <- tail(polys, input$umap_density_dense_core_size)
                        for(poly in dense_core_polys) {
                          poly_df <- as.data.frame(poly)
                          outline_width <- input$umap_density_dense_core_outline_width %||% 0
                          if (outline_width > 0) {
                            p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=darken_color(group_color, input$umap_density_dense_core_darkening), alpha=input$umap_density_dense_core_opacity, color=darken_color(group_color, 0.2), linewidth=outline_width)
                          } else {
                            p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=darken_color(group_color, input$umap_density_dense_core_darkening), alpha=input$umap_density_dense_core_opacity, color=NA)
                          }
                        }
                      }
                    }; p
                  },
                  "Points" = {
                    if(is_categorical) {
                      p + geom_point(data=plot_data, aes(x=.data[[x_col]], y=.data[[y_col]], color = .data[[color_col]], shape = data_source), size = input$umap_points_dot_size, alpha = input$umap_points_point_alpha)
                    } else {
                      sorted_data <- plot_data %>% arrange(.data[[color_col]])
                      p + geom_point(data=sorted_data, aes(x=.data[[x_col]], y=.data[[y_col]], color = .data[[color_col]], shape = data_source), size = input$umap_points_dot_size, alpha = input$umap_points_point_alpha)
                    }
                  })
    }
    
    if(is_categorical) {
      custom_labels <- if(color_col == 'display_cluster') dynamic_behavior_labels() else waiver()
      p <- p + scale_color_manual(values = palette, name = "Group", labels = custom_labels, drop = FALSE)
    } else {
      q_range <- c(input$umap_quantile_lower / 100, input$umap_quantile_upper / 100)
      color_limits <- quantile(plot_data[[color_col]], probs = q_range, na.rm = TRUE)
      p <- p + scale_color_viridis_c(option = "D", limits = color_limits, oob = scales::squish, name = "Value")
    }
    
    p <- p + scale_shape_manual(name = "Data Source", values = c("Base" = 16, "Projected" = 17))
    
    p + coord_cartesian(xlim=padded_limits[1:2], ylim=padded_limits[3:4], expand = FALSE)
  })
  
  output$umap_plot <- renderPlot({ req(umap_plot_object()); umap_plot_object() })
  output$umap_data_table <- DT::renderDataTable({ req(app_state$processed_data); DT::datatable(app_state$processed_data, options = list(scrollX = TRUE, pageLength = 5), filter = 'top', rownames = FALSE) })
  
  filtered_data_interactive <- reactive({
    req(app_state$processed_data, !is.null(input$cumulative_filter), filtered_samples_reactive())
    
    data <- app_state$processed_data %>% filter(experiment_id %in% filtered_samples_reactive())
    
    if (isFALSE(input$cumulative_filter)) {
      data_to_filter <- data %>% mutate(highlight_group = "unselected")
      for (filter_col_name in rev(filter_variables)) {
        selected_choices <- input[[paste0("filter_", filter_col_name)]]
        if (!is.null(selected_choices)) {
          data_to_filter <- data_to_filter %>% mutate(highlight_group = if_else(.data[[filter_col_name]] %in% selected_choices, as.character(.data[[filter_col_name]]), highlight_group))
        }
      }
      return(data_to_filter)
    } else {
      active_filters <- list()
      for (var in filter_variables) {
        selected <- input[[paste0("filter_", var)]]
        if (!is.null(selected) && length(selected) > 0) {
          active_filters[[var]] <- selected
        }
      }
      if (length(active_filters) == 0) {
        return(data %>% mutate(highlight_group = "unselected"))
      }
      highlight_info <- data
      for (col in names(active_filters)) {
        highlight_info <- highlight_info %>% filter(.data[[col]] %in% active_filters[[col]])
      }
      if (nrow(highlight_info) > 0) {
        highlight_info <- highlight_info %>% unite("highlight_group", all_of(names(active_filters)), sep = " & ", remove = FALSE) %>% dplyr::select(cell_id, highlight_group)
        data %>% left_join(highlight_info, by = "cell_id") %>% mutate(highlight_group = replace_na(highlight_group, "unselected"))
      } else {
        data %>% mutate(highlight_group = "unselected")
      }
    }
  })
  
  observe({
    req(filtered_data_interactive())
    data <- filtered_data_interactive()
    
    cumulative_groups <- unique(data$highlight_group)
    cumulative_groups <- cumulative_groups[grepl(" & ", cumulative_groups)]
    
    current_styles <- app_state$cumulative_group_styles
    used_colors <- c(app_state$condition_styles$color, app_state$behavior_styles$color, current_styles$color)
    
    updated_styles <- current_styles %>% filter(group %in% cumulative_groups)
    
    new_groups <- setdiff(cumulative_groups, updated_styles$group)
    if (length(new_groups) > 0) {
      available_colors <- setdiff(DISTINCT_COLORS, used_colors)
      if (length(available_colors) < length(new_groups)) {
        available_colors <- rep(DISTINCT_COLORS, length.out = length(new_groups))
      }
      
      new_styles_df <- tibble(
        group = new_groups,
        color = available_colors[seq_len(length(new_groups))]
      )
      updated_styles <- bind_rows(updated_styles, new_styles_df)
    }
    
    if (!identical(app_state$cumulative_group_styles, updated_styles)) {
      app_state$cumulative_group_styles <- updated_styles
    }
  })
  
  filter_color_map_reactive <- reactive({
    plot_data <- filtered_data_interactive()
    highlight_groups <- unique(plot_data$highlight_group[plot_data$highlight_group != "unselected"])
    
    color_map <- c("unselected"="grey85")
    
    b_palette <- dynamic_behavior_palette()
    c_palette <- dynamic_condition_palette()
    
    cumulative_styles <- app_state$cumulative_group_styles
    cumulative_palette <- if (!is.null(cumulative_styles) && nrow(cumulative_styles) > 0) {
      setNames(cumulative_styles$color, cumulative_styles$group)
    } else {
      list()
    }
    
    full_palette <- c(b_palette, c_palette, cumulative_palette)
    
    groups_to_color <- intersect(highlight_groups, names(full_palette))
    if (length(groups_to_color) > 0) {
      color_map <- c(color_map, full_palette[groups_to_color])
    }
    
    return(color_map)
  })
  
  generate_filter_plot <- function(plot_df, main_title) {
    req(nrow(plot_df) > 0, !is.null(input$filter_points_dot_size),
        !is.null(input$filter_points_point_alpha),
        !is.null(input$filter_density_padding_factor), !is.null(input$filter_layer_order),
        !is.null(input$filter_grey_points_dot_size), !is.null(input$filter_grey_points_point_alpha),
        !is.null(input$filter_color_continuous))
    
    x_col <- "feature_umap_x"; y_col <- "feature_umap_y"; padding <- input$filter_density_padding_factor / 100
    x_range <- range(app_state$processed_data[[x_col]], na.rm=TRUE); x_pad <- (x_range[2] - x_range[1]) * padding
    y_range <- range(app_state$processed_data[[y_col]], na.rm=TRUE); y_pad <- (y_range[2] - y_range[1]) * padding
    padded_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad, y_range[1] - y_pad, y_range[2] + y_pad)
    
    highlighted_df <- plot_df %>% filter(highlight_group != "unselected")
    unselected_df <- plot_df %>% filter(highlight_group == "unselected")
    
    p <- ggplot() + theme_minimal(base_size = 14) + labs(title = main_title, x="UMAP 1", y="UMAP 2") + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
    
    p <- p + geom_point(data = unselected_df, aes(x=.data[[x_col]], y=.data[[y_col]], shape = data_source),
                        color = "grey85",
                        size = input$filter_grey_points_dot_size,
                        alpha = input$filter_grey_points_point_alpha)
    
    if (input$filter_color_continuous == 'categorical') {
      color_map <- filter_color_map_reactive()
      polygon_data <- NULL; if(isTRUE(input$filter_density_show_density)) { polygon_data <- calculate_density_polygons(highlighted_df, x_col, y_col, "highlight_group", input$filter_density_smoothing, padded_limits, normalize = isTRUE(input$filter_density_normalize_density)) }
      layer_order <- input$filter_layer_order
      
      for (layer in layer_order) {
        p <- switch(layer,
                    "Glow" = {
                      if (isTRUE(input$filter_density_show_glow) && !is.null(polygon_data)) {
                        for(group_name in names(polygon_data)) {
                          polys <- polygon_data[[group_name]]; if(is.null(polys)) next;
                          group_color <- color_map[as.character(group_name)];
                          fill_color <- if (isTRUE(input$filter_density_greyscale)) "grey70" else group_color;
                          outline_color <- if (isTRUE(input$filter_density_greyscale)) "grey40" else darken_color(group_color, 0.4);
                          glow_polys <- head(polys, input$filter_density_glow_extent);
                          for(poly in glow_polys) {
                            poly_df <- as.data.frame(poly)
                            outline_width <- input$filter_density_glow_outline_width %||% 0
                            if (outline_width > 0) {
                              p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=fill_color, alpha=input$filter_density_glow_opacity, color=outline_color, linewidth=outline_width)
                            } else {
                              p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=fill_color, alpha=input$filter_density_glow_opacity, color=NA)
                            }
                          }
                        }
                      }; p
                    },
                    "Core" = {
                      if (isTRUE(input$filter_density_show_core) && !is.null(polygon_data)) {
                        for(group_name in names(polygon_data)) {
                          polys <- polygon_data[[group_name]]; if(is.null(polys)) next;
                          group_color <- color_map[as.character(group_name)];
                          fill_color <- if (isTRUE(input$filter_density_greyscale)) "grey40" else darken_color(group_color, input$filter_density_darken_factor);
                          outline_color <- if (isTRUE(input$filter_density_greyscale)) "grey20" else darken_color(group_color, 0.4);
                          core_polys <- tail(polys, input$filter_density_core_size);
                          for(poly in core_polys) {
                            poly_df <- as.data.frame(poly)
                            outline_width <- input$filter_density_core_outline_width %||% 0
                            if (outline_width > 0) {
                              p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=fill_color, alpha=input$filter_density_core_opacity, color=outline_color, linewidth=outline_width)
                            } else {
                              p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=fill_color, alpha=input$filter_density_core_opacity, color=NA)
                            }
                          }
                        }
                      }; p
                    },
                    "Dense Core" = {
                      if (isTRUE(input$filter_density_show_dense_core) && !is.null(polygon_data)) {
                        for(group_name in names(polygon_data)) {
                          polys <- polygon_data[[group_name]]; if(is.null(polys)) next;
                          group_color <- color_map[as.character(group_name)];
                          fill_color <- if (isTRUE(input$filter_density_greyscale)) "grey10" else darken_color(group_color, input$filter_density_dense_core_darkening);
                          outline_color <- if (isTRUE(input$filter_density_greyscale)) "black" else darken_color(group_color, 0.2);
                          dense_core_polys <- tail(polys, input$filter_density_dense_core_size);
                          for(poly in dense_core_polys) {
                            poly_df <- as.data.frame(poly)
                            outline_width <- input$filter_density_dense_core_outline_width %||% 0
                            if (outline_width > 0) {
                              p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=fill_color, alpha=input$filter_density_dense_core_opacity, color=outline_color, linewidth=outline_width)
                            } else {
                              p <- p + geom_polygon(data=poly_df, aes(x=x,y=y), fill=fill_color, alpha=input$filter_density_dense_core_opacity, color=NA)
                            }
                          }
                        }
                      }; p
                    },
                    "Points" = { if(nrow(highlighted_df) > 0) {
                      p + geom_point(data = highlighted_df, aes(x=.data[[x_col]], y=.data[[y_col]], color=highlight_group, shape = data_source), size = input$filter_points_dot_size, alpha = input$filter_points_point_alpha) } else { p } })
      }
      
      all_groups <- names(color_map); dynamic_labels <- dynamic_behavior_labels(); label_map <- setNames(all_groups, all_groups); kmeans_groups_in_legend <- intersect(all_groups, names(dynamic_labels)); if(length(kmeans_groups_in_legend) > 0) label_map[kmeans_groups_in_legend] <- dynamic_labels[kmeans_groups_in_legend]
      
      p <- p + scale_color_manual(values = color_map, name = "Selected Groups", labels = label_map, drop = FALSE) +
        guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1))
      
    } else {
      req(input$filter_quantile_lower, input$filter_quantile_upper)
      color_col <- input$filter_color_continuous
      
      if (nrow(highlighted_df) > 0 && color_col %in% names(highlighted_df)) {
        sorted_highlighted_df <- highlighted_df %>% arrange(.data[[color_col]])
        
        p <- p + geom_point(data = sorted_highlighted_df,
                            aes(x=.data[[x_col]], y=.data[[y_col]], color = .data[[color_col]], shape = data_source),
                            size = input$filter_points_dot_size,
                            alpha = input$filter_points_point_alpha)
        
        q_range <- c(input$filter_quantile_lower / 100, input$filter_quantile_upper / 100)
        color_limits <- quantile(highlighted_df[[color_col]], probs = q_range, na.rm = TRUE)
        legend_title <- FEATURE_NAME_LOOKUP[color_col]
        
        p <- p + scale_color_viridis_c(option = "D", limits = color_limits, oob = scales::squish, name = legend_title)
      }
    }
    
    p <- p + scale_shape_manual(name = "Data Source", values = c("Base" = 16, "Projected" = 17), drop = FALSE)
    p + coord_cartesian(xlim=padded_limits[1:2], ylim=padded_limits[3:4], expand = FALSE)
  }
  
  filter_plot_object <- reactive({
    req(filtered_data_interactive())
    plot_df <- filtered_data_interactive()
    if (nrow(plot_df) == 0) return(ggplot() + theme_void())
    generate_filter_plot(plot_df, "Interactively Filtered UMAP (Merged View)")
  })
  
  output$filter_plot <- renderPlot({ req(filter_plot_object()); filter_plot_object() })
  
  output$individual_sample_plots_ui <- renderUI({
    req(app_state$processed_data, input$filter_active_group, filtered_samples_reactive())
    selected_conditions <- input$filter_active_group
    
    if (is.null(selected_conditions) || length(selected_conditions) == 0) {
      return(wellPanel(style = "text-align: center; color: #777; margin-top: 20px;",
                       p("No conditions selected in the sidebar."),
                       p("Please select one or more conditions to view individual sample plots.")))
    }
    
    samples_to_plot <- app_state$processed_data %>%
      filter(active_group %in% selected_conditions,
             experiment_id %in% filtered_samples_reactive()) %>%
      pull(experiment_id) %>%
      unique() %>%
      sort()
    
    if (length(samples_to_plot) == 0) {
      return(wellPanel(style = "text-align: center; color: #777; margin-top: 20px;",
                       p("No samples found for the selected condition(s) and sample filter.")))
    }
    
    plot_outputs <- lapply(samples_to_plot, function(sample_id) {
      plot_id <- make.names(paste0("individual_plot_", sample_id))
      column(width = 6, style='padding:15px;',
             plotOutput(plot_id, height = "400px")
      )
    })
    
    fluidRow(plot_outputs)
  })
  
  observe({
    req(app_state$processed_data, input$filter_active_group, filtered_data_interactive(), filtered_samples_reactive())
    
    selected_conditions <- input$filter_active_group
    if (is.null(selected_conditions) || length(selected_conditions) == 0) return()
    
    full_filtered_data <- filtered_data_interactive()
    
    samples_to_plot <- app_state$processed_data %>%
      filter(active_group %in% selected_conditions,
             experiment_id %in% filtered_samples_reactive()) %>%
      pull(experiment_id) %>%
      unique()
    
    for (sample_id in samples_to_plot) {
      local({
        current_sample_id <- sample_id
        plot_id <- make.names(paste0("individual_plot_", current_sample_id))
        
        output[[plot_id]] <- renderPlot({
          foreground_df <- full_filtered_data %>% filter(experiment_id == current_sample_id)
          
          background_df <- app_state$processed_data %>%
            filter(experiment_id != current_sample_id) %>%
            mutate(highlight_group = "unselected")
          
          plot_data_for_sample <- bind_rows(foreground_df, background_df)
          
          if (nrow(plot_data_for_sample) > 0) {
            generate_filter_plot(plot_data_for_sample, main_title = current_sample_id)
          } else {
            ggplot() + theme_void() +
              labs(title = current_sample_id) +
              annotate("text", x = 1, y = 1, label = "No data points after filtering.")
          }
        })
      })
    }
  })
  
  output$download_umap_pdf <- downloadHandler(filename = function() { paste0("ACME_UMAP_Plot_", Sys.Date(), ".pdf") }, content = function(file) { plot_size <- input$umap_plot_size; pdf_width_in <- if (!is.null(plot_size$width)) plot_size$width / 96 else 11; pdf_height_in <- if (!is.null(plot_size$height)) plot_size$height / 96 else 8.5; ggsave(file, plot = umap_plot_object(), width = pdf_width_in, height = pdf_height_in, device = "pdf", units = "in", limitsize = FALSE) })
  output$download_heatmap_pdf <- downloadHandler(filename = function() { paste0("ACME_Global_Heatmap_Plot_", Sys.Date(), ".pdf") }, content = function(file) { req(app_state$heatmap_plot); plot_size <- input$heatmap_plot_size; pdf_width_in <- if (!is.null(plot_size$width)) plot_size$width / 96 else 10; pdf_height_in <- if (!is.null(plot_size$height)) plot_size$height / 96 else 14; ggsave(file, plot = app_state$heatmap_plot$gtable, width = pdf_width_in, height = pdf_height_in, device = "pdf", units = "in", limitsize = FALSE) })
  output$download_filter_plot_pdf <- downloadHandler(filename = function() { paste0("ACME_Filter_Plot_", Sys.Date(), ".pdf") }, content = function(file) { plot_size <- input$filter_plot_size; pdf_width_in <- if (!is.null(plot_size$width)) plot_size$width / 96 else 11; pdf_height_in <- if (!is.null(plot_size$height)) plot_size$height / 96 else 8.5; ggsave(file, plot = filter_plot_object(), width = pdf_width_in, height = pdf_height_in, device = "pdf", units = "in", limitsize = FALSE) })
  
  output$download_filtered_heatmap_pdf <- downloadHandler(
    filename = function() { paste0("ACME_Filtered_Heatmap_Plot_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(filtered_heatmap_object())
      plot_obj <- filtered_heatmap_object()
      plot_size <- input$filtered_heatmap_plot_size
      pdf_width_in <- if (!is.null(plot_size$width)) plot_size$width / 96 else 10
      pdf_height_in <- if (!is.null(plot_size$height)) plot_size$height / 96 else 14
      ggsave(file, plot = plot_obj$gtable, width = pdf_width_in, height = pdf_height_in, device = "pdf", units = "in", limitsize = FALSE)
    }
  )
  
  output$download_individual_plots_pdf <- downloadHandler(
    filename = function() {
      paste0("ACME_Individual_Sample_UMAPs_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(app_state$processed_data, input$filter_active_group, filtered_data_interactive(), filtered_samples_reactive())
      
      selected_conditions <- input$filter_active_group
      full_filtered_data <- filtered_data_interactive()
      
      samples_to_plot <- app_state$processed_data %>%
        filter(active_group %in% selected_conditions,
               experiment_id %in% filtered_samples_reactive()) %>%
        pull(experiment_id) %>%
        unique() %>%
        sort()
      
      if (length(samples_to_plot) == 0) {
        pdf(file, width = 8, height = 5)
        plot.new()
        text(0.5, 0.5, "No samples to plot for selected conditions.")
        dev.off()
        return()
      }
      
      plot_list <- lapply(samples_to_plot, function(sample_id) {
        foreground_df <- full_filtered_data %>% filter(experiment_id == sample_id)
        background_df <- app_state$processed_data %>%
          filter(experiment_id != sample_id) %>%
          mutate(highlight_group = "unselected")
        plot_data_for_sample <- bind_rows(foreground_df, background_df)
        
        if (nrow(plot_data_for_sample) > 0) {
          generate_filter_plot(plot_data_for_sample, main_title = sample_id)
        } else {
          ggplot() + theme_void() + labs(title = sample_id) +
            annotate("text", x = 1, y = 1, label = "No data points after filtering.")
        }
      })
      
      ncol <- 2
      nrow <- ceiling(length(plot_list) / ncol)
      pdf_width <- 8.5 * ncol
      pdf_height <- 7 * nrow
      
      arranged_plots <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow)
      
      ggsave(file, plot = arranged_plots, width = pdf_width, height = pdf_height, device = "pdf", units = "in", limitsize = FALSE)
    }
  )
  
  output$download_behavior_dist_pdf <- downloadHandler(
    filename = function() { paste0("ACME_Behavior_Distribution_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(app_state$behavior_pie_plots)
      plots_to_save <- app_state$behavior_pie_plots
      if (length(plots_to_save) == 0) {
        pdf(file, width = 8, height = 5)
        plot.new(); text(0.5, 0.5, "No plots to display."); dev.off()
        return()
      }
      
      dist_level <- isolate(input$dist_summary_level)
      ncol <- if (dist_level == "Group by Condition") 2 else 3
      nrow <- ceiling(length(plots_to_save) / ncol)
      pdf_width <- 6 * ncol
      pdf_height <- 5 * nrow
      arranged_plots <- ggpubr::ggarrange(plotlist = plots_to_save, ncol = ncol, nrow = nrow, common.legend = TRUE, legend = "bottom")
      
      pdf(file, width = pdf_width, height = pdf_height)
      print(arranged_plots)
      dev.off()
    }
  )
  
  output$download_condition_dist_pdf <- downloadHandler(
    filename = function() { paste0("ACME_Condition_Distribution_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(app_state$condition_in_behavior_pie_plots)
      plots_to_save <- app_state$condition_in_behavior_pie_plots
      if (length(plots_to_save) == 0) {
        pdf(file, width = 8, height = 5)
        plot.new(); text(0.5, 0.5, "No plots to display."); dev.off()
        return()
      }
      ncol <- 2
      nrow <- ceiling(length(plots_to_save) / ncol)
      pdf_width <- 6 * ncol
      pdf_height <- 5.5 * nrow
      arranged_plots <- ggpubr::ggarrange(plotlist = plots_to_save, ncol = ncol, nrow = nrow, common.legend = TRUE, legend = "bottom")
      
      pdf(file, width = pdf_width, height = pdf_height)
      print(arranged_plots)
      dev.off()
    }
  )
  
  output$behavior_pie_charts_ui <- renderUI({
    req(input$main_tabs == "interactive_filter", input$dist_summary_level, filtered_samples_reactive())
    selected_conditions <- input$filter_active_group
    if (is.null(selected_conditions) || length(selected_conditions) == 0) { return(wellPanel(style = "text-align: center; color: #777; margin-top: 20px;", p("No conditions selected."), p("Select conditions in the sidebar to view summary charts."))) }
    
    if (input$dist_summary_level == "Group by Condition") {
      groups_to_plot <- selected_conditions
      col_width <- 6
    } else {
      req(app_state$processed_data)
      groups_to_plot <- app_state$processed_data %>%
        filter(active_group %in% selected_conditions, experiment_id %in% filtered_samples_reactive()) %>%
        pull(experiment_id) %>%
        unique() %>%
        sort()
      col_width <- 4
    }
    
    if (length(groups_to_plot) == 0) { return(wellPanel(style = "text-align: center; color: #777; margin-top: 20px;", p("No samples found for the selected filters."))) }
    plot_outputs <- lapply(groups_to_plot, function(group_name) { column(width = col_width, style='padding:15px;', plotOutput(paste0("pie_behavior_", make.names(group_name)), height = "280px")) }); fluidRow(plot_outputs)
  })
  
  observe({
    req(input$main_tabs == "interactive_filter", app_state$processed_data, !is.null(input$dist_plot_mode), !is.null(input$dist_summary_level), filtered_samples_reactive())
    selected_conditions <- input$filter_active_group; if (is.null(selected_conditions) || length(selected_conditions) == 0) { app_state$behavior_pie_plots <- list(); return() }
    
    plot_data <- app_state$processed_data %>% filter(active_group %in% selected_conditions, experiment_id %in% filtered_samples_reactive())
    
    if (input$dist_summary_level == "Group by Condition") {
      groups_to_iterate <- selected_conditions
      grouping_var <- "active_group"
    } else {
      groups_to_iterate <- unique(plot_data$experiment_id)
      grouping_var <- "experiment_id"
    }
    generated_plots <- list()
    for (group_val in groups_to_iterate) {
      
      p <- tryCatch({
        group_data <- plot_data %>% filter(!!sym(grouping_var) == group_val & !is.na(display_cluster))
        
        if (nrow(group_data) == 0) {
          ggplot() + theme_void() + labs(title = group_val) + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
        } else if (input$dist_plot_mode == "Relative Percentage (Donut Chart)") {
          summary_data <- group_data %>% count(display_cluster, name = "n_cells")
          hole_size <- input$donut_hole_size %||% 0
          summary_data <- summary_data %>%
            mutate(
              proportion = n_cells / sum(n_cells),
              Label = paste0(round(100 * proportion, 1), "%")
            ) %>%
            arrange(desc(display_cluster)) %>%
            mutate(
              ymax = cumsum(proportion),
              ymin = lag(ymax, default = 0),
              label_pos = (ymax + ymin) / 2
            )
          outer_radius <- 2; inner_radius <- hole_size * 2
          ggplot(summary_data, aes(ymax = ymax, ymin = ymin, xmax = outer_radius, xmin = inner_radius, fill = display_cluster)) +
            geom_rect(color = "white") +
            coord_polar(theta = "y", start = 0) +
            geom_label_repel(aes(y = label_pos, x = (outer_radius + inner_radius) / 2, label = Label), show.legend = FALSE, max.overlaps = Inf, size = 3.5) +
            scale_fill_manual(values = dynamic_behavior_palette(), name = "Cluster", drop = FALSE, labels = dynamic_behavior_labels()) +
            theme_void() + labs(title = group_val) +
            theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), legend.position = "bottom") +
            xlim(0, outer_radius + 0.5)
        } else { # Frequency (Bar Chart)
          summary_data <- group_data %>%
            count(display_cluster, name = "n_cells") %>%
            mutate(freq = n_cells / sum(n_cells)) %>%
            left_join(app_state$behavior_styles, by = c("display_cluster" = "cluster_id"))
          
          ggplot(summary_data, aes(x = label, y = freq, fill = display_cluster)) +
            geom_col(show.legend = FALSE) +
            geom_text(aes(label = scales::percent(freq, accuracy = 0.1)), vjust = -0.3) +
            scale_fill_manual(values = dynamic_behavior_palette(), drop = FALSE) +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
            scale_y_continuous(labels = scales::percent_format()) +
            theme_minimal() +
            labs(title = group_val, x = NULL, y = "Frequency") +
            theme(
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank()
            )
        }
      }, error = function(e) {
        ggplot() + theme_void() + labs(title = group_val) + annotate("text", x=1, y=1, label="Plotting Error", color="red")
      })
      generated_plots[[make.names(group_val)]] <- p
    }
    app_state$behavior_pie_plots <- generated_plots
    for (group_val in groups_to_iterate) { local({ my_group_name <- make.names(group_val); output[[paste0("pie_behavior_", my_group_name)]] <- renderPlot({ req(app_state$behavior_pie_plots[[my_group_name]]); print(app_state$behavior_pie_plots[[my_group_name]]) }) }) }
  })
  
  output$behavior_condition_dist_ui <- renderUI({
    req(app_state$processed_data, app_state$behavior_styles); selected_conditions <- input$filter_active_group
    if (is.null(selected_conditions) || length(selected_conditions) == 0) { return(wellPanel(style = "text-align: center; color: #777; margin-top: 20px;", p("No conditions selected."), p("Select one or more conditions in the sidebar to view this analysis."))) }
    all_behaviors <- app_state$behavior_styles$cluster_id
    plot_outputs <- lapply(all_behaviors, function(b_id) { plotOutput(paste0("pie_cond_in_", make.names(b_id)), height="300px") })
    k_val <- length(all_behaviors)
    rows <- tagList(); for (i in seq(1, k_val, 2)) { row_plots <- list(); if (i <= k_val) row_plots[[1]] <- column(6, style='padding:15px;', plot_outputs[[i]]); if ((i+1) <= k_val) row_plots[[2]] <- column(6, style='padding:15px;', plot_outputs[[i+1]]); rows <- tagList(rows, fluidRow(row_plots)) }; return(rows)
  })
  
  observe({
    req(input$main_tabs == "interactive_filter", app_state$processed_data, app_state$behavior_styles, !is.null(input$dist_plot_mode), filtered_samples_reactive())
    selected_conditions <- input$filter_active_group; if (is.null(selected_conditions) || length(selected_conditions) == 0) { app_state$condition_in_behavior_pie_plots <- list(); return() }
    
    data_subset <- app_state$processed_data %>% filter(active_group %in% selected_conditions, experiment_id %in% filtered_samples_reactive())
    
    all_behaviors <- app_state$behavior_styles$cluster_id
    generated_plots <- list()
    for (behavior_id in all_behaviors) {
      title_label <- dynamic_behavior_labels()[behavior_id]
      p <- tryCatch({
        summary_data <- data_subset %>% filter(display_cluster == behavior_id) %>% count(active_group, name = "n_cells")
        if (nrow(summary_data) == 0) {
          ggplot() + theme_void() + labs(title = title_label) + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        } else if (input$dist_plot_mode == "Relative Percentage (Donut Chart)") {
          hole_size <- input$donut_hole_size %||% 0
          summary_data <- summary_data %>%
            mutate(
              proportion = n_cells / sum(n_cells),
              Label = paste0(round(100 * proportion, 1), "%")
            ) %>%
            arrange(desc(active_group)) %>%
            mutate(
              ymax = cumsum(proportion),
              ymin = lag(ymax, default = 0),
              label_pos = (ymax + ymin) / 2
            )
          outer_radius <- 2; inner_radius <- hole_size * 2
          ggplot(summary_data, aes(ymax = ymax, ymin = ymin, xmax = outer_radius, xmin = inner_radius, fill = active_group)) +
            geom_rect(color = "white") +
            coord_polar(theta = "y", start = 0) +
            geom_label_repel(aes(y = label_pos, x = (outer_radius + inner_radius) / 2, label = Label), show.legend = FALSE, max.overlaps = Inf, size = 3.5) +
            scale_fill_manual(values = dynamic_condition_palette(), name = "Condition", drop = FALSE) +
            theme_void() +
            labs(title = title_label) +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), legend.position = "bottom") +
            xlim(0, outer_radius + 0.5)
        } else { # Frequency (Bar Chart)
          summary_data <- summary_data %>% mutate(freq = n_cells / sum(n_cells))
          ggplot(summary_data, aes(x = active_group, y = freq, fill = active_group)) +
            geom_col(show.legend = FALSE) +
            geom_text(aes(label = scales::percent(freq, accuracy = 0.1)), vjust = -0.3) +
            scale_fill_manual(values = dynamic_condition_palette(), drop = FALSE) +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
            scale_y_continuous(labels = scales::percent_format()) +
            theme_minimal() +
            labs(title = title_label, x = "Condition", y = "Frequency") +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank()
            )
        }
      }, error = function(e) {
        ggplot() + theme_void() + labs(title = title_label) + annotate("text", x=1, y=1, label="Plotting Error", color="red")
      })
      generated_plots[[make.names(behavior_id)]] <- p
    }
    app_state$condition_in_behavior_pie_plots <- generated_plots
    for (behavior_id in all_behaviors) { local({ my_b_id <- make.names(behavior_id); output[[paste0("pie_cond_in_", my_b_id)]] <- renderPlot({ req(app_state$condition_in_behavior_pie_plots[[my_b_id]]); print(app_state$condition_in_behavior_pie_plots[[my_b_id]]) }) }) }
  })
}

# --- 6. Run Application ---
shinyApp(ui = ui, server = server)