library(shiny)
library(tidyverse)
library(Seurat)

jaccard_similarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}

standardize_bootstrap_cols <- function(df) {
  nm <- names(df)
  nm_std <- nm %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")

  names(df) <- nm_std

  if ("maxjaccard" %in% names(df) && !("max_jaccard" %in% names(df))) {
    df <- dplyr::rename(df, max_jaccard = maxjaccard)
  }

  if ("boostrap_number" %in% names(df) && !("bootstrap_number" %in% names(df))) {
    df <- dplyr::rename(df, bootstrap_number = boostrap_number)
  }

  if ("boostrapnumber" %in% names(df) && !("bootstrap_number" %in% names(df))) {
    df <- dplyr::rename(df, bootstrap_number = boostrapnumber)
  }

  if ("bootstrapnumber" %in% names(df) && !("bootstrap_number" %in% names(df))) {
    df <- dplyr::rename(df, bootstrap_number = bootstrapnumber)
  }

  required <- c("clusterid", "max_jaccard", "bootstrap_number")
  missing <- setdiff(required, names(df))

  if (length(missing) > 0) {
    stop(
      "Bootstrap TSV is missing required columns after normalization: ",
      paste(missing, collapse = ", ")
    )
  }

  df
}

MakeInterVsIntraStablePlot <- function(meta1, meta2,
                                       bootstraps1, bootstraps2,
                                       threshold, label1, label2) {

  clusters1 <- tibble::tibble(
    cellbarcode = rownames(meta1),
    clusterid1 = as.character(meta1$seurat_clusters)
  )

  clusters2 <- tibble::tibble(
    cellbarcode = rownames(meta2),
    clusterid2 = as.character(meta2$seurat_clusters)
  )

  clusters_merged <- dplyr::full_join(clusters1, clusters2, by = "cellbarcode")

  indices1 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid1)
  indices2 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid2)

  indices1 <- indices1[!is.na(names(indices1))]
  indices2 <- indices2[!is.na(names(indices2))]

  jaccard_list <- list()
  for (i in names(indices1)) {
    for (j in names(indices2)) {
      set1 <- indices1[[i]]
      set2 <- indices2[[j]]
      similarity <- jaccard_similarity(set1, set2)
      jaccard_list[[length(jaccard_list) + 1]] <- list(
        name1 = i,
        name2 = j,
        jaccard_similarity = similarity
      )
    }
  }

  jaccard_df <- do.call(rbind, lapply(jaccard_list, as.data.frame))
  jaccard_df <- type.convert(jaccard_df, as.is = TRUE)

  jaccard_df$name1 <- factor(jaccard_df$name1,
                             levels = sort(unique(jaccard_df$name1)))
  jaccard_df$name2 <- factor(jaccard_df$name2,
                             levels = sort(unique(jaccard_df$name2)))
  jaccard_tibble <- tibble::as_tibble(jaccard_df)

  stable_cluster_count <- jaccard_tibble %>%
    dplyr::filter(jaccard_similarity >= threshold) %>%
    dplyr::summarise(n = dplyr::n_distinct(name1)) %>%
    dplyr::pull(n)

  bootstrap1_summary <- bootstraps1 %>%
    dplyr::group_by(bootstrap_number) %>%
    dplyr::filter(max_jaccard >= threshold) %>%
    dplyr::summarise(n_clusters = dplyr::n_distinct(clusterid), .groups = "drop") %>%
    dplyr::mutate(method = "n_stable_name1")

  bootstrap2_summary <- bootstraps2 %>%
    dplyr::group_by(bootstrap_number) %>%
    dplyr::filter(max_jaccard >= threshold) %>%
    dplyr::summarise(n_clusters = dplyr::n_distinct(clusterid), .groups = "drop") %>%
    dplyr::mutate(method = "n_stable_name2")

  boot_stable_merged <- dplyr::bind_rows(bootstrap1_summary, bootstrap2_summary) %>%
    dplyr::mutate(method = factor(method,
                                  levels = c("n_stable_name1", "n_stable_name2")))

  method_colors <- c(
    n_stable_name1 = "forestgreen",
    n_stable_name2 = "dodgerblue3"
  )

  max_y <- max(c(boot_stable_merged$n_clusters, stable_cluster_count), na.rm = TRUE)
  upper_y <- max(1, max_y) * 1.15
  y_break_step <- max(1, floor(max_y / 8))
  label_y <- -max(1, max_y) * 0.08
  label_data <- tibble::tibble(
    method = factor(c("n_stable_name1", "n_stable_name2"),
                    levels = c("n_stable_name1", "n_stable_name2")),
    x_pos = c(1.08, 1.92),
    y_pos = c(label_y * 0.85, label_y * 1.25),
    label = c(label1, label2)
  )

  boot_stable_merged %>%
    ggplot2::ggplot(ggplot2::aes(x = method, y = n_clusters, fill = method, color = method)) +
    ggplot2::geom_hline(yintercept = -0.5, color = "black", linewidth = 0.5) +
    ggplot2::geom_violin(width = 0.7, trim = FALSE, alpha = 0.35, linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = stable_cluster_count,
                        color = "black", linetype = "dashed", linewidth = 1) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = stable_cluster_count,
      label = paste0(
        "number of inter-method stable clusters, min. Jaccard similarity threshold = ",
        as.character(threshold)
      ),
      vjust = -0.6,
      size = 5.2,
      color = "black"
    ) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = x_pos, y = y_pos, label = label, color = method),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 1,
      size = 3.4
    ) +
    ggplot2::scale_x_discrete(
      labels = NULL,
      expand = c(0.20, 0.20)
    ) +
    ggplot2::scale_fill_manual(values = method_colors, guide = "none") +
    ggplot2::scale_color_manual(values = method_colors, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max_y, by = y_break_step),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.5, upper_y), clip = "off") +
    ggplot2::labs(
      x = "",
      y = "# stable clusters"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 46, 10)
    )
}

find_method_pairs <- function(data_dir = "data") {
  rds_paths <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

  if (length(rds_paths) == 0) {
    return(tibble::tibble(label = character(), rds = character(), boot = character()))
  }

  tibble::tibble(rds = rds_paths) %>%
    dplyr::mutate(
      boot = purrr::map_chr(rds, infer_bootstrap_path),
      exists_boot = file.exists(boot),
      label = basename(rds)
    ) %>%
    dplyr::filter(exists_boot) %>%
    dplyr::select(label, rds, boot) %>%
    dplyr::arrange(label)
}

infer_bootstrap_path <- function(rds_path) {
  base_no_ext <- stringr::str_replace(rds_path, "\\.rds$", "")
  candidates <- c(
    paste0(base_no_ext, "_bootstraps.tsv"),
    paste0(base_no_ext, "_clusterbootstraps.tsv")
  )

  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    return(candidates[[1]])
  }
  existing[[1]]
}

read_seurat_meta <- function(rds_path) {
  obj <- readRDS(rds_path)

  if (!inherits(obj, "Seurat")) {
    stop("RDS is not a Seurat object: ", basename(rds_path))
  }

  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("Missing `seurat_clusters` in Seurat metadata: ", basename(rds_path))
  }

  obj@meta.data
}

read_bootstrap_tsv <- function(tsv_path) {
  readr::read_tsv(tsv_path, show_col_types = FALSE) %>%
    standardize_bootstrap_cols()
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Make method labels responsive to viewport width */
      .selectize-control.single .selectize-input,
      .selectize-dropdown .option {
        font-size: clamp(9px, 0.95vw, 12px);
      }
      .selectize-control.single .selectize-input > div.item {
        font-size: clamp(9px, 0.95vw, 12px);
      }
    "))
  ),
  titlePanel("Inter-method vs bootstrap cluster stability"),
  sidebarLayout(
    sidebarPanel(
      width = 5,
      helpText("Select two methods (.rds) with matching _bootstraps.tsv files."),
      selectInput("method1", "Method 1", choices = NULL),
      selectInput("method2", "Method 2", choices = NULL),
      numericInput("threshold", "Minimum Jaccard threshold", value = 0.6, min = 0, max = 1, step = 0.01),
      actionButton("refresh", "Refresh file list")
    ),
    mainPanel(
      width = 7,
      plotOutput("stability_plot", height = "500px"),
      verbatimTextOutput("status")
    )
  )
)

server <- function(input, output, session) {
  method_pairs <- reactiveVal(find_method_pairs("data"))

  refresh_choices <- function() {
    pairs <- find_method_pairs("data")
    method_pairs(pairs)

    if (nrow(pairs) > 0) {
      choice_map <- stats::setNames(pairs$rds, pairs$label)
      updateSelectInput(session, "method1", choices = choice_map, selected = pairs$rds[[1]])
      second_idx <- if (nrow(pairs) >= 2) 2 else 1
      updateSelectInput(session, "method2", choices = choice_map, selected = pairs$rds[[second_idx]])
    } else {
      updateSelectInput(session, "method1", choices = c())
      updateSelectInput(session, "method2", choices = c())
    }
  }

  observeEvent(input$refresh, refresh_choices(), ignoreInit = TRUE)
  observe(refresh_choices())

  loaded_data <- reactive({
    req(input$method1, input$method2)

    validate(
      need(input$method1 != input$method2, "Pick two different methods."),
      need(file.exists(input$method1), "Method 1 file does not exist."),
      need(file.exists(input$method2), "Method 2 file does not exist.")
    )

    boot1 <- infer_bootstrap_path(input$method1)
    boot2 <- infer_bootstrap_path(input$method2)

    validate(
      need(file.exists(boot1), paste("Missing bootstrap TSV for method 1:", basename(boot1))),
      need(file.exists(boot2), paste("Missing bootstrap TSV for method 2:", basename(boot2)))
    )

    list(
      meta1 = read_seurat_meta(input$method1),
      meta2 = read_seurat_meta(input$method2),
      b1 = read_bootstrap_tsv(boot1),
      b2 = read_bootstrap_tsv(boot2),
      l1 = stringr::str_remove(basename(input$method1), "\\.rds$"),
      l2 = stringr::str_remove(basename(input$method2), "\\.rds$")
    )
  })

  output$stability_plot <- renderPlot({
    dat <- loaded_data()

    MakeInterVsIntraStablePlot(
      meta1 = dat$meta1,
      meta2 = dat$meta2,
      bootstraps1 = dat$b1,
      bootstraps2 = dat$b2,
      threshold = input$threshold,
      label1 = dat$l1,
      label2 = dat$l2
    )
  })

  output$status <- renderText({
    pairs <- method_pairs()

    if (nrow(pairs) == 0) {
      return("No valid method pairs found in ./data. Expected .rds plus matching _bootstraps.tsv file.")
    }

    paste0("Detected ", nrow(pairs), " method file pairs in ./data")
  })
}

shinyApp(ui = ui, server = server)
