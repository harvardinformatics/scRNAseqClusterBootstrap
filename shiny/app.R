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
    levels = sort(unique(jaccard_df$name1))
  )
  jaccard_df$name2 <- factor(jaccard_df$name2,
    levels = sort(unique(jaccard_df$name2))
  )
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
      levels = c("n_stable_name1", "n_stable_name2")
    ))

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
      levels = c("n_stable_name1", "n_stable_name2")
    ),
    x_pos = c(1.08, 1.92),
    y_pos = c(label_y * 0.85, label_y * 1.25),
    label = c(label1, label2)
  )

  boot_stable_merged %>%
    ggplot2::ggplot(ggplot2::aes(x = method, y = n_clusters, fill = method, color = method)) +
    ggplot2::geom_hline(yintercept = -0.5, color = "black", linewidth = 0.5) +
    ggplot2::geom_violin(width = 0.7, trim = FALSE, alpha = 0.35, linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = stable_cluster_count,
      color = "black", linetype = "dashed", linewidth = 1
    ) +
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

extract_feature_table <- function(obj, assay_name) {
  assay_obj <- obj[[assay_name]]
  feature_ids <- rownames(SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts"))
  assay_symbols <- row.names(assay_obj)

  gene_symbols <- if (
    !is.null(assay_symbols) &&
    length(assay_symbols) == length(feature_ids)
  ) {
    assay_symbols
  } else {
    feature_ids
  }

  tibble::tibble(
    feature_id = feature_ids,
    gene_symbol = gene_symbols
  )
}

select_expression_assay <- function(obj) {
  assay_names <- names(obj@assays)

  if ("RNA" %in% assay_names) {
    return("RNA")
  }

  SeuratObject::DefaultAssay(obj)
}

read_method_bundle <- function(rds_path, boot_path) {
  obj <- readRDS(rds_path)

  if (!inherits(obj, "Seurat")) {
    stop("RDS is not a Seurat object: ", basename(rds_path))
  }

  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("Missing `seurat_clusters` in Seurat metadata: ", basename(rds_path))
  }

  assay_name <- select_expression_assay(obj)
  counts <- SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts")
  lib_size <- Matrix::colSums(counts)
  feature_table <- extract_feature_table(obj, assay_name)

  list(
    label = stringr::str_remove(basename(rds_path), "\\.rds$"),
    rds = rds_path,
    boot = boot_path,
    assay = assay_name,
    meta = obj@meta.data,
    bootstraps = read_bootstrap_tsv(boot_path),
    counts = counts,
    lib_size = lib_size,
    barcodes = colnames(counts),
    features = rownames(counts),
    feature_table = feature_table
  )
}

resolve_feature_id <- function(method, gene_symbol) {
  match_tbl <- method$feature_table %>%
    dplyr::filter(.data$gene_symbol == .env$gene_symbol)

  if (nrow(match_tbl) == 0) {
    return(NULL)
  }

  match_tbl$feature_id[[1]]
}

make_expression_heatmap_data <- function(method_data, gene_symbol) {
  barcode_union <- sort(unique(unlist(purrr::map(method_data, "barcodes"))))

  purrr::map_dfr(method_data, function(method) {
    values <- rep(NA_real_, length(barcode_union))
    names(values) <- barcode_union
    barcode_present <- barcode_union %in% method$barcodes
    cell_status <- ifelse(barcode_present, "expression_missing", "barcode_missing")

    feature_id <- resolve_feature_id(method, gene_symbol)

    if (!is.null(feature_id)) {
      gene_counts <- as.numeric(method$counts[feature_id, , drop = TRUE])
      normalized <- log1p((gene_counts / pmax(method$lib_size, 1)) * 10000)
      values[method$barcodes] <- normalized
      cell_status[barcode_present] <- "expression_present"
    }

    tibble::tibble(
      method = method$label,
      barcode = factor(barcode_union, levels = barcode_union),
      expression = values,
      cell_status = factor(
        cell_status,
        levels = c("barcode_missing", "expression_missing", "expression_present")
      )
    )
  }) %>%
    dplyr::mutate(
      method = factor(method, levels = purrr::map_chr(method_data, "label"))
    )
}

make_expression_heatmap_plot <- function(plot_data, gene_symbol) {
  legend_data <- tibble::tibble(
    status_label = c("Barcode missing", "Expression missing/undefined"),
    x = levels(plot_data$barcode)[1],
    y = levels(plot_data$method)[1]
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = barcode, y = method, fill = expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      na.value = "transparent",
      name = "log1p(norm counts)"
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(plot_data, cell_status == "barcode_missing"),
      fill = "white",
      inherit.aes = FALSE,
      ggplot2::aes(x = barcode, y = method)
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(plot_data, cell_status == "expression_missing"),
      fill = "gray75",
      inherit.aes = FALSE,
      ggplot2::aes(x = barcode, y = method)
    ) +
    ggplot2::geom_point(
      data = legend_data,
      ggplot2::aes(x = x, y = y, color = status_label),
      inherit.aes = FALSE,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Barcode missing" = "white",
        "Expression missing/undefined" = "gray75"
      ),
      name = "Cell status"
    ) +
    ggplot2::labs(
      title = paste("Per-method expression heatmap for", gene_symbol),
      x = "Cell barcode",
      y = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(fill = "white", color = "gray85")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 2,
        override.aes = list(alpha = 1, shape = 22, size = 5, fill = c("white", "gray75"))
      ),
      fill = ggplot2::guide_colorbar(order = 1)
    )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .selectize-control.single .selectize-input,
      .selectize-dropdown .option {
        font-size: clamp(9px, 0.95vw, 12px);
      }
      .selectize-control.single .selectize-input > div.item {
        font-size: clamp(9px, 0.95vw, 12px);
      }
    "))
  ),
  titlePanel("scRNA-seq method comparison"),
  tabsetPanel(
    id = "analysis_tabs",
    tabPanel(
      "Cluster stability",
      sidebarLayout(
        sidebarPanel(
          width = 5,
          helpText("Select two methods (.rds) with matching bootstrap TSV files."),
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
    ),
    tabPanel(
      "Gene heatmap",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          helpText("Choose a gene to compare log-normalized expression across all methods and the union of all cell barcodes."),
          shiny::tagAppendAttributes(
            textInput("heatmap_gene", "Gene symbol", value = ""),
            list = "gene-symbol-options",
            autocomplete = "off"
          ),
          uiOutput("gene_symbol_datalist"),
          actionButton("refresh_heatmap", "Refresh gene list")
        ),
        mainPanel(
          width = 8,
          plotOutput("expression_heatmap", height = "420px"),
          verbatimTextOutput("heatmap_status")
        )
      )
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
      current_method1 <- isolate(input$method1)
      current_method2 <- isolate(input$method2)

      selected_method1 <- if (!is.null(current_method1) && current_method1 %in% pairs$rds) {
        current_method1
      } else {
        pairs$rds[[1]]
      }

      fallback_method2 <- if (nrow(pairs) >= 2) pairs$rds[[2]] else pairs$rds[[1]]
      selected_method2 <- if (!is.null(current_method2) && current_method2 %in% pairs$rds) {
        current_method2
      } else {
        fallback_method2
      }

      if (identical(selected_method1, selected_method2) && nrow(pairs) >= 2) {
        selected_method2 <- pairs$rds[[which(pairs$rds != selected_method1)[1]]]
      }

      updateSelectInput(session, "method1", choices = choice_map, selected = selected_method1)
      updateSelectInput(session, "method2", choices = choice_map, selected = selected_method2)
    } else {
      updateSelectInput(session, "method1", choices = c())
      updateSelectInput(session, "method2", choices = c())
    }
  }

  observeEvent(input$refresh, refresh_choices(), ignoreInit = TRUE)
  observeEvent(input$refresh_heatmap, {
    current_gene <- isolate(input$heatmap_gene)
    method_pairs(find_method_pairs("data"))

    if (!is.null(current_gene) && nzchar(current_gene)) {
      updateTextInput(session, "heatmap_gene", value = current_gene)
    }
  }, ignoreInit = TRUE)
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

  all_method_data <- reactive({
    pairs <- method_pairs()

    validate(
      need(nrow(pairs) > 0, "No valid method pairs found in ./data.")
    )

    purrr::map2(pairs$rds, pairs$boot, read_method_bundle)
  })

  output$gene_symbol_datalist <- renderUI({
    methods <- all_method_data()
    gene_choices <- methods %>%
      purrr::map("feature_table") %>%
      dplyr::bind_rows() %>%
      dplyr::distinct(gene_symbol) %>%
      dplyr::arrange(gene_symbol) %>%
      dplyr::pull(gene_symbol)

    tags$datalist(
      id = "gene-symbol-options",
      lapply(gene_choices, function(gene) {
        tags$option(value = gene)
      })
    )
  })

  heatmap_data <- reactive({
    methods <- all_method_data()
    req(input$heatmap_gene)
    validate(
      need(nzchar(input$heatmap_gene), "Type a gene symbol to draw the heatmap.")
    )

    make_expression_heatmap_data(methods, input$heatmap_gene)
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

  output$expression_heatmap <- renderPlot({
    plot_data <- heatmap_data()
    validate(
      need(nrow(plot_data) > 0, "No heatmap data available for the selected gene."),
      need(any(!is.na(plot_data$expression)), "Selected gene symbol was not found in the loaded methods.")
    )

    make_expression_heatmap_plot(plot_data, input$heatmap_gene)
  }, res = 110)

  output$status <- renderText({
    pairs <- method_pairs()

    if (nrow(pairs) == 0) {
      return("No valid method pairs found in ./data. Expected .rds plus matching bootstrap TSV file.")
    }

    paste0("Detected ", nrow(pairs), " method file pairs in ./data")
  })

  output$heatmap_status <- renderText({
    methods <- all_method_data()
    barcode_union <- sort(unique(unlist(purrr::map(methods, "barcodes"))))
    feature_union <- methods %>%
      purrr::map("feature_table") %>%
      dplyr::bind_rows() %>%
      dplyr::distinct(gene_symbol) %>%
      nrow()

    paste0(
      "Loaded ", length(methods), " methods. Heatmap rows are methods, columns are the union of ",
      length(barcode_union), " cell barcodes, and gene choices cover ",
      feature_union, " gene symbols."
    )
  })
}

shinyApp(ui = ui, server = server)
