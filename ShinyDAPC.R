# app.R

# Load required libraries
library(shiny)  # - interactive web applications with R
library(adegenet)  # tools for the exploration of genetic and genomic data.
library(ade4)  # Exploratory and Euclidean Methods in Environmental Sciences
library(ggplot2)  # - data visualisation
library(RColorBrewer)  # - colour palettes especially for thematic maps
library(shinybusy)  # - indicators to show the user that the server is busy
library(dplyr)  # tool for working with data frame like objects
library(ape)  # - exploration of phylogenetic data,
library(tibble)  # simple data frames

# compoplot dimensions helper
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}


# UI
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "cerulean"),
  titlePanel("DNA barcode exploration using adegenet (Jombart, T. et al.)"),
  add_busy_bar(timeout = 2500,
               col = "yellow",
               centered = TRUE,
               height = "8px"),
  tabsetPanel(
    # tab 1: designing the DAPC
    tabPanel("1.a. DAPC parameters",
             sidebarLayout(
               sidebarPanel(
                 fileInput("fasta9", "Upload FASTA file (DNA barcodes)", accept = c(".fa", ".fasta")),
                 helpText("Upload an aligned FASTA of equal length sequences, with no outgroup. Ensure location label is the word after a single underscore in each header (e.g. >id_Region)."),
                 hr(),
                 numericInput("maxnclust", code("maximum number of clusters to try"), value = 5, min = 1, step = 1),
                 actionButton("runClusters", strong("Step 1.  Find clusters")),
                 br(),
                 textOutput("clustersFound"),
                 br(),
                 downloadButton("downloadClustersCsv", "Download assigned clusters"),
                 hr(),
                 numericInput("npca", code("number of principal components to retain for DAPC"), value = 10, min = 1, step = 1),  
                 # After PCA reduces your genetic data, you must decide how many PCs to keep. These PCs are then used in the discriminant analysis step to maximize group separation. Retaining too few or too many PCs can distort your results.  
                 # Retain enough PCs to capture most of the genetic variation, but not so many that you overfit.  Too many PCs may separate individuals based on noise rather than true genetic structure.  Too few PCs may miss subtle differentiation.  
                 # So choose npca that explain a substantial proportion of variance (often 80–90%). 
                 numericInput("nda", code("number of discriminant functions to retain for DAPC"), value = 2, min = 1, step = 1),
                 #helpText("maximum = number of groups - 1 "), 
                 #br(), br(),
                 # The maximum number of discriminant functions is K – 1, where K = number of groups/clusters.  Keeping all possible discriminant functions (K – 1) maximizes discrimination, but may overfit noise. 
                 # Check the proportion of variance explained by each discriminant axis (shown in the scree plot). Often, retaining fewer (e.g., 2–3) is enough to capture the main structure.  If the first two axes explain most of the separation, you don’t need more.
                 actionButton("runDAPC", strong("Step 2.  Run DAPC")),
                 width = 3
               ),
               mainPanel(    
                 plotOutput("dapcScatter", height = "600px"),
                 verbatimTextOutput("dapcSummary")
               ))),
    
    # tab 2:  DAPC scatterplot
    tabPanel("1.b. DAPC Scatterplot",     
             sidebarLayout(
               sidebarPanel(
                 numericInput("xaxis", code("DAPC axis (x)"), value = 1, min = 1, step = 1),
                 numericInput("yaxis", code("DAPC axis (y)"), value = 2, min = 1, step = 1),
                 br(),
                 downloadButton("downloadDapcByLocationPlot", "Download DAPC scatterplot"),
                 width = 3),
               mainPanel(
                 verbatimTextOutput("dapcByLocInfo"),
                 plotOutput("dapcByLocation", height = "600px")
               ))),
    
    # Tab 3: Compoplot
    # shows how individuals are assigned to genetic clusters and the proportion of their membership in each cluster.  Different colours correspond to different genetic clusters (populations or groups).  
    # The length of each colored segment within a bar shows the probability (or proportion) that the individual belongs to that cluster.  Shows genetic assignment, not evolutionary history.  
    # Membership proportions are probabilities, not absolute truths.
    # - Solid colors = clear assignment.
    # - Mixed colors = admixture or uncertainty.
    # - Group-level patterns reveal how populations are structured genetically.
    tabPanel("2. Composition plot",
             sidebarLayout(  
               sidebarPanel(
                 fileInput("fasta11", "Upload FASTA file (DNA barcodes)", accept = c(".fa", ".fasta")),
                 helpText("Upload an aligned FASTA of equal length sequences, with no outgroup. Ensure location label is the word after a single underscore in each header (e.g. >id_Region)."),
                 hr(),
                 numericInput("nPcDAPC", code("number of principal components to retain for compoplot"), value = 10, min = 1, step = 1),
                 numericInput("nDa2", code("number of discriminant functions to retain for compoplot"), value = 2, min = 1, step = 1),
                 actionButton("go", strong("Run DAPC")),
                 hr(),
                 selectInput("palette2", "Colour palette:",
                             choices = c("Set1", "Set2", "Set3", "Paired", "Dark2", "Accent"),
                             selected = "Set3"),
                 numericInput("plot_w", "Download width (pixels)", value = 1500, min = 200, step = 50),
                 numericInput("plot_h", "Download height (pixels)", value = 600, min = 200, step = 50),
                 width = 3
               ),
               mainPanel(    
                 plotOutput("compoplot", height = "600px"),
                 downloadButton("downloadPlot", "Download compoplot (png)"),
                 br(), br(), br(),
                 plotOutput("dapcScatter2", height = "600px"),
                 br(), br(),
                 verbatimTextOutput("dapcSummary")
               )))
  )
)


# Server
server <- function(input, output, session) {
  
  # 1.a. design the DAPC
  fasta_raw <- reactive({
    req(input$fasta9)
    tryCatch(read.dna(file = input$fasta9$datapath, format = "fasta"),
             error = function(e) { showNotification("Error reading fasta", type = "error"); NULL })
  })
  
  output$fasta9Info <- renderPrint({
    dna4 <- fasta_raw(); req(dna4)
    headers <- rownames(dna4)
    cat("Number of sequences:", length(headers), "\n")
    cat("Example headers (first 10):\n"); print(head(headers, 10))
    seq_mat <- as.character(dna4)
    seq_lengths <- apply(seq_mat, 1, function(x) sum(nchar(x) > 0))
    cat("\nSequence lengths summary:\n"); print(summary(seq_lengths))
  })
  
  samples_and_groups <- reactive({
    dna4 <- fasta_raw(); req(dna4)
    headers <- rownames(dna4)
    group_from_header <- function(h) {
      parts <- unlist(strsplit(h, "_"))
      if (length(parts) >= 2) tail(parts, 1) else NA_character_
    }
    groups <- sapply(headers, group_from_header, USE.NAMES = FALSE)
    data.frame(sample = headers, group = groups, stringsAsFactors = FALSE)
  })
  
  genind_obj <- reactive({
    dna4 <- fasta_raw(); req(dna4)
    gi <- tryCatch(DNAbin2genind(dna4), error = function(e) { showNotification("DNAbin2genind error", type = "error"); NULL })
    sampgrp <- samples_and_groups()
    if (!is.null(gi) && !all(is.na(sampgrp$group))) {
      ind_order <- indNames(gi)
      grp_vec <- sampgrp$group[match(ind_order, sampgrp$sample)]
      if (!all(is.na(grp_vec))) pop(gi) <- as.factor(grp_vec)
    }
    gi
  })
  
  clusters_res <- eventReactive(input$runClusters, {
    gi <- genind_obj(); req(gi)
    maxn <- as.integer(input$maxnclust); validate(need(!is.na(maxn) && maxn >= 1, "Set max.n.clust >=1"))
    npca_auto <- min(50, max(1, nInd(gi) - 1), nLoc(gi))
    fc_raw <- tryCatch(find.clusters(gi, max.n.clust = maxn, n.pca = npca_auto, choose.n.clust = FALSE),
                       error = function(e) { showNotification(paste("find.clusters failed:", e$message), type = "error"); NULL })
    # find.clustersreduces dimensionality with Principal Component Analysis (PCA), then applies k-means clustering on the retained principal components.  
    # The function computes the Bayesian Information Criterion (BIC) for different values of k. 
    # The optimal number of clusters is usually the k where BIC reaches its minimum or levels off.  
    # These clusters can then be used as the grouping factor in DAPC.
    fc_raw_str <- if (!is.null(fc_raw)) paste(capture.output(str(fc_raw, max.level = 2)), collapse = "\n") else "NULL"
    fc <- NULL
    if (!is.null(fc_raw) && is.list(fc_raw) && !is.null(fc_raw$grp)) {
      fc <- fc_raw
    } else if (!is.null(fc_raw) && is.atomic(fc_raw) && length(fc_raw) == nInd(gi)) {
      fc <- list(grp = as.integer(fc_raw), Kstat = NULL)
      showNotification("find.clusters returned atomic vector; using as cluster assignments.", type = "warning")
    } else {
      showNotification("find.clusters returned unexpected object. See diagnostics.", type = "error")
      message("find.clusters raw output:"); message(fc_raw_str)
    }
    clusters_table <- NULL
    if (!is.null(fc) && !is.null(fc$grp) && length(fc$grp) == nInd(gi)) {
      clusters_table <- data.frame(sample = indNames(gi), cluster = as.integer(fc$grp), stringsAsFactors = FALSE)
    }
    list(fc_raw = fc_raw, fc_raw_str = fc_raw_str, fc = fc, npca_used = npca_auto, clusters_table = clusters_table)
  })
  
  
  # server: show number of clusters found after find.clusters finishes
  output$clustersFound <- renderText({
    res <- clusters_res()
    if (is.null(res) || is.null(res$fc)) {
      return("Clusters: not run")
    }
    fc <- res$fc
    
    # Prefer explicit cluster size vector if present
    if (!is.null(fc$size) && length(fc$size) > 0) {
      nclust <- length(fc$size)
      sizes_text <- paste0(" (sizes: ", paste(fc$size, collapse = ", "), ")")
      return(paste0("Clusters: ", nclust, sizes_text))
    }
    
    # Fallback: derive from fc$grp if present
    if (!is.null(fc$grp)) {
      nclust <- length(unique(as.integer(fc$grp)))
      return(paste0("Clusters: ", nclust))
    }
    
    # Fallback: derive from Kstat names if present
    if (!is.null(fc$Kstat) && !is.null(names(fc$Kstat))) {
      # choose the K with minimum BIC as the "found" cluster (optional heuristic)
      Kvals <- suppressWarnings(as.integer(gsub("K=", "", names(fc$Kstat))))
      if (!any(is.na(Kvals))) {
        bestK <- Kvals[which.min(as.numeric(fc$Kstat))]
        return(paste0("Clusters (best by BIC): ", bestK))
      }
    }
    
    "Clusters: result available but number not determined"
  })
  
  observeEvent(clusters_res(), {
    res <- clusters_res()
    if (is.null(res) || is.null(res$fc)) {
      showNotification("find.clusters finished with no result", type = "warning")
    } else {
      # attempt to get number
      ncl <- if (!is.null(res$fc$size)) length(res$fc$size) else if (!is.null(res$fc$grp)) length(unique(as.integer(res$fc$grp))) else NA
      if (!is.na(ncl)) showNotification(paste0("find.clusters finished: ", ncl, " clusters found"), type = "message")
    }
  })
  
  output$clustersTable <- renderTable({
    res <- clusters_res(); req(res)
    res$clusters_table
  }, striped = TRUE, hover = TRUE)
  
  output$downloadClustersCsv <- downloadHandler(
    filename = function() {
      paste0("assigned_clusters_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # Try normalized find.clusters table first
      cres <- tryCatch(clusters_res(), error = function(e) NULL)
      out_df <- NULL
      
      if (!is.null(cres) && !is.null(cres$clusters_table)) {
        out_df <- cres$clusters_table
      }
      
      # Fallback: if no clusters_table, try DAPC posterior assignments
      if (is.null(out_df)) {
        dres <- tryCatch(dapc_res(), error = function(e) NULL)
        if (!is.null(dres) && !is.null(dres$dapc) && !is.null(dres$dapc$assign)) {
          assign_mat <- dres$dapc$assign
          if (is.null(dim(assign_mat)) && length(assign_mat) > 0) assign_mat <- matrix(assign_mat, ncol = 1)
          if (!is.null(dim(assign_mat)) && nrow(assign_mat) > 0) {
            sample_names <- rownames(assign_mat)
            gi <- tryCatch(genind_obj(), error = function(e) NULL)
            # fallback to genind / fasta order if rownames missing
            if (is.null(sample_names) || length(sample_names) != nrow(assign_mat)) {
              sg <- tryCatch(samples_and_groups(), error = function(e) NULL)
              if (!is.null(sg) && nrow(sg) == nrow(assign_mat)) sample_names <- as.character(sg$sample)
              else if (!is.null(gi) && inherits(gi, "genind") && length(indNames(gi)) == nrow(assign_mat)) sample_names <- indNames(gi)
              else sample_names <- paste0("ind", seq_len(nrow(assign_mat)))
            }
            assigned <- apply(assign_mat, 1, function(r) {
              if (all(is.na(r))) return(NA_integer_)
              which.max(as.numeric(r))
            })
            out_df <- data.frame(sample = sample_names, cluster = as.integer(assigned), stringsAsFactors = FALSE)
          }
        }
      }
      
      # Fallback: use parsed FASTA header groups (samples_and_groups)
      if (is.null(out_df)) {
        sg <- tryCatch(samples_and_groups(), error = function(e) NULL)
        if (!is.null(sg) && nrow(sg) > 0) {
          out_df <- data.frame(sample = as.character(sg$sample), cluster = as.character(sg$group), stringsAsFactors = FALSE)
        }
      }
      
      # If still nothing, write a small explanatory CSV
      if (is.null(out_df)) {
        out_df <- data.frame(message = "No cluster information available. Run find.clusters or DAPC first.")
        write.csv(out_df, file, row.names = FALSE)
        return()
      }
      
      # Coerce types sensibly
      out_df$sample <- as.character(out_df$sample)
      if (is.factor(out_df$cluster)) out_df$cluster <- as.character(out_df$cluster)
      cluster_int <- suppressWarnings(as.integer(out_df$cluster))
      if (!all(is.na(cluster_int)) && sum(!is.na(cluster_int)) / nrow(out_df) >= 0.5) {
        out_df$cluster <- cluster_int
      } else {
        out_df$cluster <- as.character(out_df$cluster)
      }
      
      # If possible order rows to match FASTA header order
      sg_all <- tryCatch(samples_and_groups(), error = function(e) NULL)
      if (!is.null(sg_all) && all(sg_all$sample %in% out_df$sample)) {
        out_df <- out_df[match(sg_all$sample, out_df$sample), , drop = FALSE]
      }
      
      write.csv(out_df, file, row.names = FALSE)
    }
  )
  
  
  # 1.b. DAPC Scatterplot
  dapc_res <- eventReactive(input$runDAPC, {
    gi <- genind_obj(); req(gi)
    fc <- NULL; if (!is.null(clusters_res())) fc <- clusters_res()$fc
    if (!is.null(fc)) grp <- fc$grp else if (!is.null(pop(gi))) grp <- pop(gi) else { showNotification("No groups available", type = "error"); return(NULL) }
    npca <- as.integer(input$npca); nda <- as.integer(input$nda)
    validate(need(!is.null(npca) && npca >= 1, "n.pca >=1"), need(!is.null(nda) && nda >= 1, "n.da >=1"))
    nind <- nInd(gi); ngroups3 <- length(unique(as.character(grp)))
    max_npca_allowed <- max(1, nind - ngroups3)
    if (npca > max_npca_allowed) { npca <- max_npca_allowed; showNotification(paste("n.pca reduced to", npca), type = "message") }
    max_nda_allowed <- min(ngroups3 - 1, npca)
    if (nda > max_nda_allowed) { nda <- max_nda_allowed; showNotification(paste("n.da reduced to", nda), type = "message") }
    dapc_obj <- tryCatch(dapc(gi, pop = grp, n.pca = npca, n.da = nda), error = function(e) { showNotification(paste("dapc error:", e$message), type = "error"); NULL })
    list(dapc = dapc_obj, grp = grp, npca = npca, nda = nda)
  })
  
  output$dapcSummary <- renderPrint({
    res <- dapc_res(); req(res)
    if (is.null(res$dapc)) { cat("DAPC failed or not run.\n"); return() }
    print(res$dapc)
  })
  
  observe({
    res <- dapc_res()
    if (is.null(res) || is.null(res$nda)) return()
    nda <- as.integer(res$nda)
    updateNumericInput(session, "xaxis", value = min(input$xaxis, nda), min = 1, max = nda, step = 1)
    updateNumericInput(session, "yaxis", value = ifelse(input$yaxis <= nda, input$yaxis, ifelse(nda>=2, 2, 1)), min = 1, max = nda, step = 1)
  })
  
  # Restored DAPC scatter: uses adegenet::scatter for default plotting
  output$dapcScatter <- renderPlot({
    res <- dapc_res(); req(res)
    dapc_obj <- res$dapc; req(dapc_obj)
    # Use adegenet's scatter function which shows individuals and groups
    # This will use the DAs retained; we keep defaults for appearance
    scatter(dapc_obj, posi.da = "bottomright", scree.pca = TRUE, cex = 2.5, clab = 0, pch=20, ratio.da=0.2, ratio.pca=0.2, leg=TRUE,)  # , cstar=0
  })
  
  output$dapcAssign <- renderTable({
    res <- dapc_res(); if (is.null(res) || is.null(res$dapc)) return(data.frame(message = "No DAPC result"))
    dapc_obj <- res$dapc
    assign_mat <- dapc_obj$assign
    if (is.null(assign_mat)) return(data.frame(message = "DAPC has no 'assign' matrix"))
    if (is.null(dim(assign_mat))) {
      if (length(assign_mat) == 0) return(data.frame(message = "DAPC 'assign' is empty"))
      assign_mat <- matrix(assign_mat, ncol = 1); colnames(assign_mat) <- "score"
    }
    if (nrow(assign_mat) < 1) return(data.frame(message = "DAPC 'assign' has no rows"))
    sample_names <- rownames(assign_mat)
    if (is.null(sample_names) || any(sample_names == "")) {
      gi <- tryCatch(genind_obj(), error = function(e) NULL)
      if (!is.null(gi) && inherits(gi, "genind") && length(indNames(gi)) == nrow(assign_mat)) sample_names <- indNames(gi) else sample_names <- paste0("ind", seq_len(nrow(assign_mat)))
    }
    assigned <- apply(assign_mat, 1, function(r) if (all(is.na(r)) ) NA_integer_ else which.max(as.numeric(r)))
    data.frame(sample = sample_names, assigned = as.integer(assigned), stringsAsFactors = FALSE)
  })
  
  output$downloadDapcPlot <- downloadHandler(
    filename = function() {
      paste0("dapc_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      res <- dapc_res()
      if (is.null(res) || is.null(res$dapc)) {
        # create a simple PNG telling the user no DAPC available
        png(file, width = 1000, height = 800)
        plot.new()
        title("No DAPC result available. Run DAPC first.")
        dev.off()
        return()
      }
      
      dapc_obj <- res$dapc
      
      # Determine width/height/resolution
      png(file, width = 1200, height = 900, res = 150)
      
      # attempt to extract chosen axes and draw a scatter of those axes if possible
      use_axes <- FALSE
      xax <- NULL; yax <- NULL
      if (!is.null(input$xaxis) && !is.null(input$yaxis)) {
        xax <- as.integer(input$xaxis); yax <- as.integer(input$yaxis)
        coords <- tryCatch(as.data.frame(dapc_obj$ind.coord), error = function(e) NULL)
        if (!is.null(coords) && ncol(coords) >= max(xax, yax)) use_axes <- TRUE
      }
      
      if (use_axes) {
        # Draw a simple base R scatter of the selected DAPC axes to match app choices
        xvals <- dapc_obj$ind.coord[, xax]
        yvals <- dapc_obj$ind.coord[, yax]
        samp_names <- rownames(dapc_obj$ind.coord)
        if (is.null(samp_names)) samp_names <- paste0("ind", seq_len(length(xvals)))
        
        # Colouring: try to use pop info if present, otherwise draw plain points
        grp <- NULL
        if (!is.null(res$grp)) {
          grp <- res$grp
          if (!is.null(names(grp))) grp <- grp[match(samp_names, names(grp))]
        } else {
          gi <- tryCatch(genind_obj(), error = function(e) NULL)
          if (!is.null(gi) && inherits(gi, "genind")) grp <- pop(gi)
        }
        if (is.null(grp)) grp <- factor(rep("unknown", length(xvals))) else grp <- factor(grp)
        
        palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, length(levels(grp)))), "Set2"))(length(levels(grp)))
        cols3 <- palette[as.integer(grp)]
        
        par(mar = c(5,5,4,8))
        plot(xvals, yvals, pch = 21, bg = cols3, col = "black", cex = 2,
             xlab = paste0("DAPC axis ", xax), ylab = paste0("DAPC axis ", yax),
             main = paste0("DAPC scatter: axis ", xax, " vs ", yax))
        # add legend on right
        legend("topright", inset = c(-0.25,0), legend = levels(grp), pt.bg = palette, pch = 21, title = "Group", xpd = TRUE)
      } else {
        # fallback: use adegenet::scatter to reproduce the original app plot
        tryCatch({
          scatter(dapc_obj, posi.da = "bottomright", scree.pca = TRUE, cex = 2, clab = 0.8)
        }, error = function(e) {
          plot.new()
          title("Error drawing DAPC scatter: see console for details")
        })
      }
      
      dev.off()
    }
  )
  
  output$dapcByLocInfo <- renderPrint({
    res <- dapc_res()
    if (is.null(res) || is.null(res$dapc)) {
      cat("No DAPC available. Run DAPC first.")
      return()
    }
    #cat("DAPC scatter coloured by FASTA locations (substring after underscore).")
  })
  
  output$dapcByLocation <- renderPlot({
    res <- dapc_res()
    req(res)
    dapc_obj <- res$dapc
    req(dapc_obj)
    
    coords <- tryCatch(as.data.frame(dapc_obj$ind.coord), error = function(e) NULL)
    if (is.null(coords) || ncol(coords) < 2) {
      plot.new(); title("DAPC coordinates not available or fewer than 2 axes."); return()
    }
    
    # use axis selectors if present, otherwise default to first two columns
    xax <- if (!is.null(input$xaxis)) as.integer(input$xaxis) else 1
    yax <- if (!is.null(input$yaxis)) as.integer(input$yaxis) else 2
    if (xax < 1) xax <- 1
    if (yax < 1) yax <- 2
    if (ncol(coords) < max(xax, yax)) {
      plot.new(); title("Selected axes not available"); return()
    }
    
    plot_df <- data.frame(Axis1 = coords[[xax]], Axis2 = coords[[yax]], stringsAsFactors = FALSE)
    plot_df$sample <- rownames(coords)
    if (is.null(plot_df$sample)) plot_df$sample <- paste0("ind", seq_len(nrow(plot_df)))
    
    # get FASTA-parsed locations
    sg <- tryCatch(samples_and_groups(), error = function(e) NULL)
    if (!is.null(sg) && nrow(sg) > 0) {
      matched <- sg$group[match(plot_df$sample, sg$sample)]
      matched[is.na(matched)] <- "unknown"
      plot_df$location <- factor(matched)
    } else {
      plot_df$location <- factor(rep("unknown", nrow(plot_df)))
    }
    
    # palette
    nloc <- length(levels(plot_df$location))
    if (nloc <= 8) {
      pal <- RColorBrewer::brewer.pal(max(3, nloc), "Set2")[seq_len(nloc)]
    } else {
      pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nloc)
    }
    
    # ggplot filled circles shape 21 with legend
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Axis1, y = Axis2, fill = location)) +
      ggplot2::geom_point(shape = 21, color = "black", size = 4, alpha = 0.95) +
      ggplot2::scale_fill_manual(values = pal, name = "Location") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(x = paste0("DAPC axis ", xax), y = paste0("DAPC axis ", yax)) +
      ggplot2::theme(legend.position = "right",
                     legend.title = ggplot2::element_text(size = 11),
                     legend.text = ggplot2::element_text(size = 10))
    
    print(p)
  })
  
  
  output$downloadDapcByLocationPlot <- downloadHandler(
    filename = function() {
      paste0("dapc_by_location_", Sys.Date(), ".png")
    },
    content = function(file) {
      res <- tryCatch(dapc_res(), error = function(e) NULL)
      if (is.null(res) || is.null(res$dapc)) {
        png(file, width = 1200, height = 900, res = 150)
        plot.new()
        title("No DAPC result available. Run DAPC first.")
        dev.off()
        return()
      }
      
      dapc_obj <- res$dapc
      coords <- tryCatch(as.data.frame(dapc_obj$ind.coord), error = function(e) NULL)
      if (is.null(coords) || ncol(coords) < 2) {
        png(file, width = 1200, height = 900, res = 150)
        plot.new()
        title("DAPC coordinates not available or fewer than 2 axes.")
        dev.off()
        return()
      }
      
      # determine axes to plot
      xax <- if (!is.null(input$xaxis)) as.integer(input$xaxis) else 1
      yax <- if (!is.null(input$yaxis)) as.integer(input$yaxis) else 2
      if (xax < 1) xax <- 1
      if (yax < 1) yax <- 2
      if (ncol(coords) < max(xax, yax)) {
        png(file, width = 1200, height = 900, res = 150)
        plot.new()
        title("Selected axes not available")
        dev.off()
        return()
      }
      
      plot_df <- data.frame(Axis1 = coords[[xax]], Axis2 = coords[[yax]], stringsAsFactors = FALSE)
      plot_df$sample <- rownames(coords)
      if (is.null(plot_df$sample)) plot_df$sample <- paste0("ind", seq_len(nrow(plot_df)))
      
      sg <- tryCatch(samples_and_groups(), error = function(e) NULL)
      if (!is.null(sg) && nrow(sg) > 0) {
        matched <- sg$group[match(plot_df$sample, sg$sample)]
        matched[is.na(matched)] <- "unknown"
        plot_df$location <- factor(matched)
      } else {
        plot_df$location <- factor(rep("unknown", nrow(plot_df)))
      }
      
      nloc <- length(levels(plot_df$location))
      if (nloc <= 8) {
        pal <- RColorBrewer::brewer.pal(max(3, nloc), "Set2")[seq_len(nloc)]
      } else {
        pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nloc)
      }
      
      # create PNG and draw ggplot 
      png(file, width = 1200, height = 900, res = 150)
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Axis1, y = Axis2, fill = location)) +
        ggplot2::geom_point(shape = 21, color = "black", size = 3, alpha = 0.95) +
        ggplot2::scale_fill_manual(values = pal, name = "Location") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(x = paste0("DAPC axis ", xax), y = paste0("DAPC axis ", yax)) +
        ggplot2::theme(axis.title = ggplot2::element_text(size = 10),   # << smaller axis label text
                       legend.position = "right",
                       legend.title = ggplot2::element_text(size = 11),
                       legend.text = ggplot2::element_text(size = 10))
      print(p)
      dev.off()
    }
  ) 
  
  # 2. Compoplot
  seqs_dnabin <- reactive({
    req(input$fasta11)
    inFile <- input$fasta11$datapath
    dna7 <- tryCatch({
      read.dna(inFile, format = "fasta", as.character = FALSE)
    }, error = function(e2){
      stop("Failed to read FASTA: ensure file is a proper DNA FASTA.")
    })
    dna7
  })
  
  genind_obj2 <- reactive({
    dna7 <- seqs_dnabin()
    
    gi7 <- tryCatch({
      DNAbin2genind(dna7)
    }, error = function(e2) {
      # Fallback: collapse each full sequence into a single categorical locus
      dna_char <- as.character(dna7)
      seqs7 <- apply(dna_char, 1, paste0, collapse = "")
      df7 <- data.frame(seq = as.factor(seqs7), stringsAsFactors = TRUE)
      gi7 <- df2genind(df7, ploidy = 1, ncode = 1)
      indNames(gi7) <- rownames(dna_char)
      return(gi7)
    })
    
    # Ensure individual names
    dna_char <- as.character(dna7)
    rn <- rownames(dna_char)
    if (!is.null(rn)) indNames(gi7) <- rn
    
    gi7
  })
  
  # Extract groups from header names: the word immediately following the first underscore
  fasta_groups <- reactive({
    dna7 <- seqs_dnabin()
    rn <- rownames(as.character(dna7))
    if (is.null(rn)) {
      stop("FASTA headers not found; ensure headers are present")
    }
    # For each header, split at first underscore and take the token after it.
    groups7 <- sapply(rn, function(x7) {
      parts7 <- unlist(strsplit(x7, "_"))
      if (length(parts7) >= 2) {
        # take the first token after the first underscore
        token <- parts7[2]
        # if token contains whitespace or additional separators, take the first word
        token2 <- unlist(strsplit(token, "[[:space:]:;,|]"))[1]
        token2
      } else {
        "unknown"
      }
    }, USE.NAMES = FALSE)
    groups7
  })
  
  # Run DAPC using original groups when Go pressed
  dapc_result <- eventReactive(input$go, {
    gi7 <- genind_obj2()
    groups7 <- fasta_groups()
    
    if (length(groups7) != nInd(gi7)) stop("Number of groups extracted does not match number of sequences")
    
    pop(gi7) <- as.factor(groups7)
    
    nind7 <- nInd(gi7)
    npop7 <- nPop(gi7)
    # safe upper limit for n.pca
    max_possible_pcs <- max(1, min(nind7 - npop7, nLoc(gi7)))
    n_pca_dapc <- min(as.integer(input$nPcDAPC), max_possible_pcs)
    
    n_da <- min(as.integer(input$nDa2), npop7 - 1)
    if (n_da < 1) n_da <- 1
    
    dapc_obj7 <- dapc(gi7, pop = groups7, n.pca = n_pca_dapc, n.da = n_da)
    
    list(dapc = dapc_obj7, groups = groups7, gi = gi7)
  })
  
  output$compoplot <- renderPlot({
    req(dapc_result())
    res7 <- dapc_result()
    dapc_obj7 <- res7$dapc
    if (is.null(dapc_obj7) || is.null(dapc_obj7$posterior)) stop("DAPC result or posterior matrix missing")
    
    post <- dapc_obj7$posterior
    ind_names <- rownames(post)
    if (is.null(ind_names)) {
      indn <- indNames(genind_obj2())
      if (!is.null(indn)) rownames(post) <- indn
      ind_names <- rownames(post)
    }
    if (is.null(ind_names)) stop("Cannot determine individual names for ordering")
    
    prior_group <- sapply(ind_names, function(x) {
      parts7 <- strsplit(x, "_", fixed = TRUE)[[1]]
      if (length(parts7) >= 2) {
        token <- parts7[2]
        strsplit(token, "[[:space:]:;,|]")[[1]][1] # \\s:;,|
      } else {
        "unknown"
      }
    }, USE.NAMES = FALSE)
    
    ord <- order(prior_group, ind_names)
    post_ord <- post[ord, , drop = FALSE]
    ind_names_ord <- ind_names[ord]
    prior_group_ord <- prior_group[ord]
    
    K7 <- ncol(post_ord)
    basePal7 <- RColorBrewer::brewer.pal(max(3, min(12, K7)), input$palette2)
    if (K7 > length(basePal7)) basePal7 <- colorRampPalette(basePal7)(K7)
    cols7 <- setNames(basePal7[seq_len(K7)], levels(prior_group))
    
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    
    # increase right margin for legend
    par(mar = c(10, 4, 4, 12) + 0.1, las = 2, xpd = NA)
    
    n_bars <- nrow(post_ord)
    names_arg <- rep("", n_bars)
    
    bp <- barplot(t(post_ord),
                  col = cols7,
                  border = NA,
                  space = 0,
                  axes = FALSE,
                  ylim = c(0, 1),
                  names.arg = names_arg)
    
    axis(2)
    text(x = bp, y = par("usr")[3] - 0.02 * (par("usr")[4] - par("usr")[3]),
         labels = ind_names_ord, srt = 90, adj = c(1, 0.5), cex = 0.7, xpd = NA)
    
    group_ends <- tapply(seq_along(prior_group_ord), prior_group_ord, max)
    group_end_idxs <- as.integer(group_ends)
    if (length(group_end_idxs) > 1 && length(bp) >= 2) {
      for (i in seq_len(length(group_end_idxs)-1)) {
        vpos <- (bp[group_end_idxs[i]] + bp[group_end_idxs[i] + 1]) / 2
        abline(v = vpos, col = "white", lwd = 1)
      }
    }
    
    box()
    
  },
  width = reactive({ max(200, input$plot_w %||% 1400) }),
  height = reactive({ max(200, input$plot_h %||% 1000) })
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("compoplot_dapc_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(dapc_result())
      # validate and coerce width/height
      w <- as.integer(max(200, input$plot_w %||% 2000))
      h <- as.integer(max(200, input$plot_h %||% 600))
      
      res7 <- dapc_result()
      dapc_obj7 <- res7$dapc
      if (is.null(dapc_obj7) || is.null(dapc_obj7$posterior)) stop("DAPC result or posterior matrix missing")
      
      post <- dapc_obj7$posterior
      ind_names <- rownames(post)
      if (is.null(ind_names)) {
        indn <- indNames(genind_obj2())
        if (!is.null(indn)) rownames(post) <- indn
        ind_names <- rownames(post)
      }
      if (is.null(ind_names)) stop("Cannot determine individual names for ordering")
      
      # extract prior groups and order as before
      prior_group <- sapply(ind_names, function(x8) {
        parts7 <- strsplit(x8, "_", fixed = TRUE)[[1]]
        if (length(parts7) >= 2) {
          token <- parts7[2]
          strsplit(token, "[[:space:]:;,|]")[[1]][1]  #  \\s:;,|
        } else {
          "unknown"
        }
      }, USE.NAMES = FALSE)
      ord <- order(prior_group, ind_names)
      post_ord <- post[ord, , drop = FALSE]
      ind_names_ord <- ind_names[ord]
      prior_group_ord <- prior_group[ord]
      
      K7 <- ncol(post_ord)
      basePal7 <- RColorBrewer::brewer.pal(max(3, min(12, K7)), input$palette2)
      if (K7 > length(basePal7)) basePal <- colorRampPalette(basePal7)(K7)
      cols7 <- setNames(basePal7[seq_len(K7)], levels(prior_group))
      
      png(filename = file, width = w, height = h, res = 150)
      oldpar <- par(no.readonly = TRUE)
      on.exit({ par(oldpar); dev.off() }, add = TRUE)
      
      par(mar = c(10, 4, 4, 12) + 0.1, las = 2, xpd = NA)
      
      n_bars <- nrow(post_ord)
      names_arg <- rep("", n_bars)
      
      bp <- barplot(t(post_ord),
                    col = cols7,
                    border = NA,
                    space = 0,
                    axes = FALSE,
                    ylim = c(0, 1),
                    names.arg = names_arg)
      
      axis(2)
      text(x = bp, y = par("usr")[3] - 0.02 * (par("usr")[4] - par("usr")[3]),
           labels = ind_names_ord, srt = 90, adj = c(1, 0.5), cex = 0.7, xpd = NA)
      
      group_ends <- tapply(seq_along(prior_group_ord), prior_group_ord, max)
      group_end_idxs <- as.integer(group_ends)
      if (length(group_end_idxs) > 1 && length(bp) >= 2) {
        for (i in seq_len(length(group_end_idxs)-1)) {
          vpos <- (bp[group_end_idxs[i]] + bp[group_end_idxs[i] + 1]) / 2
          abline(v = vpos, col = "white", lwd = 1)
        }
      }
      
      box()
    }
  )
  
  
  #2.a. design the DAPC
  output$dapcSummary2 <- renderPrint({
    res7 <- dapc_result(); req(res7)
    if (is.null(res7$dapc)) { cat("DAPC failed or not run.\n"); return() }
    print(res7$dapc)
  })
  
  observe({
    res7 <- dapc_result()
    if (is.null(res7) || is.null(res7$n.Da)) return()
    n.Da <- as.integer(res7$n.Da)
    updateNumericInput(session, "xaxis", value = min(input$xaxis, n.Da), min = 1, max = n.Da, step = 1)
    updateNumericInput(session, "yaxis", value = ifelse(input$yaxis <= n.Da, input$yaxis, ifelse(n.Da>=2, 2, 1)), min = 1, max = n.Da, step = 1)
  })
  
  # Restored DAPC scatter. uses adegenet::scatter for default plotting
  output$dapcScatter2 <- renderPlot({
    res7 <- dapc_result(); req(res7)
    dapc_obj7 <- res7$dapc; req(dapc_obj7)
    # Use adegenet's scatter function which shows individuals and groups
    # This will use the DAs retained; we keep defaults for appearance
    scatter(dapc_obj7, posi.da = "bottomright", scree.pca = TRUE, cex = 2.5, clab = 0, pch=20, ratio.da=0.3, ratio.pca=0.3, leg=TRUE, cellipse=0, cstar=0)  # , leg=FALSE   # , label = levels(grp)
    #legend("topright", legend = levels(grp), cex = 0.8, text.width = max(strwidth(levels(grp))) * 1.2, xpd = TRUE)
  }) 
  
  
}

shinyApp(ui, server)
