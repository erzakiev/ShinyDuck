nextflowModuleUI <- function(id) {
  ns <- NS(id)

  tagList(
    fileInput(ns("files"), "Select FASTQ files", multiple = TRUE),

    actionButton(ns("run"), "Run Nextflow"),
    actionButton(ns("stop"), "Stop"),

    tags$hr(),

    textOutput(ns("status")),
    verbatimTextOutput(ns("log"), placeholder = TRUE)
  )
}

nextflowModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    rv <- reactiveValues(
      process = NULL,
      log = "",
      status = "idle",
      progress = NULL,
      output_dir = NULL
    )

    # ---------------------------
    # 🧠 helper: append log safely
    # ---------------------------
    append_log <- function(old, new, max_lines = 500) {
      combined <- c(strsplit(old, "\n")[[1]], new)
      paste(tail(combined, max_lines), collapse = "\n")
    }

    # ---------------------------
    # 🧠 helper: create sample sheet
    # ---------------------------
    create_samplesheet <- function(files) {
      df <- data.frame(
        sample = tools::file_path_sans_ext(basename(files)),
        fastq = files
      )

      tmp <- tempfile(fileext = ".csv")
      write.csv(df, tmp, row.names = FALSE)
      tmp
    }

    # ---------------------------
    # 🧠 helper: parse progress
    # ---------------------------
    parse_progress <- function(lines) {
      pattern <- "\\[(\\d+)/(\\d+)\\]"

      for (line in rev(lines)) {
        m <- regmatches(line, regexec(pattern, line))[[1]]
        if (length(m) == 3) {
          return(as.numeric(m[2]) / as.numeric(m[3]))
        }
      }
      NULL
    }

    # ---------------------------
    # ▶️ RUN
    # ---------------------------
    observeEvent(input$run, {
      req(input$files)

      # 🚫 already running
      if (!is.null(rv$process) && rv$process$is_alive()) {
        showNotification("Process already running", type = "error")
        return()
      }

      samplesheet <- create_samplesheet(input$files$datapath)
      outdir <- file.path(tempdir(), paste0("nf_", Sys.time()))
      dir.create(outdir, recursive = TRUE)

      cmd <- "nextflow"
      args <- c(
        "run", "your_pipeline.nf",
        "--input", samplesheet,
        "--outdir", outdir,
        "-ansi-log", "false"
      )

      rv$process <- processx::process$new(
        command = cmd,
        args = args,
        stdout = "|",
        stderr = "|"
      )

      rv$log <- "[Starting Nextflow]\n"
      rv$status <- "running"
      rv$output_dir <- outdir
    })

    # ---------------------------
    # ⏹ STOP
    # ---------------------------
    observeEvent(input$stop, {
      if (!is.null(rv$process) && rv$process$is_alive()) {
        rv$process$kill()
        rv$status <- "killed"
        rv$log <- paste(rv$log, "\n[Process killed]")
      }
    })

    # ---------------------------
    # 🔄 LOG POLLING
    # ---------------------------
    observe({
      invalidateLater(500)

      req(rv$process)

      if (rv$process$is_alive()) {

        out <- rv$process$read_output_lines()
        err <- rv$process$read_error_lines()
        new_lines <- c(out, err)

        if (length(new_lines)) {
          rv$log <- append_log(rv$log, new_lines)

          # 📊 update progress
          p <- parse_progress(new_lines)
          if (!is.null(p)) {
            rv$progress <- p
          }
        }

      } else if (rv$status == "running") {
        # финал
        out <- rv$process$read_all_output_lines()
        err <- rv$process$read_all_error_lines()

        rv$log <- append_log(rv$log, c(out, err, "[Process finished]"))
        rv$status <- "finished"
      }
    })

    # ---------------------------
    # 📤 OUTPUTS
    # ---------------------------
    output$log <- renderText({
      rv$log
    })

    output$status <- renderText({
      if (rv$status == "running" && !is.null(rv$progress)) {
        paste0("Running: ", round(rv$progress * 100), "%")
      } else {
        paste("Status:", rv$status)
      }
    })

    # ---------------------------
    # 🔁 RETURN VALUE (important!)
    # ---------------------------
    return(
      reactive({
        list(
          status = rv$status,
          progress = rv$progress,
          output_dir = rv$output_dir
        )
      })
    )

  })
}
