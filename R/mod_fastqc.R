mod_fastqc_ui <- function(id) {
  ns <- NS(id)
  tagList(
  textOutput(ns("Fastqc_Stats_Name")),
  tabsetPanel(id='fastqcstatspanel', type = "tabs", 
              tabPanel(title = "General stats", DT::DTOutput(ns("fastqc_stats"))),
              tabPanel(title = "Per base sequence quality", plotOutput(ns("fastqc_sequence_q"), height = '2000px')),
              tabPanel(title = "Overrepresented sequences", DT::DTOutput(ns("fastqc_overrepseq")))),
  textOutput(ns("Aftertrimming")),
  tabsetPanel(id='fastqcstatspanel_AfterTrimming', type = "tabs",
              tabPanel(title = "General stats", DT::DTOutput(ns("fastqc_stats_AfterTrimming"))),
              tabPanel(title = "Per base sequence quality", plotOutput(ns("fastqc_sequence_q_AfterTrimming"), height = '2000px')) ,
              tabPanel(title = "Overrepresented sequences", DT::DTOutput(ns("fastqc_overrepseq_AfterTrimming"))))
  )
}


mod_fastqc_server <- function(
    id,
    rv,
    fastqc_report_before_trimming,
    fastqc_report_after_trimming
) {
  
  moduleServer(id, function(input, output, session) {
    
    output$fastqc_stats <- DT::renderDT({
      req(rv$fastqc_finished == 1)
      
      frbt <- fastqc_report_before_trimming()
      
      basicstats <- lapply(frbt, function(x){
        dedup <- round(100 - x$total_deduplicated_percentage, 2)
        c(as.data.frame(x$basic_statistics)[c(1,4,5,6,7),2], dedup)
      })
      
      basicstats <- do.call(rbind, basicstats)
      
      colnames(basicstats) <- c(
        as.data.frame(frbt[[1]]$basic_statistics)[c(1,4,5,6,7),1],
        "%Duplicated"
      )
      
      basicstats
    })
    
    output$fastqc_overrepseq <- DT::renderDT({
      req(rv$fastqc_finished == 1)
      frbt <- fastqc_report_before_trimming()
      samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(frbt))
      
      overseq<-list(); for (i in 1:length(frbt)){
        overseq[[i]] <- qc_plot(frbt[[i]], "Overrepresented sequences")
        if(is.null(dim(overseq[[i]]))) overseq[[i]] <- tibble(Sequence="none", Count=NULL,Percentage=NULL,"Possible Source"="none", Sample=samplenames[i]) else {
          overseq[[i]] <- overseq[[i]] %>% add_column(Sample=samplenames[i])
        } 
      }
      
      overseq <- overseq %>% bind_rows()
      
      d1 <- DT::datatable(overseq,
                          rownames = FALSE, 
                          extensions = 'RowGroup', 
                          options = list(rowGroup = list(dataSrc=c(2)),
                                         pageLength = 100,
                                         columnDefs = list(list(visible=FALSE, targets=c(2)))
                          ))
      return(d1)
    })
    
    output$fastqc_sequence_q <- renderPlot({
      req(rv$fastqc_finished == 1)
      p_pbsq <- list(); 
      frbt <- fastqc_report_before_trimming()
      samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(frbt))
      p_pbsq[[1]] <- qc_plot(frbt[[1]], 
                             modules = 'Per base sequence quality') + ggtitle(label="Per base sequence quality", 
                                                                              subtitle =samplenames[1] ); 
      for (i in 2:length(frbt))
        p_pbsq[[i]] <- qc_plot(frbt[[i]], modules = 'Per base sequence quality') + ggtitle(label="",subtitle = samplenames[i]); 
      return(Reduce(f = '+', p_pbsq))
    })
    
    
    fastqc_report_after_trimming <- reactive({
      req(rv$salmon_finished == 1)
      print("line 895 was ok")
      fastqc.zip_files <- list.files(path = trimmed_folder(), pattern = "fastqc.zip", full.names = T)
      print("line 897 was ok")
      print(fastqc.zip_files)
      print("line 898 was ok")
      fastqc.zip_files_read <-list()
      for (i in 1:length(fastqc.zip_files)){
        fastqc.zip_files_read[[i]] <- suppressMessages(qc_read(file = fastqc.zip_files[[i]]))
      }
      saveRDS(fastqc.zip_files_read, file = paste0(ProjFolderFull(),'/fastqc_report_after_trimming.RDS'))
      return(fastqc.zip_files_read)
    })
    
    output$fastqc_stats_AfterTrimming <- DT::renderDT({
      req(rv$fastqc_finished == 1)
      req(rv$bbduk_finished == 1)
      basicstats <- data.frame(NULL)
      
      frat <- fastqc_report_after_trimming()
      for (i in 1:length(frat)){
        dedup <- round(100-frat[[i]]$total_deduplicated_percentage,2)
        basicstats <- rbind(basicstats, c(as.data.frame(frat[[i]]$basic_statistics)[c(1,4,5,6,7),2], dedup))
      }
      colnames(basicstats) <- c(as.data.frame(frat[[1]]$basic_statistics)[c(1,4,5,6,7),1], "%Duplicated")
      
      return(basicstats)
    })
    output$fastqc_overrepseq_AfterTrimming <- DT::renderDT({
      req(rv$fastqc_finished == 1)
      req(rv$bbduk_finished == 1)
      frat <- fastqc_report_after_trimming()
      samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(frat))
      
      overseq<-list(); for (i in 1:length(frat)){
        overseq[[i]] <- qc_plot(frat[[i]], "Overrepresented sequences")
        if(is.null(dim(overseq[[i]]))) overseq[[i]] <- tibble(Sequence="none", Count=NULL,Percentage=NULL,"Possible Source"="none", Sample=samplenames[i]) else {
          overseq[[i]] <- overseq[[i]] %>% add_column(Sample=samplenames[i])
        } 
      }
      
      overseq <- overseq %>% bind_rows()
      
      d1 <- DT::datatable(overseq,
                          rownames = FALSE, 
                          extensions = 'RowGroup', 
                          options = list(rowGroup = list(dataSrc=c(2)),
                                         pageLength = 100,
                                         columnDefs = list(list(visible=FALSE, targets=c(2)))
                          )
      )
      return(d1)
    })
    
    output$fastqc_sequence_q_AfterTrimming <- renderPlot({
      req(rv$fastqc_finished == 1)
      req(rv$bbduk_finished == 1)
      p_pbsq <- list(); 
      frat <- fastqc_report_after_trimming()
      samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(frat))
      p_pbsq[[1]] <- qc_plot(frat[[1]], 
                             modules = 'Per base sequence quality') + ggtitle(label="Per base sequence quality", 
                                                                              subtitle =samplenames[1] ); 
      for (i in 2:length(frat))
        p_pbsq[[i]] <- qc_plot(frat[[i]], modules = 'Per base sequence quality') + ggtitle(label="",subtitle = samplenames[i]); 
      return(Reduce(f = '+', p_pbsq))
    })
  })
}



# mod_fastqc_server <- function(
#     id,
#     rv,
#     fastqc_report_before_trimming,
#     fastqc_report_after_trimming
# ) {
#   
#   moduleServer(id, function(input, output, session) {
#     
#     make_basicstats <- function(fr) {
#       basicstats <- lapply(fr, function(x) {
#         dedup <- round(100 - x$total_deduplicated_percentage, 2)
#         c(as.data.frame(x$basic_statistics)[c(1,4,5,6,7),2], dedup)
#       })
#       
#       basicstats <- do.call(rbind, basicstats)
#       
#       colnames(basicstats) <- c(
#         as.data.frame(fr[[1]]$basic_statistics)[c(1,4,5,6,7),1],
#         "%Duplicated"
#       )
#       
#       basicstats
#     }
#     
#     make_overrep <- function(fr) {
#       samplenames <- tools::file_path_sans_ext(
#         extractNamesFromFastqcrReports(fr)
#       )
#       
#       overseq <- lapply(seq_along(fr), function(i) {
#         tmp <- qc_plot(fr[[i]], "Overrepresented sequences")
#         
#         if (is.null(dim(tmp))) {
#           tibble(
#             Sequence = "none",
#             Count = NA,
#             Percentage = NA,
#             `Possible Source` = "none",
#             Sample = samplenames[i]
#           )
#         } else {
#           tmp %>% dplyr::add_column(Sample = samplenames[i])
#         }
#       })
#       
#       dplyr::bind_rows(overseq)
#     }
#     
#     make_sequence_plot <- function(fr) {
#       samplenames <- tools::file_path_sans_ext(
#         extractNamesFromFastqcrReports(fr)
#       )
#       
#       plots <- lapply(seq_along(fr), function(i) {
#         qc_plot(fr[[i]], modules = "Per base sequence quality") +
#           ggtitle(
#             label = ifelse(i == 1, "Per base sequence quality", ""),
#             subtitle = samplenames[i]
#           )
#       })
#       
#       Reduce(`+`, plots)
#     }
#     
#     # ---- BEFORE TRIMMING ----
#     
#     output$fastqc_stats <- DT::renderDT({
#       req(rv$fastqc_finished == 1)
#       fr <- fastqc_report_before_trimming()
#       make_basicstats(fr)
#     })
#     
#     output$fastqc_overrepseq <- DT::renderDT({
#       req(rv$fastqc_finished == 1)
#       fr <- fastqc_report_before_trimming()
#       
#       DT::datatable(
#         make_overrep(fr),
#         rownames = FALSE,
#         extensions = "RowGroup",
#         options = list(
#           rowGroup = list(dataSrc = 2),
#           pageLength = 100,
#           columnDefs = list(list(visible = FALSE, targets = 2))
#         )
#       )
#     })
#     
#     output$fastqc_sequence_q <- renderPlot({
#       req(rv$fastqc_finished == 1)
#       fr <- fastqc_report_before_trimming()
#       make_sequence_plot(fr)
#     })
#     
#     # ---- AFTER TRIMMING ----
#     
#     output$fastqc_stats_AfterTrimming <- DT::renderDT({
#       req(rv$fastqc_finished == 1)
#       req(rv$bbduk_finished == 1)
#       fr <- fastqc_report_after_trimming()
#       make_basicstats(fr)
#     })
#     
#     output$fastqc_overrepseq_AfterTrimming <- DT::renderDT({
#       req(rv$fastqc_finished == 1)
#       req(rv$bbduk_finished == 1)
#       fr <- fastqc_report_after_trimming()
#       
#       DT::datatable(
#         make_overrep(fr),
#         rownames = FALSE,
#         extensions = "RowGroup",
#         options = list(
#           rowGroup = list(dataSrc = 2),
#           pageLength = 100,
#           columnDefs = list(list(visible = FALSE, targets = 2))
#         )
#       )
#     })
#     
#     output$fastqc_sequence_q_AfterTrimming <- renderPlot({
#       req(rv$fastqc_finished == 1)
#       req(rv$bbduk_finished == 1)
#       fr <- fastqc_report_after_trimming()
#       make_sequence_plot(fr)
#     })
#     
#   })
# }
