#' @import shiny
#' @import dplyr
#' @import ggplot2

app_ui <- function() {

  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = span("ShinyDucks🦆🦆",
                   style = "color: white; font-size: 25px; font-weight: bold")
    ),

    shinydashboard::dashboardSidebar(
      tags$head(tags$style(HTML(".content-wrapper {overflow-x: scroll;}"))),
      shinydashboard::sidebarMenu(id = "sidebar",
                                  shinydashboard::menuItem("Setup",
                                                           tabName = "setup",
                                                           icon = icon("play")),
                                  shinydashboard::menuItem("Design & Subsetting",
                                                           tabName = "coldata",
                                                           icon = icon("bezier-curve")),
                                  shinydashboard::menuItem("QC",
                                                           tabName = "qc",
                                                           icon = icon("medal")),
                                  shinydashboard::menuItem("TPM Table",
                                                           tabName = "tpm",
                                                           icon = icon("table")),
                                  shinydashboard::menuItem("DEGs",
                                                           tabName = "deg",
                                                           icon = icon("yin-yang")),
                                  shinydashboard::menuItem("PCA",
                                                           tabName = "pca",
                                                           icon = icon("arrow-up-right-dots")),
                                  shinydashboard::menuItem("Volcano",
                                                           tabName = "volcano",
                                                           icon = icon("mountain")),
                                  shinydashboard::menuItem("GO enrichment",
                                                           tabName = "enrich",
                                                           icon = icon("sitemap")),
                                  shinydashboard::menuItem("GSVA",
                                                           tabName = "GSVA",
                                                           icon = icon("chart-column")),
                                  shinydashboard::menuItem("Downloads",
                                                           tabName = "download",
                                                           icon = icon("download"))

      )
    ),

    shinydashboard::dashboardBody(

      tags$script(HTML('
                           $(document).ready(function() {
                           $("header").find("nav").append(\'<div id="pageHeader" class="myClass"></div>\');
                           })
                           ')),
      tags$script(HTML('$("body").addClass("fixed");')),


      tags$head(
        tags$style(HTML("
        .disabled-menu-item {
          pointer-events: none;
          opacity: 0.5;
          background-color: #222d33;
        }
        .disabled-menu-item a {
          color: #999 !important;
        }
      ")),
        tags$style(HTML("table.dataTable tbody tr {height: 40px;}")),
      ),

      tags$head(tags$style(HTML(
        '.myClass {
            font-size: 20px;
            line-height: 50px;
            text-align: left;
            font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
            padding: 0 15px;
            overflow: hidden;
            color: white;
            };
        .small-box {height: 250px}
            '))),




      shinyjs::useShinyjs(),

      shinydashboard::tabItems(

        shinydashboard::tabItem(tabName = "setup",
                mod_project_setup_ui("setup")
        ),

        shinydashboard::tabItem(tabName = "coldata",
                                mod_coldata_ui("coldata")
        ),

        shinydashboard::tabItem(tabName = "qc",
                mod_fastqc_ui("fastqc")
        ),

        shinydashboard::tabItem(tabName = "tpm",
                mod_tpm_table_ui("tpm")
        ),

        shinydashboard::tabItem(tabName = "deg",
                mod_deseq_deg_ui("deg")
        ),

        shinydashboard::tabItem(tabName = "pca",
                mod_pca_ui("pca")
        ),

        shinydashboard::tabItem(tabName = "volcano",
                mod_volcano_ui("volcano")
        ),
        shinydashboard::tabItem(tabName = "enrich",
                mod_enrich_ui("enrich")
        ),
        shinydashboard::tabItem(tabName = "GSVA",
                mod_gsva_ui("gsva")
        ),
        shinydashboard::tabItem(tabName = "download",
                                        mod_download_ui("download")
                                )
    ),

    tags$script(HTML("
      Shiny.addCustomMessageHandler('toggle_menu_items', function(message) {
        var menuItems = document.querySelectorAll('.sidebar-menu li');
        menuItems.forEach(function(item, index) {
          // Skip the first item or specific indices
          if(index !== 0) {  // Disable all except first item
            if(message.disable) {
              item.classList.add('disabled-menu-item');
            } else {
              item.classList.remove('disabled-menu-item');
            }
          }
        });
      });
    "))
  )

  )
}
