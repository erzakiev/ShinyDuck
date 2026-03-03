#' @import shiny
#' @import dplyr

app_ui <- function() {

  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = span("ShinyDuck🦆",
                   style = "color: white; font-size: 22px; font-weight: bold")
    ),

    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(id = "sidebar",
                                  shinydashboard::menuItem("Setup",
                                                           tabName = "setup",
                                                           icon = icon("play")),
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
                                  shinydashboard::menuItem("GSVA",
                                                           tabName = "GSVA",
                                                           icon = icon("chart-column")),
                                  shinydashboard::menuItem("GO enrichment",
                                                           tabName = "enrich",
                                                           icon = icon("sitemap")),
                                  shinydashboard::menuItem("Downloads",
                                                           tabName = "download",
                                                           icon = icon("download"))

      )
    ),

    shinydashboard::dashboardBody(
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
      "))
      ),


      shinyjs::useShinyjs(),

      shinydashboard::tabItems(

        shinydashboard::tabItem(tabName = "setup",
                mod_project_setup_ui("setup")
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
        )#,
        #shinydashboard::tabItem(tabName = "GSVA",
        #        mod_gsva_ui("gsva")
        #),
        #shinydashboard::tabItem(tabName = "enrich",
        #        mod_enrich_ui("enrich")
        #),
        #shinydashboard::tabItem(tabName = "download",
        #        mod_download_ui("download")
        #)
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
