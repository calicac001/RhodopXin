# This example is adapted from
# RStudio Inc. (2013). Tabsets. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/tabsets.html

library(shiny)
library(shinythemes)
library(RhodopXin)

# Define UI for random distribution app ----
ui <- navbarPage(
  title = "RhodopXin",
  theme = shinythemes::shinytheme("united"),

  tabPanel("Helix Alignments",
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      sidebarPanel(
        tags$p("Performs pairwise alignment between each of the template
               rhodopsin's helices and each of the query sequences given as input
               assessing the conservation of residues across a given helix.
               See ?createHelixAlignments for more info."),

        # Header for template rhodopsin and instructions
        tags$h3("Template Rhodopsin"),
        tags$p("Choose a template rhodopsin below (see Sample Data for more
               info) or provide the RCSB PDB id of the rhodopsin you want to use."),

        # Radio buttons for selecting a template rhodopsin
        radioButtons(
          inputId = "template_type",
          label = "Choose an option:",
          choices = c("Sample Templates", "Choose Own Template")
        ),

        # Conditional dropdown for the sample rhodopsins
        conditionalPanel(
          condition = "input.template_type == 'Sample Templates'",  # Show only when 'Option B' is selected
          selectInput(
            inputId = "template",
            label = "Choose a template rhodopsin from samples:",
            choices = c("Bacteriorhodopsin (1QHJ)",
                        "Channelrhodopsin (3UG9)",
                        "Halorhodopsin (3A7K)",
                        "Proteorhodopsin (4JQ6)",
                        "Xanthorhodopsin (3DDL)")
          )
        ),

        # Conditional text input for selecting own template
        conditionalPanel(
          condition = "input.template_type == 'Choose Own Template'",  # Show only when 'Option B' is selected
          textInput(
            inputId = "template",
            label = "Enter RCSB PDB id:",
            placeholder = " A 4-character alphanumeric identifier"
          )
        ),

        # Header for query rhodopsin and instructions
        tags$h3("Query Rhodopsin"),
        tags$p("Attach a FASTA file containing rhodopsin sequences or use sample
               query rhodopsins provided in package (see Sample Data for more info)."),

        # Radio buttons for selecting query rhodopsins
        radioButtons(
          inputId = "query_type",
          label = "Choose an option:",
          choices = c("Sample Rhodopsins", "Choose Own Rhodopsin/s")
        ),


        # Conditional text input for selecting own queries
        conditionalPanel(
          condition = "input.query_type == 'Choose Own Rhodopsin/s'",
          fileInput(
            inputId = "file",
            label = "Choose FASTA File",
            accept = c("text/csv", "text/comma-separated-values","text/plain",
                       ".csv", ".txt"))
        ),

        # actionButton
        actionButton(inputId = "run_alignment", label = "Run Helix Alignments",
                     class = "btn-primary", style = "width: 100%;")
      ),

      # Main panel for displaying outputs
      mainPanel(
        fluidRow(
          # Two halves on top
          column(6, uiOutput("query_and_helix")),
          column(6, uiOutput("helix_struct"))
        ),
        fluidRow(
          # Full-width bottom section
          column(12, uiOutput("helix_plot"))
        )

      )
    )
  )
)


# Define server logic for random distribution app ----
server <- function(input, output) {
  startAlignment <- eventReactive(input$run_alignment, {
    template <- NULL
    if (input$template_type == "Sample Templates"){
      template <- switch(input$template,
                         "Bacteriorhodopsin (1QHJ)" = RhodopXin::sample_rhodopsins[1],
                         "Channelrhodopsin (3UG9)" = RhodopXin::sample_rhodopsins[2],
                         "Halorhodopsin (3A7K)" = RhodopXin::sample_rhodopsins[3],
                         "Proteorhodopsin (4JQ6)" = RhodopXin::sample_rhodopsins[4],
                         "Xanthorhodopsin (3DDL)" = RhodopXin::sample_rhodopsins[5])
    } else if (input$template_type == "Choose Own Template"){
      template <- RhodopXin::loadFromRCSB(input$template)
    }

    query <- NULL
    if (input$query_type == "Sample Rhodopsins"){
      query <- RhodopXin::sample_rhodopsins
    } else if (input$query_type == "Choose Own Rhodopsin/s"){
      query <- RhodopXin::loadSequence(input$file)
    }

    rcsb_id <- NULL
    if (input$template_type == "Sample Templates"){
      rcsb_id <- switch(input$template,
                         "Bacteriorhodopsin (1QHJ)" = "1QHJ",
                         "Channelrhodopsin (3UG9)" = "3UG9",
                         "Halorhodopsin (3A7K)" = "3A7K",
                         "Proteorhodopsin (4JQ6)" = "4JQ6",
                         "Xanthorhodopsin (3DDL)" = "3DDL")
    } else if (input$template_type == "Choose Own Template"){
      rcsb_id <- input$template
    }

    RhodopXin::createHelixAlignments(template = template,
                                     sequences = query,
                                     rcsb_id = rcsb_id)
  })

  output$query_and_helix <- renderPrint({
    if(!is.null(startAlignment))
      startAlignment()

  })
}

# Create Shiny app ----
shinyApp(ui, server)

#[END]
