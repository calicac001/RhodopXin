# References

# Attali D (2021). _shinyjs: Easily Improve the User Experience of Your Shiny Apps in
# Seconds_. R package version 2.1.0, <https://CRAN.R-project.org/package=shinyjs>.

# Attali D, Sali A (2024). _shinycssloaders: Add Loading Animations to a 'shiny' Output While
# It's Recalculating_. R package version 1.1.0,
# <https://CRAN.R-project.org/package=shinycssloaders>.

# Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A,
# Borges B (2024). _shiny: Web Application Framework for R_. R package version 1.9.1,
# <https://CRAN.R-project.org/package=shiny>.

# Chang W (2021). _shinythemes: Themes for Shiny_. R package version 1.2.0,
# <https://CRAN.R-project.org/package=shinythemes>.

# OpenAI. (2024). ChatGPT (November 2024 version). Retrieved from
# https://chat.openai.com

library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinythemes)
library(RhodopXin)

# Define UI for random distribution app ----
ui <- fluidPage(
  title = "RhodopXin",
  theme = shinythemes::shinytheme("united"),
  h1("RhodopXin",
     style = "text-align: center;
              background-color: #e95420;
              color: white;
              padding: 20px;
              border-radius: 5px;"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    sidebarPanel(
      tags$p("Performs pairwise alignment between each of the template
             rhodopsin's helices and each of the query sequences given as input
             assessing the conservation of residues across a given helix.
             See ?createHelixAlignments for more info."),

      # Header for template rhodopsin and instructions
      tags$h3("Template Rhodopsin"),
      tags$p("Choose a template rhodopsin below (see ?template_rhodopsins for more
             info) or provide the RCSB PDB id of the rhodopsin you want to use."),

      # Radio buttons for selecting a template rhodopsin
      radioButtons(
        inputId = "template_type",
        label = "Choose an option:",
        choices = c("Sample Templates", "Choose Own Template"),
        selected = "Sample Templates"
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
                      "Xanthorhodopsin (3DDL)"),
          selected = "Bacteriorhodopsin (1QHJ)"
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
             query rhodopsins provided in package (see ?sample_rhodopsins for more info)."),

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
          accept = c(".fasta", ".fa", ".txt"))
      ),

      # Button to run alignment
      actionButton(inputId = "run_alignment", label = "Run Helix Alignments",
                   class = "btn-primary", style = "width: 100%;")
    ),

    # Main panel for displaying outputs
    mainPanel(
      tabsetPanel(
        tabPanel("Alignments",
                 div(
                   id = "loading",
                   withSpinner(plotOutput("helix_plot"), type = 4, color = "#e95420"))
                 ),
        tabPanel("Mapping", uiOutput("query_mapping"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Setup reactive values for input
  template_rv <- reactiveVal(NULL)
  observe({
    if (input$template_type == "Sample Templates"){
      template_rv(switch(input$template,
                         "Bacteriorhodopsin (1QHJ)" = RhodopXin::template_rhodopsins[1],
                         "Channelrhodopsin (3UG9)" = RhodopXin::template_rhodopsins[2],
                         "Halorhodopsin (3A7K)" = RhodopXin::template_rhodopsins[3],
                         "Proteorhodopsin (4JQ6)" = RhodopXin::template_rhodopsins[4],
                         "Xanthorhodopsin (3DDL)" = RhodopXin::template_rhodopsins[5],
                         RhodopXin::template_rhodopsins[1]))
    } else if (input$template_type == "Choose Own Template"){
      req(input$template)
      template_rv(RhodopXin::loadFromRCSB(input$template))
    }
  })

  query_rv <- reactiveVal(NULL)
  observe( {
    if (input$query_type == "Sample Rhodopsins"){
      query_rv(RhodopXin::sample_rhodopsins)
    } else if (input$query_type == "Choose Own Rhodopsin/s"){
      req(input$file)
      query_rv(RhodopXin::loadSequence(input$file))
    }
  })

  rcsb_id_rv <- reactiveVal(NULL)
  observe( {
    if (input$template_type == "Sample Templates"){
      rcsb_id_rv(switch(input$template,
                        "Bacteriorhodopsin (1QHJ)" = "1QHJ",
                        "Channelrhodopsin (3UG9)" = "3UG9",
                        "Halorhodopsin (3A7K)" = "3A7K",
                        "Proteorhodopsin (4JQ6)" = "4JQ6",
                        "Xanthorhodopsin (3DDL)" = "3DDL",
                        "1QHJ"))
    } else if (input$template_type == "Choose Own Template"){
      rcsb_id_rv(input$template)
    }
  })

  # Expression for when alignment buttion is clicked
  startAlignment <- eventReactive(input$run_alignment, {
    req(template_rv())
    req(query_rv())
    req(rcsb_id_rv())

    shinyjs::show("loading")

    RhodopXin::createHelixAlignments(template = template_rv(),
                                     sequences = query_rv(),
                                     rcsb_id = rcsb_id_rv())
  })

  # Determine height dynamically
  plot_height <- reactive({
    if (!is.null(startAlignment())) {
      results <- startAlignment()
      all_pwa <- results$all_pwa
      return(ceiling(length(all_pwa) / 2) * 200)  # Calculate height dynamically
    }
    return(500)  # Default height if no data
  })

  # Render the helix alignment plot
  output$helix_plot <- renderPlot({
    if(!is.null(startAlignment())){
      results <- startAlignment()
      all_pwa <- results$all_pwa
      plot <- RhodopXin::visualizeHelixAlignments(all_pwa = all_pwa)
      shinyjs::hide("loading")
      return(plot)
    }
  }, height = function() {
    plot_height()  # Pass the reactive plot height
  }
  )

  # Dropdown to choose query for 3d mapping
  output$query_mapping <- renderUI({
    if (!is.null(startAlignment())) {
      results <- startAlignment()

      colnames_tr <- colnames(results$template_ranges)
      query_columns <- grep("^Query", colnames_tr, value = TRUE)
      num_query <- length(query_columns)

      tagList(
        selectInput(inputId = "query_to_map",
                    label = "Choose Query to Map",
                    choices = 1:num_query,
                    selected = 1),
        r3dmol::r3dmolOutput("mapping_structure")
      )

    }
  })

  # Render the highlighted 3d structure depending on query input
  output$mapping_structure <- r3dmol::renderR3dmol({
    if (!is.null(startAlignment())){
      results <- startAlignment()
      template_ranges <- results$template_ranges

      RhodopXin::visualizeHelixMapping(template_ranges = template_ranges,
                                       template = template_rv(),
                                       rcsb_id = rcsb_id_rv(),
                                       query_num = as.numeric(input$query_to_map))
    }
  })
}

# Create Shiny app
shinyApp(ui, server)

#[END]
