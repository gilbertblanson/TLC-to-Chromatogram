library(shiny)
library(bslib)
library(DT)

# Define UI
ui <- fluidPage(# Page title
  title = "Chromatogram Simulator",
  
  # Navbar
  page_fluid(
    title = "Chromatogram Simulator",
    
    # Content
    layout_columns(
      # Chromatogram
      card(card_header(h4("Chromatogram")),
           plotOutput(outputId = "chromatogram")),
      
      # Other information card
      navset_card_underline(
        #Panel with TLC plot
        nav_panel("TLC", plotOutput(outputId = "tlc")),
        
        #Panel with Resolution Data (Rs)
        nav_panel(HTML(paste0("R", tags$sub(
          "s"
        ))), verbatimTextOutput("resolution")),
        
        #Panel with Elution Information
        nav_panel("Peak Data"),
        
        #Panel with Fraction information
        nav_panel("Fraction Data")
      ),
      
      # Data Input
      navset_card_underline(
        title = "Data entry",
        
        # If the user wants to plot experimental Rfs ----
        nav_panel(
          "Empirical Rf",
          id = "empirical_isocratic",
          DTOutput("rf_table"),
          layout_columns(
            actionButton("add_row_rf_table", "Add Row"),
            actionButton("delete_row_rf_table", "Delete Row")
          )
        ),
        
        # If the user wants to use predicted Rfs ----
        nav_panel(
          "Predicted Rf",
          id = "predicted",
          sliderInput(
            "strong_solvent_percent",
            "Eluent, Strong solvent%",
            min = 0,
            max = 100,
            value = 1,
            ticks = FALSE,
            animate = TRUE
          ),
          p(
            "Step 1: Provide the necessary Rf data for each spot. Each Rf column represents a different strong solvent%. Add at least three data points for each spot for accurate prediction!"
          ),
          DTOutput("prediction_table"),
          layout_columns(
            actionButton("add_row_pred_table", "Add Row"),
            actionButton("delete_row_pred_table", "Delete Row")
          ),
          p("Step 2: Provide corresponding strong solvent% values for the Rfs."),
          layout_columns(
            numericInput("rf_1", "Rf_1", value = 0),
            numericInput("rf_2", "Rf_2", value = 0),
            numericInput("rf_3", "Rf_3", value = 0),
            numericInput("rf_4", "Rf_4", value = 0),
            numericInput("rf_5", "Rf_5", value = 0)
          )
        )
      ),
      
      # Optional Settings
      card(
        card_header(h4("Optional Settings")),
        p(
          "Change the values below if you wish to set your own column size or fraction size. Otherwise, set to 0."
        ),
        numericInput("user_silica", "Silica (g):", value = 0),
        numericInput("user_fraction_size", "Fraction size (mL):", value = 0)
      ),
      col_widths = c(8, 4, 8, 4)
    )
  ))

# Define server logic
server <- function(input, output, session) {
  # Chromatogram
  output$chromatogram <- renderPlot({
    plot(c(1, 2, 3, 4), c(2, 3, 4, 5))
  })
  
  ##USER TLC DATA
  # Create a reactive dataframe to store the table data
  user_tlc_data <- reactiveVal(data.frame(
    "Spot_Number" = c(1, 2, 3),
    "Rf" = c(0.5, 0.6, 0.7),
    "Xa" = c(0.2, 0.2, 0.6)
  ))
  
  # Render the editable Rf table
  output$rf_table <- renderDT({
    # Display the reactive dataframe as an editable table
    datatable(
      user_tlc_data(),
      editable = TRUE,
      rownames = FALSE,
      selection = 'single',
      options = list(dom = 't', ordering = F)
    )
  })
  
  # Update the dataframe when the user edits the table
  observeEvent(input$rf_table_cell_edit, {
    info <- input$rf_table_cell_edit
    
    # Update the reactive dataframe based on user edits
    newData <- user_tlc_data()
    newData[info$row, info$col] <- info$value
    user_tlc_data(newData)
  })
  
  ##USER DATA FOR RF PREDICION
  user_prediction_data <- reactiveVal(
    data.frame(
      "Spot_Number" = c(1, 2, 3),
      "Rf_1" = c(0.1, 0.1, 0.2),
      "Rf_2" = rep(NA, 3),
      "Rf_3" = rep(NA, 3),
      "Rf_4" = rep(NA, 3),
      "Rf_5" = rep(NA, 3),
      "Xa" = c(0.2, 0.2, 0.4)
    )
  )
  
  output$prediction_table <- renderDT({
    # Display the reactive dataframe as an editable table
    datatable(
      user_prediction_data(),
      editable = TRUE,
      rownames = FALSE,
      selection = 'single',
      options = list(dom = 't', ordering = F)
    )
  })
}

# Call the shinyapp
shinyApp(ui = ui, server = server)
