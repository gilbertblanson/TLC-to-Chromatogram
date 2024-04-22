library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(ggplot2)

##FUNCTIONS##
#Easy: when deltaRf >= 0.2; Hard: when deltaRf < 0.2
#Efficiency, N
N_easy <- function(mass_loading){
  return(33.64 * mass_loading ^ (-0.44))
}
N_hard <- function(mass_loading){
  return(51.7 * mass_loading ^ (-0.44) )
}
# Total silica mass (SiO_g)
siog_easy <- function(crude_mass){
  return(59.8 * crude_mass)
}
siog_hard <- function(crude_mass){
  return(151.2 * crude_mass + 0.5)
}

# Column void volume (Vm)
calc_Vm <- function(siog){
  return(1.81*siog+0.31)
}
# Retention volume (Vr)
calc_Vr <- function(Vm, rf) {
  return(Vm * (1 + ((1 - rf) / rf)) * 0.64)
}
# bandsize (Vb)
calc_bandsize <- function(Vm, efficiency){
  Vb <- 4 * Vm/sqrt(efficiency)
  return(Vb)
}
#Elution curve calculation (y)
elution_curve <- function(x, Vr, Xa, bandsize) {
  Xa * exp(-(x - Vr)^2/(bandsize/2.8)^2)
}

# Define UI
ui <- fluidPage(
  # Page title
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
          numericInput("crude_mass", "Crude Mass (g)", 0.5, min = 0),
          DTOutput("rf_table"),
          layout_columns(
            actionButton("add_row_rf_table", "Add Row"),
            actionButton("delete_row_rf_table", "Delete Row")
          ),
          verbatimTextOutput("results")
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
  
  # ##USER DATA FOR RF PREDICION
  # user_prediction_data <- reactiveVal(
  #   data.frame(
  #     "Spot_Number" = c(1, 2, 3),
  #     "Rf_1" = c(0.1, 0.1, 0.2),
  #     "Rf_2" = rep(NA, 3),
  #     "Rf_3" = rep(NA, 3),
  #     "Rf_4" = rep(NA, 3),
  #     "Rf_5" = rep(NA, 3),
  #     "Xa" = c(0.2, 0.2, 0.4)
  #   )
  # )
  
  # output$prediction_table <- renderDT({
  #   # Display the reactive dataframe as an editable table
  #   datatable(user_prediction_data(), editable = TRUE, rownames = FALSE, selection = 'single', options = list(dom = 't', ordering = F))
  # })
  
  ##USER TLC DATA INPUT
  # Create a reactive dataframe to store the table data
  user_tlc_data <- reactiveVal(data.frame(
    "Rf" = c(0.7, 0.5, 0.25),
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
  
  #CALCULATIONS AND PLOTTING---
  
  plotdf <- reactive({
    isolate({
      #make a copy of the user-inputted data to perform manipulations on
      new_data <- isolate(user_tlc_data())
      
      # Calculate if the separation is hard
      is_hard <- reactiveVal(any(diff(sort(user_tlc_data()$Rf)) < 0.2))

      # Calculate efficiency of each spot (N) and column silica mass
      if (is_hard() == TRUE) {
        new_data$efficiency <- N_hard(new_data$Rf)
        silica_mass <- siog_hard(input$crude_mass)
      } else {
        new_data$efficiency <- N_easy(new_data$Rf)
        silica_mass <- siog_easy(input$crude_mass)
      }

      # Calculate void volume
      Vm <- calc_Vm(silica_mass)

      # Calculate retention volume
      new_data$Vr <- calc_Vr(Vm, new_data$Rf)

      # Calculate band size
      new_data$bandsize <- calc_bandsize(new_data$Vr, new_data$efficiency)

      # Create X variables
      plotdf <- data.frame("solvent" = seq(0, Vm * 5, length = 500))

      # Calculate Spot columns
      for (i in 1:nrow(user_tlc_data())) {
        # Create column names
        spot_col <- paste("Spot", i, sep = "_")

        # Calculate values using the formula for current row
        values <- elution_curve(plotdf$solvent, new_data$Vr[i], new_data$Xa[i], new_data$bandsize[i])

        # Assign values to new column
        plotdf[[spot_col]] <- values
      }

      # Convert to long format so that it's ggplot compatible
      #Should turn it from
      # solvent spot 1 spot 2 spot 3
      #   0       0     0.1     0
      #   1       0     0.2     0.1 
      # To
      # solvent name  value
      #   0     spot1   0
      #   0     spot2   0.1
      #   0     spot3   0   
      #   1     spot1   0 
      # etc etc
      plotdf <- pivot_longer(plotdf, cols = starts_with("Spot"))
      
      return(plotdf)
    })
  })
  
  #Output to see what the datatable looks like. For some reason, renderDT doesn't work on plotdf(), not sure what's happening
  output$results <- renderText({paste(head(plotdf()))})
  
  # View the resulting dataframe
  output$chromatogram <- renderPlot({
    #Calculations
    
    ggplot(data=plotdf(),
           aes(x=solvent, y=value, colour=name)) +
      geom_line() +
      labs( x = "Eluent (mL)", y = "Relative peak intensity")
  })
  
}

# Call the shinyapp
shinyApp(ui = ui, server = server)
