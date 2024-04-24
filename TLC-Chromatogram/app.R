library(shiny)
library(bslib)
library(rhandsontable)
library(DT)
library(tidyverse)
library(ggplot2)


##FUNCTIONS-------------
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

#-----------
#UI---------
ui <- fluidPage(
  title = "Chromatogram Simulator",
  page_fluid(
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
        
        # Experimental Isocratic Mode ----
        nav_panel(
          "Empirical Rf",
          id = "empirical_isocratic",
          numericInput("crude_mass", "Crude Mass (g)", 0.5, min = 0),
          rHandsontableOutput("rftable"),
          verbatimTextOutput("debugger")
        ),
        
        # Prediction Isocratic Mode ----
        nav_panel(
          "Predicted Rf",
          id = "predicted"
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
  )
)


#-------------
#SERVER-------
server <- function(input, output, session) {
  
  ##USER TLC DATA INPUT
  # Default TLC data for editable table
  default_tlc <- reactive(data.frame(
    "Rf" = c(0.71, 0.4, 0.2),
    "Xa" = c(0.2, 0.5, 0.3)
  ))
  # Render the editable Rf table
  output$rftable <- renderRHandsontable({
    rhandsontable(default_tlc(), selectCallback=TRUE, readOnly = FALSE)
  })
  
  #--------------------
  #CALCULATIONS -------
  
  
  rv <- reactiveValues(tlc_data = NULL,
                       is_hard = NULL,
                       efficiency = NULL,
                       silica_mass = NULL,
                       Vm = NULL,
                       Vr = NULL,
                       bandsize = NULL,
                       plotdata = data.frame("solvent" = rep(0, length = 500))
                       )
  
  observe({ rv$tlc_data <- hot_to_r(input$rftable)
    rv$is_hard <- any(diff(sort(rv$tlc_data$Rf)) < 0.2)
    #Calculate efficiency of each spot (N) and column silica mass
    if (rv$is_hard == TRUE) {
      rv$efficiency <- N_hard(rv$tlc_data$Rf)
      rv$silica_mass <- siog_hard(input$crude_mass)
    } else {
      rv$efficiency <- N_easy(rv$tlc_data$Rf)
      rv$silica_mass <- siog_easy(input$crude_mass)
    }
    #Calculate void volume (Vm), retention volume (Vr), bandsize
    rv$Vm <- calc_Vm(rv$silica_mass)
    rv$Vr <- calc_Vr(rv$Vm, rv$tlc_data$Rf)
    rv$bandsize <- calc_bandsize(rv$Vr, rv$efficiency)
    
    rv$plotdata <- data.frame("solvent" = seq(0, rv$Vm * 5, length = 500))
    if (!is.null(nrow(rv$tlc_data))){
      for (i in 1:nrow(rv$tlc_data)) { 
        # Create column names
        spot_col <- paste("Spot", i, sep = "_")
        
        # Assign values to new column
        rv$plotdata[spot_col] <- elution_curve(rv$plotdata$solvent, rv$Vr[i], rv$tlc_data$Xa[i], rv$bandsize[i])
      }
    #Convert to long format for ggplot2
    rv$plotdata <- pivot_longer(rv$plotdata, cols = starts_with("Spot"))
    }
  }) %>%
  bindEvent({ #make it update only when the table or crude_mass is updated
    input$rftable
    input$crude_mass
    }, label = "calculations")
  
  # output$debugger <- renderText(
  #   paste(str(rv$plotdata))
  # )
  #--------------
  #PLOTTING------
  output$chromatogram <- renderPlot({
    ggplot(data=rv$plotdata, aes(x = solvent, y = value, colour = name)) +
        geom_line() +
        labs( x = "Eluent (mL)", y = "Relative peak intensity")
  })
  
}

# Call the shinyapp
shinyApp(ui = ui, server = server)
