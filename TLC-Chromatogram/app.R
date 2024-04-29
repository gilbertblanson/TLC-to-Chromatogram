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
  return(33.64 * mass_loading ^ (-0.44) )
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
calc_Vm <- function(SiOg){
  return(1.8*SiOg+0.3)
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
  Xa * exp( -(x - Vr)^2 / (bandsize/2.8)^2 )
}
#gen layers for plot
line_layer_add <- function(plot, funct, label){
  plot + geom_function(fun = funct, aes(colour = label), n = 500)
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
        nav_panel( "TLC", plotOutput(outputId = "tlc") ),
        nav_panel( #HTML( paste0("R", tags$sub("s") ) )
                  "Peak Info", 
                  rHandsontableOutput("resolution")
                  )
      ),
      
      # Data Input
      navset_card_underline(
        title = "Data entry",
        
        # Experimental Isocratic Mode ----
        nav_panel(
          "Empirical Rf",
          id = "empirical_isocratic",
          numericInput("crude_mass", "Crude Mass (g)", 0.5, min = 0),
          rHandsontableOutput("rftable")
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
        verbatimTextOutput("silica"),
        numericInput("user_fraction_size", "Fraction size (mL):", value = 0, min = 0),
        verbatimTextOutput("fraction_size")
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
                       silica_mass = NULL,
                       Vm = NULL,
                       plotdata = data.frame("solvent" = rep(0, length = 500)),
                       sep = NULL,
                       comb_plot = NULL,
                       fraction_size = NULL
                       )
  
  observe({ rv$tlc_data <- hot_to_r(input$rftable)
    rv$is_hard <- any(diff(sort(rv$tlc_data$Rf)) < 0.2)
    #Calculate efficiency of each spot (N) and column silica mass
    if (rv$is_hard == TRUE) {
      rv$tlc_data$efficiency <- N_hard(rv$tlc_data$Xa)
      rv$silica_mass <- siog_hard(input$crude_mass)
    } else {
      rv$tlc_data$efficiency <- N_easy(rv$tlc_data$Xa)
      rv$silica_mass <- siog_easy(input$crude_mass)
    }
    #if user supplies silica amt
    if (!is.na(input$user_silica) & input$user_silica != 0){
      rv$silica_mass <- input$user_silica
    }
    #Calculate void volume (Vm), retention volume (Vr), bandsize
    rv$Vm <- calc_Vm(rv$silica_mass)
    rv$tlc_data$Vr <- calc_Vr(rv$Vm, rv$tlc_data$Rf)
    rv$tlc_data$bandsize <- calc_bandsize(rv$tlc_data$Vr, rv$tlc_data$efficiency)
    
    #make sure the rows are calculated before this is done
    if (!is.null(nrow(rv$tlc_data))){
      rv$sep <- data.frame("Rs" = 2*( rv$tlc_data$Vr - lag(rv$tlc_data$Vr) )/( rv$tlc_data$bandsize + lag(rv$tlc_data$bandsize) )) %>%
        mutate( Separation =  case_when(
          is.na(Rs) == TRUE ~ NA, 
          Rs > 1.5 ~ "Good", 
          Rs > 0.8 ~ "Moderate", 
          .default = "Bad")
          ) 
      rv$sep <- data.frame(
        "Peak_Start_mL" = rv$tlc_data$Vr - rv$tlc_data$bandsize/2,
        "Peak_Middle_mL" = rv$tlc_data$Vr,
        "Peak_End_mL" = rv$tlc_data$Vr + rv$tlc_data$bandsize/2
        ) %>%
        bind_cols(rv$sep, .)
      #to-do: start middle end peak
      if (!is.na(input$user_fraction_size)) {
        rv$fraction_size <- case_when(input$user_fraction_size != 0 ~ input$user_fraction_size, .default = rv$Vm/3)
      } else { rv$fraction_size <- rv$Vm/3 }
      
      rv$sep <- mutate(rv$sep, Fraction_Start = Peak_Start_mL/rv$fraction_size, 
                       Fraction_Mid = Peak_Middle_mL/rv$fraction_size, 
                       Fraction_End = Peak_End_mL/rv$fraction_size,
                       Fraction_sep = Fraction_Start - lag(Fraction_End))
    }
    
  }) %>%
  bindEvent({ #make it update only when the table/crude_mass/user_silica is updated
    input$rftable
    input$crude_mass
    input$user_silica
    input$user_fraction_size
    }, label = "calculations")

  output$resolution <- renderRHandsontable(
    #calculate fractions input$user_fraction_size
    #output start middle end fractions
    rhandsontable(rv$sep, readOnly = TRUE)
  )
  output$silica <- renderText(rv$silica_mass)
  output$fraction_size <- renderText(rv$fraction_size)
  #--------------
  #PLOTTING------
  output$tlc <- renderPlot({
    no_spots <- nrow(rv$tlc_data)
    tlc <- ggplot(data=rv$tlc_data, aes(x = rep(0, no_spots), y = Rf, colour = as.factor(1:no_spots) )) +
      geom_point(size=5) +
      ylim(0,1) +
      scale_x_discrete(labels = NULL, breaks = NULL) +
      labs(x = NULL) +
      theme(legend.position="none", panel.grid.minor=element_blank())
    tlc + scale_colour_brewer(palette = "Set1")
  })
  
  output$chromatogram <- renderPlot({
    if (!is.null(rv$tlc_data$bandsize)){
      # Generate list of functions to plot with all but x already applied
      partials_list <- mapply(FUN=partial,
                              Vr=rv$tlc_data$Vr,
                              Xa=rv$tlc_data$Xa,
                              bandsize=rv$tlc_data$bandsize,
                              MoreArgs=list(`.f`=elution_curve))
      
      # Data to plot. Defining the domain for the plotted functions.
      plot_data <- data.frame(x = seq(0, rv$Vm * 10, length.out = 10))
      
      # Initial plot object
      init_plot <- ggplot(plot_data, aes(x=x))
      
      # Generate layers and accumulate to plot object with label
      comb_plot <- reduce(.x = seq_along(partials_list),
                          .f = function(p, i) line_layer_add(p, partials_list[[i]], paste0("Spot ", i)),
                          .init = init_plot)
      comb_plot +
        labs( x = "Eluent (mL)", y = "Relative peak intensity")
    }
  })
  
}

# Call the shinyapp
shinyApp(ui = ui, server = server)
