library(shiny)
library(bslib)
library(bsicons)
library(rhandsontable)
library(DT)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(plotly)
library(gslnls)
library(purrr)
cbbPalette <- c( '#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB')
theme_set(theme_bw())

#--------------------------------------------------------------------------------------------------------------------
##FUNCTIONS----------------------------------------------------------------------------------------------------------
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
  return(1.81*SiOg+0.31)
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
  plot + geom_function(fun = funct, aes(colour = label), n = 400)
}
#regression function for predicted Rfs
perform_regression <- function(data) {
  gsl_nls(Rf ~ (solvent_percent/phi_nought)^(-1/k), 
          data, 
          start = list(phi_nought = 8, k = -0.5),
          algorithm = "lm",
          control = list(scale = "levenberg"),
          upper = list(k = 0),
          lower = list(phi_nought = 0)) %>%
    coef()
}

#--------------------------------------------------------------------------------------------------------
#UI------------------------------------------------------------------------------------------------------
accordion_filters <- accordion(
  accordion_panel(
    title = "Data Entry", 
    icon = bs_icon("menu-app"),
    selectInput("data_entry", "Rf values:", c("Empirical", "Predicted"), selected = "Empirical" ),
    numericInput("crude_mass", "Crude Mass (g)", 0.5, min = 0),
    numericInput("no_spots", "Number of Spots", value = 3, min = 1, step = 1),
    rHandsontableOutput("Xa"),
    conditionalPanel(
      condition = "input.data_entry == 'Empirical'",
      uiOutput("sliders")
    ),
    conditionalPanel(
      condition = "input.data_entry == 'Predicted'",
      sliderInput("chosen_solvent_percent", "Strong solvent %", 1, 100, 15, animate = animationOptions(interval = 800, loop = TRUE), post="%"),
      span("Spot data", tooltip(bs_icon("info-circle"),
                                p("Input at least 2 Rf values for each spot. You can add/remove datapoints by right clicking and navigating the pop-up menu.")
      ) 
      ),
      rHandsontableOutput("data"),
      textOutput("pred")
    )
  ),
  accordion_panel(
    title = "Optional Settings",
    icon = bs_icon("sliders"),
    card("Column size and fraction size is automatically calculated based off Still's recommendations. If you wish to set your own column size or fraction size, you can do so here. Otherwise, set to 0."),
    numericInput("user_silica", "Silica (g):", value = 0),
    verbatimTextOutput("silica_warning"),
    numericInput("user_fraction_size", "Fraction size (mL):", value = 0, min = 0),
  )
)

ui <- fluidPage(
  title = "Chromatogram Simulator",
  page_sidebar(
    sidebar = sidebar(accordion_filters, width = 350),
    
    #Main
    layout_columns(
      col_widths = c(6, 6, 9, 3, 12, 12),
      fill = FALSE,
      value_box(
        title = "Recommended Silica (g)",
        value = textOutput("rec_sil"),
        showcase = bs_icon("bar-chart"),
        p("Selected silica"),
        textOutput("user_sil")
      ),
      value_box(
        title = "Recommended Fraction Size (mL)",
        value = textOutput("rec_frac"),
        showcase = bs_icon("graph-up"),
        p("Selected fraction size"),
        textOutput("user_frac")
      ),
      card(card_header("Chromatogram",
                       popover(
                         #Plot settings go here!
                         title = "Plot settings",
                         bs_icon("gear"),
                         radioButtons("x_axis_scale", "X axis:", 
                                      c("mL Eluent" = "ml", 
                                        "Fractions" = "frac", 
                                        "Column Volumes" = "cv"), 
                                      selected = "ml"),
                         radioButtons("frac_lines", "Show fraction lines:", 
                                      c("Yes" = "y", "No" = "n"), 
                                      selected = "n")
                       ),
                       class = "d-flex justify-content-between"),
           plotlyOutput(outputId = "chromatogram"),
           max_height = 500
      ),
      
      card( card_header("TLC"), plotOutput(outputId = "tlc"), max_height = 500 ), 
      
      card( card_header("Separation details"), tableOutput("resolution") ),
      #Gutter---------------------------------------------------------------------------------------------
      tags$div("This Shiny app is based off the work of J. D. Fair and C. M. Kormos, Journal of Chromatography A, 2008, 1211, 49–54.",tags$br(),
               "Rf prediction algorithm is based off of P. Kręcisz, K. Czarnecka and P. Szymański, Journal of Chromatographic Science, 2022, 60, 472–477.",
               style = "font-size:10px;")
    )
    
  )
)


#----------------------------------------------------------------------------------------------------------
#SERVER----------------------------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  pred_values <- reactiveValues(data = tibble("spot" = factor(c(1,2,3,1,2,3,1,2,3), levels = 1:3),
                                              "solvent_percent" = c(10, 10, 10, 15, 15, 15, 20, 20, 20),
                                              "Rf" = c(0.54, 0.23, 0.1, 0.64, 0.42, 0.17, 0.7, 0.56, 0.28) ),
                                Xa = c(0.7, 0.1, 0.2),
                                emp_Rf = c(0.54, 0.23, 0.1))
  #update values
  observe({
    req(input$data)
    pred_values$data = hot_to_r(input$data)
  })
  observe({
    req(input$Xa)
    pred_values$Xa = as.vector(hot_to_r(input$Xa))
  })
  
  observe(label = "val_update", priority = 1, {
    req(input$no_spots)
    
    l_max <- max(as.numeric(levels(pred_values$data$spot)))
    spot_diff <- input$no_spots - l_max
    
    if (spot_diff > 0 ) {
      
      #increase levels of spot rows, emp_Rf, and Xa
      levels(pred_values$data$spot) <-  1:input$no_spots
      pred_values$emp_Rf <-  append(pred_values$emp_Rf, rep(0.2, spot_diff))
      pred_values$Xa <-  append(pred_values$Xa, rep(0.5, spot_diff))
      
    } else if (spot_diff < 0 ) {
      #delete rows of larger spot no's, change levels of spot rows, emp_Rf, and Xa
      pred_values$data <- filter(pred_values$data, as.numeric(spot) <= input$no_spots) %>%
        droplevels()
      levels(pred_values$data$spot) <-  1:input$no_spots
      pred_values$emp_Rf <-  pred_values$emp_Rf[1:input$no_spots]
      pred_values$Xa <-  pred_values$Xa[1:input$no_spots]
    }
  })
  
  ##USER TLC DATA INPUT -- ISOCRATIC ---------------------------------------------------------------------------------------------
  
  # Make spot i slider update the value of Rf[i] using map
  observe({
    map(1:input$no_spots, ~{
      observeEvent(input[[paste0("Spot ", .x)]], {
        pred_values$emp_Rf[.x] <- input[[paste0("Spot ", .x)]]
      })
    })
  }) %>% bindEvent(input$no_spots)
  
  output$sliders <- renderUI({ 
    # First, create a list of sliders each with a different name
    sliders <- map(1:input$no_spots, ~{
      inputName <- paste0("Spot ", .x)
      sliderInput(inputName, inputName, min = 0, max = 1, value = pred_values$emp_Rf[.x])
    })
    # Create a tagList of sliders (this is important)
    do.call(tagList, sliders)
  }) %>% bindEvent(input$no_spots)
  
  ##USER TLC DATA INPUT -- RF PREDICTION ---------------------------------------------------------------------------------------- 
  
  # Group by spot and perform regression for each group
  
  predict_data_coef <- reactive({ 
    validate(
      need(any(sapply(levels(pred_values$data$spot), function(x) sum(pred_values$data$spot == x)) < 2) == FALSE, "Error: each spot requires at least two data points"),
      need(any(is.na(pred_values$data)) == FALSE, "Error: please populate all cells")
    )
    pred_values$data %>% group_split(spot) %>% map(perform_regression) 
  }) %>% 
    bindCache(pred_values$data) %>%
    bindEvent(pred_values$data, input$data_entry == "Empirical")
  
  # Apply the function to each element of predict_data_coef
  predicted_rf_values <- reactive({
    predicted_rf_values <- try(sapply(predict_data_coef(), function(x) (input$chosen_solvent_percent / x[[1]])^(-1/x[[2]])) )
    predicted_rf_values[predicted_rf_values > 1] <- 1
    predicted_rf_values[predicted_rf_values < 0] <- 0
    return(predicted_rf_values)
  }) %>% bindCache(predict_data_coef(), input$chosen_solvent_percent)
  
  #outputs & editable tables
  output$pred <- renderText({
    c("Predicted Rf:", 
      round(predicted_rf_values(), 2)
    )
  })
  
  output$Xa <-  renderRHandsontable({
    rhandsontable(as.matrix(pred_values$Xa)%>%`colnames<-`("Mass Balance"))
  })
  
  output$data <-  renderRHandsontable({
    rhandsontable(pred_values$data) 
  })
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #CALCULATIONS ----------------------------------------------------------------------------------------------------
  
  
  rv <- reactiveValues(tlc_data = NULL,
                       is_hard = NULL,
                       silica_mass = NULL,
                       Vm = NULL,
                       plotdata = data.frame("solvent" = rep(0, length = 500)),
                       sep = NULL,
                       comb_plot = NULL,
                       fraction_size = NULL,
                       less_than_still = NULL
  )
  
  #load Rf/Xa data into calculations
  #separate these equations to make it run more efficiently, using bindcache
  observe({
    if (input$data_entry == "Empirical") {
      rv$tlc_data <- tibble("Rf" = pred_values$emp_Rf, "Xa" = pred_values$Xa)
    } else {
      req(predicted_rf_values())
      #Stop app from crashing when number of predicted values doesn't match length of Xa (because we're still inputting the values)
      req(length(predicted_rf_values()) == length(pred_values$Xa) )
      rv$tlc_data <- tibble("Rf" = predicted_rf_values(), "Xa" = pred_values$Xa) 
    }
  }) %>%
    bindEvent({ #make it update only when the table/crude_mass/user_silica is updated
      input$data_entry
      predicted_rf_values()
      pred_values$emp_Rf
      pred_values$Xa
    }, label = "TLC table")
  
  
  #do calculations
  observe({ 
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
      if(rv$is_hard == TRUE){
        rv$less_than_still <- input$user_silica < siog_hard(input$crude_mass)
      }else{
        rv$less_than_still <- input$user_silica < siog_easy(input$crude_mass)
      }
      rv$silica_mass <- input$user_silica
    } else {
      rv$less_than_still <- FALSE
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
      rv$tlc_data
      input$crude_mass
      input$user_silica
      input$user_fraction_size
    }, label = "calculations")
  
  #Resolution Data Output
  output$resolution <- renderTable(rv$sep)
  
  #Silica and Fraction Data output
  output$silica_warning <- renderText({
    if(rv$less_than_still == TRUE){
      c("Warning: Silica is below recommended amount!")
    }else ""
  })
  
  output$user_frac <- renderText(rv$fraction_size)
  output$rec_frac <- renderText(rv$Vm/3)
  
  output$user_sil <- renderText(rv$silica_mass)
  output$rec_sil <- renderText({
    if (rv$is_hard == TRUE) siog_hard(input$crude_mass)
    else siog_easy(input$crude_mass)
  })
  
  
  #-----------------------------------------------------------------------------------------------------------
  #PLOTTING---------------------------------------------------------------------------------------------------
  output$tlc <- renderPlot({
    no_spots <- nrow(rv$tlc_data)
    tlc <- ggplot(data=rv$tlc_data, aes(x = rep(0, no_spots), y = Rf, colour = as.factor(1:no_spots) )) +
      geom_point(size=5) +
      ylim(0,1) +
      scale_x_discrete(labels = NULL, breaks = NULL) +
      labs(x = NULL, y = "R<sub><i>f</i></sub>") +
      theme(legend.position="none", 
            panel.grid.minor=element_blank(),
            axis.title.y = element_markdown(),
            aspect.ratio = 2,
            text = element_text(size = 20))
    tlc + scale_colour_manual(values=cbbPalette) + geom_hline(yintercept = 0) + geom_hline(yintercept = 1)
  }) %>% bindCache(rv$tlc_data$Rf)
  
  
  
  output$chromatogram <- renderPlotly({
    req(rv$tlc_data$bandsize)
    
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
    
    #show fraction lines if selected
    if(input$frac_lines == "y"){
      comb_plot <- comb_plot + geom_vline(xintercept = seq( 0, (rv$Vm * 10), by =  rv$fraction_size), alpha = 0.5, linetype = 3, linewidth = 0.5, colour = "#009988" )
    }
    
    comb_plot <- switch(input$x_axis_scale,
                        ml = comb_plot + labs( x = "Eluent (mL)", y = "Relative peak intensity \n", colour = ""),
                        frac =  comb_plot + labs( x = "Fraction number", y = "Relative peak intensity \n", colour = "") + scale_x_continuous(labels = scales::label_number(scale = 1/rv$fraction_size), breaks = seq( 0, rv$Vm * 10, by =  rv$fraction_size*2) ),
                        cv =  comb_plot + labs( x = "Column Volume", y = "Relative peak intensity \n", colour = "") + scale_x_continuous(labels = scales::label_number(scale = 1/rv$Vm), breaks = seq( 0, rv$Vm * 10, by =  rv$Vm) )
    )
    
    #display correct scale
    ggplotly(comb_plot + 
               theme(text = element_text(size = 14)) +
               scale_colour_manual(values=cbbPalette)
    )
  }) %>% bindCache(rv$tlc_data, input$x_axis_scale, rv$Vm, rv$fraction_size, input$frac_lines)
  
}

# Call the shinyapp
shinyApp(ui = ui, server = server)
