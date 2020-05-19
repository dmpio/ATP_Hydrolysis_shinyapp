library(shiny)
library(tibble)
library(ggplot2)
library(dplyr)
library(magrittr)

source(file = "constants.R")
source(file = "functions.R")

### Reactive Function to iterate over a variable

iterVar <- function(selected_var, tempC, buffer_IS, pH, mg, adp, atp, phosphate, seq_divisions){
  
  iter_var <- switch(selected_var, 
                     "tempC" = tempC,
                     "buffer_IS" = buffer_IS, 
                     "pH" = pH, 
                     "mg" = mg, 
                     "atp" = atp, 
                     "adp" = adp, 
                     "phosphate" = phosphate
                     )
  
  iter_range <-  seq(from = iter_var[[1]],
                     to = iter_var[[2]], 
                     length.out = seq_divisions)
  
  to_plot = list()
  
  # iterate over the range of numbers
  for (x in iter_range){
    # Assign the input variables based on the plotted variable (x)
    iter_template <- switch(selected_var, 
                            "tempC" = list(tempC = x, 
                                           buffer_IS = buffer_IS,
                                           pH = pH,
                                           mg = mg,
                                           atp = atp,
                                           adp = adp,
                                           phosphate = phosphate
                            ),
                            
                            "buffer_IS" = list(tempC = tempC, 
                                               buffer_IS = x,
                                               pH = pH,
                                               mg = mg,
                                               atp = atp,
                                               adp = adp,
                                               phosphate = phosphate
                            ),
                            
                            "pH" = list(tempC = tempC, 
                                        buffer_IS = buffer_IS,
                                        pH = x,
                                        mg = mg,
                                        atp = atp,
                                        adp = adp,
                                        phosphate = phosphate
                            ),
                            
                            "mg" = list(tempC = tempC, 
                                        buffer_IS = buffer_IS,
                                        pH = pH,
                                        mg = x,
                                        atp = atp,
                                        adp = adp,
                                        phosphate = phosphate
                            ),
                            
                            "atp" = list(tempC = tempC, 
                                         buffer_IS = buffer_IS,
                                         pH = pH,
                                         mg = mg,
                                         atp = x,
                                         adp = adp,
                                         phosphate = phosphate
                            ),
                            
                            "adp" = list(tempC = tempC, 
                                         buffer_IS = buffer_IS,
                                         pH = pH,
                                         mg = mg,
                                         atp = atp,
                                         adp = x,
                                         phosphate = phosphate
                            ),
                            
                            "phosphate" = list(tempC = tempC, 
                                               buffer_IS = buffer_IS,
                                               pH = pH,
                                               mg = mg,
                                               atp = atp,
                                               adp = adp,
                                               phosphate = x
                            ),
    )
    
    # get the results of the calculation
    results <- calc_apparent_Keq(tempC = iter_template$tempC, 
                                 buffer_IS = iter_template$buffer_IS,
                                 pH = iter_template$pH,
                                 mg = iter_template$mg,
                                 atp = iter_template$atp,
                                 adp = iter_template$adp,
                                 phosphate = iter_template$phosphate,
                                 reaction = "ATP", 
                                 is0_Kref=is0_Kref_atpHyd)
    # assemble the results in a list
    # only selecting the apparent Gibbs energy
    to_plot <- c(to_plot, results$dG_ATP)
  }
  # combine the results list with the variable range for plotting as a df
  d <- as_tibble(list('X' = iter_range, 
                      'j' = unlist(to_plot, recursive=FALSE)))
  
  # calculate different unit options
  # j = Joules/mol, k = Kilojoules/mol, c = kilocals/mol
  d %<>% mutate(k = j/1000, cal = j*0.000239006)
  
  return(d)
}

# server function
function(input, output, session) {
  # Dynamic UI Generation
  # Enable Switching plotted X variable 
  # By changing widget from value to range
  output$tempC_output <- renderUI({
    if (input$x_var == "tempC"){
      val <- ui_tempC_const$default_range
    } else {
      val <- ui_tempC_const$default_value
    }
    sliderInput(inputId = "tempC",
                label = h3("Temperature [°C]"), 
                min =  ui_tempC_const$min,
                max =  ui_tempC_const$max,
                step = ui_tempC_const$step,
                value = val, 
    )
  })
  
  output$bufferIS_output <- renderUI({
    if (input$x_var == "buffer_IS"){
      val <- ui_bufferIS_const$default_range
    } else {
      val <- ui_bufferIS_const$default_value
    }
    sliderInput(inputId = "buffer_IS",
                label = h3("Ionic Strength [mM]"), 
                min =  ui_bufferIS_const$min,
                max =  ui_bufferIS_const$max,
                step = ui_bufferIS_const$step,
                value = val
    )
  })
  
  output$pH_output <- renderUI({
    if (input$x_var == "pH"){
      val <- ui_pH_const$default_range
    } else {
      val <- ui_pH_const$default_value
    }
    sliderInput(inputId = "pH",
                label = h3("pH"), 
                min = ui_pH_const$min,
                max = ui_pH_const$max,
                step = ui_pH_const$step,
                value = val
    )
  })
  
  output$mg_output <- renderUI({
    if (input$x_var == "mg"){
      val <- ui_mg_const$default_range
    } else {
      val <- ui_mg_const$default_value
    }
    sliderInput(inputId = "mg",
                label = h3("Free Mg", tags$sup("2+"), "[mM]"), 
                min = ui_mg_const$min,
                max = ui_mg_const$max,
                step = ui_mg_const$step, 
                value = val
    )
  })
  
  output$atp_output <- renderUI({
    if (input$x_var == "atp"){
      val <- ui_atp_const$default_range
    } else {
      val <- ui_atp_const$default_value
    }
    sliderInput(inputId = "atp",
                label = h3("ATP [mM]"), 
                min = ui_atp_const$min,
                max = ui_atp_const$max,
                step = ui_atp_const$step, 
                value = val
    )
  })
  
  output$adp_output <- renderUI({
    if (input$x_var == "adp"){
      val <- ui_adp_const$default_range
    } else {
      val <- ui_adp_const$default_value
    }
    sliderInput(inputId = "adp",
                label = h3("ADP [mM]"), 
                min = ui_adp_const$min,
                max = ui_adp_const$max,
                step = ui_adp_const$step, 
                value = val
    )
  })
  
  output$phosphate_output <- renderUI({
    if (input$x_var == "phosphate"){
      val <- ui_phosphate_const$default_range
    } else {
      val <- ui_phosphate_const$default_value
    }
    sliderInput(inputId = "phosphate",
                label = h3("Phosphate [mM]"), 
                min = ui_phosphate_const$min,
                max = ui_phosphate_const$max,
                step = ui_phosphate_const$step, 
                value = val
    )
  })
  
  # reactively cally the function to calculate the values
  dataInput <- reactive({
    
    iterVar(selected_var = input$x_var,
            temp = input$tempC, 
            buffer_IS = input$buffer_IS, 
            pH = input$pH, 
            mg = input$mg, 
            adp = input$adp, 
            atp = input$atp, 
            phosphate = input$phosphate,
            seq_divisions = input$num_divisions
            )
  })
  
  
  # plot the calculated values
  output$plot <- renderPlot({
    
    # to prevent the plot from attempting to generate before
    # the rest of the ui loads up
    new_var <- switch(input$x_var, 
                      "tempC" = input$tempC,
                      "buffer_IS" = input$buffer_IS, 
                      "pH" = input$pH, 
                      "mg" = input$mg, 
                      "atp" = input$atp, 
                      "adp" = input$adp, 
                      "phosphate" = input$phosphate)
    if (is.null(new_var) | length(new_var)==1) 
      return()
    
    
    d <- dataInput()
    
    x_label <- switch(input$x_var, 
                      "tempC" = "Temperature (°C)",
                      "buffer_IS" = "Ionic Strength (mM)", 
                      "pH" = "pH", 
                      "mg" = "Free Magnesium (mM)", 
                      "atp" = "ATP (mM)", 
                      "adp" = "ADP (mM)", 
                      "phosphate" = "Phosphate (mM)")
    
    x_lim <- switch(input$x_var, 
                    "tempC" = c(ui_tempC_const$min, ui_tempC_const$max),
                    "buffer_IS" = c(ui_bufferIS_const$min, ui_bufferIS_const$max),
                    "pH" = c(ui_pH_const$min, ui_pH_const$max),
                    "mg" = c(ui_mg_const$min, ui_mg_const$max), 
                    "atp" = c(ui_atp_const$min, ui_atp_const$max),
                    "adp" = c(ui_adp_const$min, ui_adp_const$max),
                    "phosphate" = c(ui_phosphate_const$min, ui_phosphate_const$max)
                    )
    
    # account for the units
    anchors <- switch(input$units, 
                      "j" = list(ymin = -63000, 
                                 ymax = -58000, 
                                 breakby=2000, 
                                 ylabel_unit = "(J/mol)"),
                      "k" = list(ymin = -63, 
                                 ymax = -58, 
                                 breakby=2, 
                                 ylabel_unit = "(kJ/mol)"),
                      "cal" = list(ymin = -15,
                                   ymax = -13.8, 
                                   breakby=0.5, 
                                   ylabel_unit = "(kCal/mol)")
                      )
    # Anchor ymin
    if (min(d[[input$units]]) < anchors$ymin){
      ymin  <-  round(min(d[[input$units]]) - anchors$breakby)
      
    } else {
      ymin <- anchors$ymin
    }
    
    # Anchor ymax
    if (max(d[[input$units]])> anchors$ymax){
      ymax <- round(max(d[[input$units]]) + anchors$breakby)
    } else {
      ymax <- anchors$ymax
    }
    
    p <- ggplot(d, aes(X, .data[[input$units]])) + 
      geom_line(size=1, color="#3F3F3F") + 
      geom_point(color="#3F3F3F") + 
      scale_x_continuous(name = x_label, limits = x_lim) +
      scale_y_continuous(name = paste("Gibbs Energy", anchors$ylabel_unit),
                         limits = c(ymin, ymax),
                         breaks = seq(from=ymin,
                                      to=ymax,
                                      by=anchors$breakby))
    print(p)
  }, 
  res=120)
}