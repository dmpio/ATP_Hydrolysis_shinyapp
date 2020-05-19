library(shiny)


fluidPage(
  
  # Application title
  titlePanel("Apparent Gibbs Energy (ΔG') of ATP Hydrolysis"),
  
  # sidebar with variables
  sidebarLayout(
    sidebarPanel(
      uiOutput("tempC_output"),
      uiOutput("bufferIS_output"),
      uiOutput("pH_output"),
      uiOutput("mg_output"),
      uiOutput('atp_output'),
      uiOutput('adp_output'),
      uiOutput('phosphate_output')
      
    ),
    
    # main panel with plot, variable selection, and number of plotted points
    mainPanel(
      selectInput(
        inputId = 'x_var',
        label = h3("X-axis Variable"),
        choices = list("Temperature" = 'tempC', 
                       "Ionic Strength" = 'buffer_IS', 
                       "pH" = "pH",
                       "Free Mg" = 'mg', 
                       "ATP" = 'atp',
                       "ADP" = 'adp', 
                       "Phosphate" = 'phosphate')
      ),
      selectInput(
        inputId = 'units',
        label = h3("Y-axis Units"),
        choices = list("Joules per mole" = 'j', 
                       "Kilojoules per mole" = 'k', 
                       "Kilocalories per mole" = "cal"
                       ), 
        selected = 'k'
      ),
      numericInput(
        inputId = "num_divisions", 
        label = "# Points", 
        value = 25, 
        min = 3, 
        max=100, 
        step = 1
      ),
      plotOutput(
        outputId = 'plot'
      )
    )
  )
)