#################################################################
##                             App                             ##
#################################################################
library(shiny)
library(bslib)
library(tidyverse)





# Large Sample Procedures (W, RS, Wilks LR)
source("large-sample-methods.R", encoding = "UTF-8")

# Strict Method Procedures (Clopper-Pearson, MST/OC)
source("strict-methods.R", encoding = "UTF-8")

# Server Functions
source("server-fns.R", encoding = "UTF-8")



# Define UI for application that takes in user's observed x and confidence level

ui <- page_sidebar(
  title = "Poisson Confidence Interval Generator",
  sidebar = sidebar(
    width = 450,
    selectInput(
      inputId = "method",
      label = "Select a Confidence Procedure",
      choices = list("Wald" = 1, "Rao's Score" = 2, "Wilks' Likelihood Ratio" = 3, 
                     "Analog to Clopper-Pearson" = 4, "Modified Stern/Optimal Coverage" = 5,
                     "Crow & Gardner" = 6)
    ),
    card(
      card_header("Confidence Level (%)"),
      sliderInput(
        inputId = "conf_level",
        label = "Set value",
        min = 80,
        max = 99,
        value = 95
      )
    ),
    card(
      card_header(withMathJax("Observed x")),
      numericInput(inputId = "obs_x", 
                   label = "Enter a number", value = 1)
    ),
    card(
      checkboxInput(inputId = "checkbox",
                    label = "Display intervals up to observed x",
                    value = FALSE), # width = NULL),
      selectInput(
        inputId = "digits",
        label = "Choose up to 4 decimal places",
        choices = list("2" = 1, "3" = 2, "4" = 3)
      )
    )
  ),
  
  # elements in the main panel
  
  card(
    full_screen = TRUE,
    card_header("Confidence Interval for x"),
    tableOutput(outputId = "intervals")
   ) #,
  # card(
  #   card_header("CPF Graph"),
  #   plotOutput("CPF")
  # )
)



# Define server logic required to draw a histogram
 server <- function(input, output) {

  digits <- reactive({
    as.numeric(input$digits) + 1
    })

  output$ci_header <- renderText({
    "blegh"
  })
  
  # output table of confidence interval(s)
  output$intervals <- renderTable({
    expr = find_ci(method = input$method,x = input$obs_x, 
                   conf.level = input$conf_level, all = input$checkbox) 
  }, digits = digits)
  
  
  # CPF graph
    output$CPF <- renderPlot({
      
    })

}
# Run the application 
shinyApp(ui = ui, server = server)
