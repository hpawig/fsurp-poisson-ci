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
    
    # drop down menu to choose method
    card(
      selectInput(
        inputId = "method",
        label = "Select a Confidence Procedure",
        choices = list("Wald" = 1, "Rao's Score" = 2, "Wilks' Likelihood Ratio" = 3, 
                       "Analog to Clopper-Pearson" = 4, "Modified Stern/Optimal Coverage" = 5,
                       "Crow & Gardner" = 6)
      ),
      
      # slider to choose confidence level
      sliderInput(
        inputId = "conf_level",
        label = "Specify Confidence Level %",
        min = 80, max = 99, value = 95
      ),
      
      # user enters their observed x
      numericInput(inputId = "obs_x", 
                   label = "Enter observed x", value = 10),

      
      # user choice of intervals to be displayed
      checkboxInput(inputId = "checkbox",
                    label = "Display intervals up to observed x",
                    value = FALSE),
      selectInput(
        inputId = "digits",
        label = "Choose up to 4 decimal places",
        choices = list("2" = 1, "3" = 2, "4" = 3)
      ),
      
      # button to submit user's changes
      submitButton("Submit", icon("refresh")),
      
      # contact information section
      card_footer(
        div("Shiny App & Base R code by" , 
            a(href = "https://www.linkedin.com/in/hannahpawig/", target = "_blank", "Hannah Pawig"),
            align = "right", style = "font-size: 8pt"),
        div("Maintenance by" , 
            a(href = "https://www.linkedin.com/in/hannahpawig/", target = "_blank", "Hannah Pawig"),
            align = "right", style = "font-size: 8pt"),
        div("Contact: hpawig@calpoly.edu",
            align = "right", style = "font-size: 8pt"),
        div("Cal Poly, San Luis Obispo",
            align = "right", style = "font-size: 8pt")
      )
    ), 
  ),
  
  # main panel elements
  
  card(
    full_screen = TRUE,
    card_header("Confidence Interval for x"),
    tableOutput(outputId = "intervals")
   )
)



# Define server logic required to output table
 server <- function(input, output) {

  digits <- reactive({
    as.numeric(input$digits) + 1
    })
  


  # output table of confidence interval(s)
    output$intervals <- renderTable({
      expr = find_ci(method = input$method, x = input$obs_x,
                     conf.level = input$conf_level, all = input$checkbox)
    }, digits = digits)
    

}
# Run the application 
shinyApp(ui = ui, server = server)
