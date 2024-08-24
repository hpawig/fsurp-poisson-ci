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
options(warn=-1) # turn warnings off temporarily
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
                       "Crow & Gardner" = 6, "Blaker" = 7, "Conditional Minimal Cardinality" = 8)
      ),
      
      # slider to choose confidence level
      sliderInput(
        inputId = "conf_level",
        label = "Specify Confidence Level %",
        min = 80, max = 99, value = 95
      ),
      
      # user enters their observed x
      numericInput(inputId = "obs_x", 
                   label = "Enter observed x", value = 10,
                   min = 0),
      
      card_footer(
        div(tags$small("*Note: " , 
                       tags$i("x"), " must be an integer representing the number of ",
                       tags$i("total"), " events over time period of interest"),
            align = "center", style = "color: red")
      )
      
      ),
    

      card(
      card_header("Options"),
      
      # user choice of intervals to be displayed
      checkboxInput(inputId = "displayAll",
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
        div("Shiny App by" , 
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
  card (
    height = 160,
    card_header("Inputs"),
    verbatimTextOutput(outputId = "inputs")
  ),
  
  card(
    full_screen = T,
    card_header("Confidence Intervals"),
    max_height = 440,
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
                     conf.level = input$conf_level, all = input$displayAll,
                     digits = input$digits) |> select(x:interval)
    }, digits = digits)
    

  # output string of user's inputs
    output$inputs <- renderText({
      find_ci(method = input$method, x = input$obs_x,
              conf.level = input$conf_level, all = input$displayAll,
              digits = input$digits)$input.str[1]
    })

 }
 options(warn=0) # turn warnings back on
# Run the application 
shinyApp(ui = ui, server = server)


