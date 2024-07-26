#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
library(tidyverse)

# Define UI for application that draws a histogram
# ui <- fluidPage(
# 
#     # Application title
#     titlePanel("Poisson Confidence Interval Generator"),
# 
#     # Sidebar with a slider input for number of bins
#     sidebarLayout(
#         sidebarPanel(
#             sliderInput(inputId = "conf.level",
#                         label = "Confidence Level (%)",
#                         min = 80,
#                         max = 100,
#                         value = 95),
#             selectInput(inputId = "method",
#                         label = "Choose a confidence procedure:",
#                         choices = list("Wald" = 1, "Rao's Score" = 2,
#                                        "Wilks' Likelihood Ratio" = 3,
#                                        "Clopper-Pearson" = 4,
#                                        "Modified Stern/Optimal Coverage" = 5),
#                         )
#         ),
# 
# 
#         # Show a plot of the generated distribution
#         mainPanel(
#            plotOutput("distPlot")
#         )
#     )
# )

ui <- page_fluid(
  titlePanel("Poisson Confidence Interval Generator"),
  layout_columns(
    card(
      card_header("Select a Confidence Procedure"),
      selectInput(
        "select",
        "Select option",
        choices = list("Wald" = 1, "Rao's Score" = 2, "Wilks' Likelihood Ratio" = 3, 
                       "Clopper-Pearson" = 4, "Modified Stern/Optimal Coverage" = 5),
        selected = 1
      )
    ),
    card(
      card_header("Confidence Level %:"),
      sliderInput(
        inputId = "conf_level",
        label = "Set value",
        min = 80,
        max = 99,
        value = 95
      )
    )
  )
)


ui <- page_sidebar(
  title = "Poisson Confidence Interval Generator",
  sidebar = sidebar(
    width = 450,
    selectInput(
      inputId = "method",
      label = "Select a Confidence Procedure",
      choices = list("Wald" = 1, "Rao's Score" = 2, "Wilks' Likelihood Ratio" = 3, 
                     "Clopper-Pearson" = 4, "Modified Stern/Optimal Coverage" = 5)
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
      numericInput(inputId = "num", 
                   label = "Enter a number", value = 1)
    )
  ),
  
  # elements in the main panel
  
  card(
    full_screen = TRUE,
    card_header("Confidence Intervals go here") #,
    #plotOutput("p")
  ),
  card(
    card_header("CPF Graph"),
    plotOutput("CPF")
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  temp <- ""
  case_when(
    (input$method == "Wald") ~ temp <- "W",
    (input$method == "Rao's Score") ~ temp <- "S",
    (input$method == "Wilks' Likelihood Ratio") ~ temp <- "LR",
    (input$method == "Clopper-Pearson") ~ temp <- "CP",
    (input$method == "Modified Stern/Optimal Coverage") ~ temp <- "OC"
  )
  
  # output table of confidence intervals
  # 
  
  
  # CPF graph
    output$CPF <- renderPlot({
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
