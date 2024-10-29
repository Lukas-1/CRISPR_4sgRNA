### 2024-09-22



# Load packages -----------------------------------------------------------

library("shiny")




# Define UI ---------------------------------------------------------------

ui <- fluidPage(
  selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
  verbatimTextOutput("summary"),
  reactable::reactableOutput("table")
)





# Define server logic -----------------------------------------------------


server <- function(input, output, session) {
  output$summary <- renderPrint({
    dataset <- get(input$dataset, "package:datasets")
    summary(dataset)
  })

  output$table <- reactable::renderReactable({
    dataset <- get(input$dataset, "package:datasets")
    reactable::reactable(dataset)
  })
}

shinyApp(ui, server)



