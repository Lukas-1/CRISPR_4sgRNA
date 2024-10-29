### 2024-09-22


# Load packages -----------------------------------------------------------

library("shiny")





# Define UI ---------------------------------------------------------------

ui <- fluidPage(
  selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
  verbatimTextOutput("summary"),
  DT::DTOutput("table")
)





# Define server logic -----------------------------------------------------


server <- function(input, output, session) {
  output$summary <- renderPrint({
    dataset <- get(input$dataset, "package:datasets")
    summary(dataset)
  })

  output$table <- DT::renderDT({
    dataset <- get(input$dataset, "package:datasets")
    dataset
  })
}

shinyApp(ui, server)



