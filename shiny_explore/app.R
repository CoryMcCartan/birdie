library(shiny)
source(here::here("R/simulation.R"))

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Imputation Bias"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("x_conf", "X confounding",
                        min=-1, max=1, value=0, step=0.1),
            sliderInput("r_bias", "R probability bias",
                        min=-1, max=1, value=0, step=0.1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("ex_plot"),
            textOutput("rmse"),
            plotOutput("sim_plot"),
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$ex_plot <- renderPlot({
        ex_plot(x_conf=input$x_conf, r_bias=input$r_bias)
    })

    res = reactive({
        run_sim(100, x_conf=input$x_conf, r_bias=input$r_bias)
    })

    output$rmse <- renderText({
        r = res()
        rmse = sqrt(mean((r$diff_act - r$diff_est)^2))
        str_glue("RMSE: {percent(rmse, 0.1)}")
    })

    output$sim_plot <- renderPlot({
        r = res()
        p1 = qplot(diff_act - diff_est, data=r, bins=20, xlim=c(-0.15, 0.15))
        p1
    })
}

# Run the application
shinyApp(ui = ui, server = server)
