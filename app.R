#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("stage_simulation.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Disease Spread Simulation with Flexible/Dynamic Interventions"),
  
  # Sidebar with input choices
  sidebarLayout(
    sidebarPanel(
      actionButton("go1", "Simulate!",
                   icon("refresh"), 
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
      numericInput(
        "N",
        "Population size:",
        min = 20,
        max = 500,
        value = 100,
        step = 1
      ),
      sliderInput(
        "p0",
        "Initial network edge density:",
        min = 0.001,
        max = 0.999,
        value = 0.05
      ),
      numericInput(
        "R0",
        "Basic reproduction number (R0):",
        value = 4,
        min = 1.5,
        max = 16,
        step = 0.5
      ),
      numericInput(
        "meanRecov",
        "Average recovery time (in days):",
        value = 7,
        min = 2,
        max = 30,
        step = 1
      ),
      selectInput(
        "normalContact",
        "Contact pattern without intervention:",
        c(
          "Static social contacts" = "static",
          "Dynamic social contacts" = "dynamic"
        )
      ),
      selectInput(
        "NPI",
        "Social intervention strategies:",
        list(
          `Population wide` = list("do nothing", "quarantine"),
          `Case isolation` = list("case isolation only",
                                  "plus population quarantine")
        )
      ),
      selectInput(
        "stageType",
        "Phase social interventions by:",
        c(
          "Days since 1st case" = "time",
          "Cumulative number of cases" = "cases"
        )
      ),
      numericInput(
        "st1",
        "Intervention starting point:",
        value = 5,
        min = 0,
        max = 500,
        step = 1
      ),
      numericInput(
        "st2",
        "Intervention end point:",
        value = 30,
        min = 0,
        max = 500,
        step = 1
      )#,
      # actionButton("go2", "Simulate!",
      #              icon("refresh"), 
      #              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
    ),
    
    # Show plots
    mainPanel(
      tabsetPanel(
        tabPanel("Prevalence curve", 
                 plotOutput("prevalPlot"),
                 actionButton("go2", "Simulate!",
                              icon("refresh"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        ),
        tabPanel("Reproduction number", 
                 plotOutput("RtPlot"),
                 actionButton("go3", "Simulate!",
                              icon("refresh"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 ),
        tabPanel("Network density",
                 plotOutput("densPlot"),
                 actionButton("go4", "Simulate!",
                              icon("refresh"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 )
      )
      # plotOutput("prevalPlot"),
      # plotOutput("RtPlot"),
      # plotOutput("densPlot")
    )
  )
)


# a function to generate link rates for use
gen_LR <- function(p0, no_intervention, intervention){
  if(no_intervention == "static"){
    default_alpha = default_omega = 0
    qua_omega = .05
  }else{
    default_omega = qua_omega = .05
    default_alpha = default_omega * p0
  }
  if(intervention == "do nothing"){
    LR = list(
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3)),
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3)),
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3))
    )
  }else if(intervention == "quarantine"){
    LR = list(
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3)),
      list(alpha.r = rep(default_alpha*.5, 3), alpha.d = rep(qua_omega, 3)),
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3))
    )
  }else if(intervention == "case isolation only"){
    LR = list(
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3)),
      list(alpha.r = c(default_alpha, 0, 0), 
           alpha.d = c(default_omega, .2, .2)),
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3))
    )
  }else{
    LR = list(
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3)),
      list(alpha.r = c(default_alpha*.5, 0, 0), 
           alpha.d = c(qua_omega, .2, .2)),
      list(alpha.r = rep(default_alpha, 3), alpha.d = rep(default_omega, 3))
    )
  }
  return(LR)
}

# pre-generate a simulation result
set.seed(42)
res = stage_coevolve(N = 100, bet = 4 / 7, gam = 1 / 7,
                     init.p = .05, stage_type = 'time',
                     stage_starts = c(0, 5, 30),
                     link_rates = gen_LR(.05, "static", "do nothing"),
                     plot = F)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  rv = reactiveValues(res = res)
  
  observeEvent({
    as.numeric(input$go1) + as.numeric(input$go2) + 
      as.numeric(input$go3) + as.numeric(input$go4)
  },{
    rv$res = stage_coevolve(
      N = input$N,
      bet = input$R0 / input$meanRecov,
      gam = 1 / input$meanRecov,
      init.p = input$p0,
      stage_type = input$stageType,
      stage_starts = c(0, input$st1, input$st2),
      link_rates = gen_LR(input$p0, input$normalContact, input$NPI),
      plot = F
    )
    while(sum(rv$res$events$event == 1) < 0.05 * input$N){
      rv$res = stage_coevolve(
        N = input$N,
        bet = input$R0 / input$meanRecov,
        gam = 1 / input$meanRecov,
        init.p = input$p0,
        stage_type = input$stageType,
        stage_starts = c(0, input$st1, input$st2),
        link_rates = gen_LR(input$p0, input$normalContact, input$NPI),
        plot = F
      )
    }
  }, ignoreInit = T)
  
  # res <- eventReactive({
  #   as.numeric(input$go1) + as.numeric(input$go2)
  #   }, {
  #   stage_coevolve(
  #     N = input$N,
  #     bet = input$R0 / input$meanRecov,
  #     gam = 1 / input$meanRecov,
  #     init.p = input$p0,
  #     stage_type = input$stageType,
  #     stage_starts = c(0, input$st1, input$st2),
  #     link_rates = gen_LR(input$p0, input$normalContact, input$NPI),
  #     plot = F
  #   )
  # }, ignoreNULL = F, ignoreInit = F)
  
  # res <- eventReactive(input$go2, {
  #   stage_coevolve(
  #     N = input$N,
  #     bet = input$R0 / input$meanRecov,
  #     gam = 1 / input$meanRecov,
  #     init.p = input$p0,
  #     stage_type = input$stageType,
  #     stage_starts = c(0, input$st1, input$st2),
  #     link_rates = gen_LR(input$p0, input$normalContact, input$NPI),
  #     plot = F
  #   )
  # }, ignoreNULL = T, ignoreInit = T)
  
  # event.log <- reactive(res()$events)
  # stage_times <- reactive(res()$stage_times)
  # total_infec <- isolate(sum(event.log()$event == 1))
  # pop_size <- isolate(input$N)
  # 
  # while(total_infec < 0.05 * pop_size){
  #   res <- reactive(
  #     stage_coevolve(
  #       N = input$N,
  #       bet = input$R0 / input$meanRecov,
  #       gam = 1 / input$meanRecov,
  #       init.p = input$p0,
  #       stage_type = input$stageType,
  #       stage_starts = stage_starts(),
  #       link_rates = LR(),
  #       plot = F
  #     )
  #   )
  #   event.log <- reactive(res()$events)
  #   stage_times <- reactive(res()$stage_times)
  #   total_infec <- isolate(sum(event.log()$event == 1))
  #   pop_size <- isolate(input$N)
  # }
  
  output$prevalPlot <- renderPlot({
    plot(
      preval ~ time,
      data = rv$res$events,
      xlab = "Days since 1st case",
      ylab = "Fraction of active (infectious) cases",
      main = "Disease prevalence vs. Days",
      ylim = c(0, 1),
      type = "l",
      lwd = 2
    )
    abline(v = rv$res$stage_times, col = "red", lwd = 2)
  })
  
  output$RtPlot <- renderPlot({
    plot(
      Rt ~ time,
      data = rv$res$events,
      xlab = "Days since 1st case",
      ylab = "R_t",
      main = "Reproduction number vs. Days",
      type = "l",
      lwd = 2
    )
    abline(v = rv$res$stage_times, col = "red", lwd = 2)
  })
  
  output$densPlot <- renderPlot({
    plot(
      dens ~ time,
      data = rv$res$events,
      xlab = "Days since 1st case",
      ylab = "Edge density",
      main = "Network edge density vs. Days",
      type = "l",
      lwd = 2
    )
    abline(v = rv$res$stage_times, col = "red", lwd = 2)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
