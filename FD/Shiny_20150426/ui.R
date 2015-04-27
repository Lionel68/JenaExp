shinyUI(fluidPage(
  titlePanel("Exploring FD index"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("S","Species richness",2,20,10,1),
      sliderInput("var","Variation in species abundance",0.1,10,1,0.5),
      sliderInput("var_t1","Variation in T1",0.1,3,1,0.1),
      sliderInput("var_t2","Variation in T2",0.1,10,1,1)),
    mainPanel(
      plotOutput("Plot1"),
      plotOutput("Plot2"),
      plotOutput("Plot3"),
      plotOutput("Plot4")
    )
  )
))