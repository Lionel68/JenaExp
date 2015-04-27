shinyUI(fluidPage(
  titlePanel("Exploring FD index"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("S","Species richness",2,20,10,1),
      sliderInput("var","Variation in species abundance",0.1,10,1,0.5),
      sliderInput("var_t1","Variation in T1",0.1,3,1,0.1),
      sliderInput("var_t2","Variation in T2",0.1,10,1,1)),
    mainPanel(
      h2("Some explanations:"),
      p("The points in the graphs represent the different species, the size of the point is proportionnal to their relative abundances"),
      p("The crossed red point is the weighted centroid, the circle in FDiv is the average distance from the vertices to the centroid"),
      p("The last graph represent the slope of the minimum spanning tree, see Ricotta & Moretti 2008 Community Ecology"),
      plotOutput("Plot1"),
      plotOutput("Plot2"),
      plotOutput("Plot3")      
    )
  )
))