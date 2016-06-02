shinyUI(fluidPage(
  titlePanel("Dynamics of early forest succession in a 1-ha patch"),
  
    sidebarLayout(
    sidebarPanel(
      sliderInput("MAXY",
                  label=div(HTML("Annual seed arrival (ha<sup>-1</sup>)")),
                  min = 0.01, max = 5000, value = 200,step=1),
          sliderInput("GROWTHY",
                  "Annual height growth of trees in canopy (m)",
                  min = 0.01,
                  max = 3,
                  value = 0.36),
      sliderInput("SURVY",
                  "Annual survival of trees in canopy",
                  min = 0.001,
                  max = 0.99,
                  value = 0.86),
      sliderInput("SURVYg",
                  "Annual survival of trees in grass layer",
                  min = 0.001,
                  max = 0.99,
                  value = 0.57),
      sliderInput("GROWTHYg",
                  "Annual height growth of trees in grass layer (m)",
                  min = 0.01,
                  max = 3,
                  value = 0.36)
          ),
    
    mainPanel(
    plotOutput("CANOPYplot"),
    textOutput("text1")
    )
  )
))

#
