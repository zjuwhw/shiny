library(shiny)
ui <- fluidPage(
  sliderInput(inputId = "pA", label = "Frequency of major allele A at locus 1 P(A)", value = 0.75, min = 0.5, max = 1),
  sliderInput(inputId = "pB", label = "Frequency of major allele B at locus 2 P(B)", value = 0.75, min = 0.5, max = 1),
  textOutput("haplotype"),
  plotOutput("Dplot"),
  plotOutput("Dprimeplot"),
  plotOutput("r2plot"),
  sliderInput(inputId = "pAB", label = "Frequency of haplotyp AB P(AB)", value = 0, min = 0, max = 1),
  tableOutput("table1"),
  tableOutput("table2"),
  textOutput("ld")
)

ld = function(pA,pB,pAB){
  D = pAB-pA*pB
  Dprime = ifelse(D>0, D/min(pA*(1-pB),pB*(1-pA)), D/min(pA*pB, (1-pA)*(1-pB)))
  r2 = D^2/(pA*pB*(1-pA)*(1-pB))
  return(list(D,Dprime,r2))
}
server <- function(input, output){

  minnum = reactive({input$pA+input$pB-1}) 
  maxnum = reactive({min(input$pA,input$pB)})
  output$haplotype = renderText({paste("pAB: ", minnum," - ", maxnum, sep="")})
  #output$Dplot = renderPlot({
  #  pABseq = seq(minnum, maxnum, 100)
  #  plot(pABseq,ld(input$pA, input$pB, pABseq),xlab="pAB",y="D")
  #})
#  output$table1 = renderTable({
#    cbind("Locus" = c("A", "a", "Frequency"),
#          B = c("P(AB)", "P(aB)", "P(B)"),
#          b = c("P(Ab)", "P(ab)", "P(b)"),
#          Frequency = c("P(A)", "P(a)", "1"))
#  })
#  output$table2 = renderTable({
#    cbind(Locus = c("A", "a", "Frequency"), 
#          B = c(input$pAB, input$pB-input$pAB, input$pB),
#          b = c(input$pA-input$pAB, 1-input$pA-input$pB+input$pAB, 1-input$pB),
#          Frequency = c(input$pA,1-input$pA,1))
#  })
#  output$note = renderText({
#    ldstat = ld(input$pA, input$pB, input$pAB)
#    paste("The frequency of haplotype with two major allele range from ", input$pA+input$pB-1," to ", min(input$pA, input$pB),"; ",
#          "D = ",ldstat$D, "; ", "D' = ", ldstat$Dprime, "; ", "r2 = ", ldstat$r2, ".", sep="")
#  })
}

shinyApp(ui = ui, server = server)
