library(shiny)

ui <- fluidPage(
    titlePanel("Allele Frequencies and the Linkage Disequilibrium"),
    sidebarLayout(
        sidebarPanel(
            sliderInput(inputId = "pA", label = "P(A)", value = 0.75, min = 0.5, max = 1),
            sliderInput(inputId = "pB", label = "P(B)", value = 0.75, min = 0.5, max = 1),
            sliderInput(inputId = "pAB", label = "P(AB)", value = 0, min = 0, max = 1)),
        mainPanel(
            plotOutput("plot"),
            textOutput("haplotype"),
            tableOutput("table1"),
            tableOutput("table2"),
            textOutput("ld"))
    )
)

server <- function(input, output){
    minnum = reactive({input$pA+input$pB-1 })
    maxnum = reactive({min(input$pA,input$pB)})
    output$haplotype = renderText({
        paste("p(AB): ", minnum()," - ", maxnum(), sep="")
    })
    output$plot = renderPlot({
        pABseq = seq(from=minnum(), to=maxnum(), length.out=100)
        D = pABseq-input$pA*input$pB
        Dprime = ifelse(D>0, D/min(input$pA*(1-input$pB),input$pB*(1-input$pA)), D/min(input$pA*input$pB, (1-input$pA)*(1-input$pB)))
        r2 = D^2/(input$pA*input$pB*(1-input$pA)*(1-input$pB))
        layout(matrix(1:3,nrow=1))
        plot(pABseq, D, type="l", xlab="p(AB)", ylab="D")
        plot(pABseq, Dprime, type="l", xlab="p(AB)",ylab="D'")
        plot(pABseq, r2, type="l", xlab="p(AB)", ylab="r2")
    })
    output$table1 = renderTable({
        cbind("Locus" = c("A", "a", "Frequency"),
              B = c("P(AB)", "P(aB)", "P(B)"),
              b = c("P(Ab)", "P(ab)", "P(b)"),
              Frequency = c("P(A)", "P(a)", "1"))
    })
    output$table2 = renderTable({
        cbind(Locus = c("A", "a", "Frequency"), 
              B = c(input$pAB, input$pB-input$pAB, input$pB),
              b = c(input$pA-input$pAB, 1-input$pA-input$pB+input$pAB, 1-input$pB),
              Frequency = c(input$pA,1-input$pA,1))
    })
    output$ld = renderText({
        pA=input$pA
        pB=input$pB
        pAB=input$pAB
        D = pAB-pA*pB
        Dprime = ifelse(D>0, D/min(pA*(1-pB),pB*(1-pA)), D/min(pA*pB, (1-pA)*(1-pB)))
        r2 = D^2/(pA*pB*(1-pA)*(1-pB))
        paste("D = ",D, "; ", "D' = ", Dprime, "; ", "r2 = ", r2, ".", sep="")
    })
}

shinyApp(ui = ui, server = server)

