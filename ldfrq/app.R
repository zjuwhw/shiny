library(shiny)

ui <- fluidPage(
    withMathJax(),
    titlePanel("Allele Frequencies and the Linkage Disequilibrium"),
    sidebarLayout(
        sidebarPanel(
            wellPanel(
                h4("Allele Frequencies"),
                sliderInput(inputId = "pA", label = "P(A)", value = 0.75, min = 0.5, max = 1),
                sliderInput(inputId = "pB", label = "P(B)", value = 0.75, min = 0.5, max = 1)
            ),
            wellPanel(
                h4("Notation for Haplotype and Allele Frequencies"),
                tableOutput("table"),
                span("A and B are the major alleles for locus 1 and 2"),
                h4("Different measurements for LD"),
                p("D=pAB-pA*pB"),
                p("D'=D/min[pA(1-pB),pB(1-pA)] if D>0"),
                p("or D/min[pApB,(1-pA)(1-pB)] if D<0"),
                p("r^2=D^2/pApB(1-pA)(1-pB)")
            ),
            wellPanel(
                p("Author: Huanwei Wang; Email: huanwei.wang@uq.edu.au")
            )
        ),
        mainPanel(
            wellPanel(
                h5("The range of haplotype frequencies and LD measures:"),
                uiOutput("rangepAB"),
                uiOutput("rangeD"),
                uiOutput("rangeDprime"),
                uiOutput("ranger2")
            ),
            plotOutput("plot")
        )
    )
)


server <- function(input, output){
    output$table = renderTable({
        cbind("Locus" = c("A", "a", "Frequency"),
              B = c("P(AB)", "P(aB)", "P(B)"),
              b = c("P(Ab)", "P(ab)", "P(b)"),
              Frequency = c("P(A)", "P(a)", "1"))
    })

    minnum = reactive({input$pA+input$pB-1 })
    maxnum = reactive({min(input$pA,input$pB)})
    output$rangepAB = renderText({paste("p(AB): from ", minnum()," to ", maxnum(), sep="")})
    output$rangeD = renderText({paste("D: from ",minnum()-input$pA*input$pB, " to ", maxnum()-input$pA*input$pB, sep="")})
    output$rangeDprime = renderText({"D': from -1 to 1"})
    output$ranger2 = renderText({paste("r^2: from 0 to ", min(input$pA*(1-input$pB)/input$pB/(1-input$pA),input$pB*(1-input$pA)/input$pA/(1-input$pB)),sep="")})
    
    output$plot = renderPlot({
        pABseq = seq(from=minnum(), to=maxnum(), length.out=100)
        D = pABseq-input$pA*input$pB
        Dprime = ifelse(D>0, D/min(input$pA*(1-input$pB),input$pB*(1-input$pA)), D/min(input$pA*input$pB, (1-input$pA)*(1-input$pB)))
        r2 = D^2/(input$pA*input$pB*(1-input$pA)*(1-input$pB))
        layout(matrix(1:3,nrow=1))
        plot(pABseq, D, type="l", xlab="p(AB)", ylab="D")
        abline(v=input$pA*input$pB,col="red")
        plot(pABseq, Dprime, type="l", xlab="p(AB)",ylab="D'")
        abline(v=input$pA*input$pB,col="red")
        plot(pABseq, r2, type="l", xlab="p(AB)", ylab="r2")
        abline(v=input$pA*input$pB,col="red")
    })
}

shinyApp(ui = ui, server = server)

