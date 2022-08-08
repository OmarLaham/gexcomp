library(shiny)
ui <- fluidPage(
    titlePanel("CountPlots"),
    sidebarLayout(
        sidebarPanel(textInput(inputId = "Gene", label = "Enter a gene name:", placeholder = "Gene name in uppercase letters"), width = 2),
    mainPanel(
      tabsetPanel(
        tabPanel(title = "CountPlot", plotOutput(outputId = "CountPlot", width = "700px", height="600px")),
        tabPanel(title = "TPM CountPlot", plotOutput(outputId = "TPMPlot", width = "700px", height="600px"))
        
      ),
      p("*** = padj <= 0.001"),
      p(" ** = padj <= 0.01"),
      p("  * = padj <= 0.05"),
      p(" ns = padj > 0.05"),
      p(" NA = removed from statistical analysis")
    )
  )
)
                
   

dds <- readRDS("002_dds.rds")
res <- readRDS("002_res.rds")
TPM <- readRDS("002_TPM_counts.rds")
label_pvs <- function(Gene){
    rn <- which(mcols(dds)$gene_symbol==Gene)
    pv <- res[rn,]$padj
    
    if(is.na(pv)){
        l <- "NA"
    } else if(pv<=0.001){
        l <- "***"
    } else if(pv<=0.01){
        l <- "**"
    } else if(pv<=0.05){
        l <- "*"
    } else {
        l <- "ns"
    }
    return(l)
}             
                
#choose_ref <- function(Ref){
#     if(Ref==1){
#         l2 <- l[c(1:3)]
#     } else if(Ref==2){
#         l2 <- l[c(4:6)]
#     } else if(Ref==3){
#         l2 <- l[c(7:9)]
#     } else if(Ref==4){
#         l2 <- l[c(10:12)]
#     } else {
#         print("Error in reference selection.")
#     }
#     return(l2)
# }
df <- data.frame("Condition"=dds@colData@listData$Condition)


server <- function(input,output) {
    library(DESeq2)
    library(ggplot2)
    set.seed(42)
    #dds@colData@listData$batch <- c("batch0","batch3","batch3","batch3","batch0","batch0","batch1","batch1","batch1","batch2","batch2","batch2")
    l <- reactive({
        label_pvs(Gene=input$Gene)
    })
    #l2 <- reactive({
    #    choose_ref(Ref=input$Ref)
    #})
    p0 <- reactive({plotCounts(dds, gene=rownames(dds)[which(mcols(dds)$gene_symbol==input$Gene)], intgroup = "Condition", returnData = T)})
    p1 <- reactive({
            ggplot(p0(), aes(x= Condition, y=count, col=dds@colData@listData$Cell_line, shape=dds@colData@listData$Batch))+ 
                geom_jitter(size=3, height=0, width = 0.1)+ ggtitle(paste0("Counts of ", input$Gene))+ scale_y_log10()+ 
                scale_x_discrete(limits=c("slow", "fast")) + theme_gray(base_size = 22)+theme(legend.title=element_blank()) +
                annotate("text",size=8, x=2, y=max(p0()$count)+0.1*max(p0()$count), label=paste(l()))+
        ylab("Normalized counts (DESeq2)")+
        xlab("")
    })
    output$CountPlot <- renderPlot({
        (p1())
        }) 
    
    df <- reactive({
      mydf <- data.frame("Condition"=dds@colData@listData$Condition)
      mydf$TPM <- t(TPM[TPM$gene_name==input$Gene,1:12])[,1]
      mydf
    })

    
    p2 <- reactive({
      ggplot(df(), aes(x=Condition, y=TPM, col=dds@colData@listData$Cell_line, shape=dds@colData@listData$Batch))+ 
        geom_jitter(size=3, height=0, width = 0.1)+ ggtitle(paste0("Counts of ", input$Gene))+ scale_y_log10()+ 
        scale_x_discrete(limits=c("slow", "fast")) + theme_gray(base_size = 22)+theme(legend.title=element_blank()) #+
        # annotate("text",size=8, x=2, y=max(p0()$count)+0.1*max(p0()$count), label=paste(l()))+
        # ylab("TPM")+
        # xlab("")
    })
    output$TPMPlot <- renderPlot({
      (p2())
    }) 
}

shinyApp(ui = ui, server = server)