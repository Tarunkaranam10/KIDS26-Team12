# KIDS26 demonstration: precomputed, approved outputs only.
library(shiny)
fixture <- data.frame(sample_id=c("SYNTH-DEMO-A","SYNTH-DEMO-B","SYNTH-DEMO-OOD"),cancer_type=c("SYNTH-A","SYNTH-B","SYNTH-OOD"),split="synthetic_fixture",
  actual=c(20,50,NA),estimate_for_display=c(23,46,NA),lower=c(10,30,NA),upper=c(36,62,NA),
  reportable=c(TRUE,TRUE,FALSE),warning=c("synthetic illustration","synthetic illustration","outside_training_distribution"),
  provenance="SYNTHETIC ENGINEERING FIXTURE - NOT MODEL RESULTS")
path <- Sys.getenv("KIDS26_DEMO_RESULTS","")
results <- if(nzchar(path)) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE) else fixture
required <- c("sample_id","estimate_for_display","lower","upper","reportable","warning","provenance")
if(!all(required %in% names(results)))stop("Demo table missing fields; see docs/08_N_OF_1_ROADMAP.md")
if(anyDuplicated(results$sample_id))stop("Duplicate demo IDs")
ui <- fluidPage(titlePanel("KIDS26: methylation prediction of genomic-scar burden"),
  p("Research demonstration. Predicted reference HRDsum is not a treatment recommendation."),
  selectInput("sample","Approved sample",choices=results$sample_id),
  verbatimTextOutput("provenance"),h3(textOutput("score")),textOutput("interval"),
  textOutput("reference"),textOutput("residual"),verbatimTextOutput("qc"),
  verbatimTextOutput("validation"),tableOutput("by_cancer"),plotOutput("scatter"),plotOutput("distribution"),
  p("Method: frozen shared-probe elastic net. Intervals require independent calibration; unseen-domain coverage is not guaranteed. Unavailable quantities are intentionally omitted."))
server <- function(input,output,session){
  s<-reactive(results[results$sample_id==input$sample,,drop=FALSE])
  output$provenance<-renderText(s()$provenance)
  output$score<-renderText(if(isTRUE(s()$reportable[1]) && is.finite(s()$estimate_for_display[1])) paste("Predicted reference HRDsum:",round(s()$estimate_for_display[1],1)) else "Prediction unreliable / withheld")
  output$interval<-renderText(if(all(is.finite(c(s()$lower[1],s()$upper[1]))))paste("Illustrated interval:",round(s()$lower[1],1),"to",round(s()$upper[1],1)) else "Individual prediction interval unavailable")
  output$reference<-renderText(if("actual"%in%names(s()) && is.finite(s()$actual[1]))paste("Independent reference HRDsum:",round(s()$actual[1],1)) else "Independent reference HRDsum unavailable")
  output$residual<-renderText(if("actual"%in%names(s()) && is.finite(s()$actual[1]) && is.finite(s()$estimate_for_display[1]))paste("Residual (prediction - reference):",round(s()$estimate_for_display[1]-s()$actual[1],1)) else "Residual unavailable")
  output$qc<-renderText(s()$warning)
  output$validation<-renderText({
    if(!"actual"%in%names(results))return("Held-out performance unavailable")
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(sum(ok)<2)return("Held-out performance unavailable")
    err<-results$estimate_for_display[ok]-results$actual[ok]
    group<-if("cancer_type"%in%names(results)) paste("; cancers:",length(unique(results$cancer_type[ok]))) else ""
    paste0("Precomputed validation N=",sum(ok),group,"; MAE=",round(mean(abs(err)),2),"; RMSE=",round(sqrt(mean(err^2)),2))
  })
  output$by_cancer<-renderTable({
    if(!all(c("actual","cancer_type")%in%names(results)))return(NULL)
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(!any(ok))return(NULL)
    d<-results[ok,,drop=FALSE];d$absolute_error<-abs(d$estimate_for_display-d$actual)
    n<-aggregate(absolute_error~cancer_type,d,length);names(n)[2]<-"N"
    m<-aggregate(absolute_error~cancer_type,d,mean);names(m)[2]<-"MAE"
    o<-merge(n,m,by="cancer_type");o$MAE<-round(o$MAE,2);o[order(o$cancer_type),]
  })
  output$scatter<-renderPlot({
    if(!"actual"%in%names(results)){plot.new();text(.5,.5,"No independent reference values supplied");return()}
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(sum(ok)<2){plot.new();text(.5,.5,"Insufficient paired values");return()}
    lim<-range(c(results$actual[ok],results$estimate_for_display[ok]))
    plot(results$actual[ok],results$estimate_for_display[ok],xlab="Reference HRDsum",ylab="Predicted HRDsum",xlim=lim,ylim=lim,pch=19,col="#147d92",main=if(!nzchar(path))"Synthetic illustration only" else "Approved precomputed results");abline(0,1,lty=2,col="grey")
  })
  output$distribution<-renderPlot({
    ok<-is.finite(results$estimate_for_display)
    if(!any(ok)){plot.new();text(.5,.5,"No reportable predicted distribution");return()}
    hist(results$estimate_for_display[ok],breaks="FD",col="#7fc8d8",border="white",xlab="Predicted reference HRDsum",main=if(!nzchar(path))"Synthetic fixture distribution" else "Approved prediction distribution")
    if(isTRUE(s()$reportable[1])&&is.finite(s()$estimate_for_display[1]))abline(v=s()$estimate_for_display[1],col="#b2182b",lwd=2)
  })
}
shinyApp(ui,server)
