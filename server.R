library(shiny)
library(gplots)
library(ggplot2)
library(colorRamps)
library(survival)
source('heatmap-mik.R')

load('shiny-SB-SCC-microarray-20151116.RData')

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

shinyServer(function(input, output) {
    
  genes <- reactive({
    if(input$list=="gps"){
      g<-c("OSTC","MCM6","RPA3","MCM7","PCNA","XRCC6","KPNA2","ANLN","RNASEH2A",
           "PBK","GMNN","RRM1","CDC45","MAD2L1","RAN","DUT","RRM2","CDK7",
           "MLH3","SMC4","SMC3","POLD2","POLE2","BCCIP","GINS2","TREX1",
           "BUB3","FEN1","DBF4B","MOB4","CCNE1","RPA1","POLE3","RFC4","MCM3",
           "CHEK1","CCND1","CDC37")
    }
    if(input$list=="emt"){
#       g<-c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2")
# These are epithelial tissue markers:      
#      g<-c("CDH1","CLDN4","CLDN7","TJP3","MUC1")
      g<-c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2",
           "CDH1","CLDN4","CLDN7","TJP3","MUC1")
    }
    if(input$list=="th1"){
      g<-c("GZMB","GNLY","IFNG","IRF1","CD3Z","CD8A","TBX21","TNFRSF10A")
    }
    if(input$list=="topg"){
      g<-rownames(ttDat)[1:input$ngenes]
    }
    if(input$list=="user"){
      g<-strsplit(gsub(":"," ",gsub(","," ",gsub(";"," ",input$genes,fixed=T),
                                   fixed=T),fixed=T)," ")[[1]]
      g<-g[g!=""]
    }
    if(input$list=="gcor") g<-toupper(input$corGene)
    g<-unlist(lapply(tolower(g),.simpleCap))
    return(g)
  })
    
  keepSamples<-reactive({
    kp<-rep(TRUE,ncol(expDat))
#    if(input$samp=="pri") kp<-clinDat$tissue=="Primary Tumor"
#    if(input$samp=="dmt") kp<-clinDat$tissue=="Distant Metastasis"
#    if(input$samp=="rln") kp<-clinDat$tissue=="Regional Lymph Node"
#    if(input$samp=="rtm") kp<-clinDat$tissue=="Regional tissue"
#    if(input$samp=="met") kp<-clinDat$tissue%in%c("Distant Metastasis","Regional Lymph Node","Regional tissue")
#    if(input$samp=="braf") kp<-clinDat$brafMut==1
#    if(input$samp=="nras") kp<-clinDat$nrasMut==1
    return(kp)
  })

  eDat <- reactive({
    g<-genes()
    if(input$list=="gcor"){
      pCor<-cbind(expDat[match(unlist(annDat[match(g,annDat[,3]),1]),rownames(expDat)),])
      if(nrow(pCor)==1) pCor<-t(pCor)
      corDat<-cor(pCor,t(expDat))
      p<-rownames(expDat)[order(abs(corDat),decreasing=T)[sort(abs(corDat),decreasing=T)>=input$corThres]]
      g<-unlist(annDat[match(p,annDat[,1]),3])
      if(length(g)>200) g<-g[1:200]
    }
    mt<-as.vector(na.omit(match(g,annDat[,3])))
    g<-g[!is.na(match(g,annDat[,3]))]
    p<-unlist(annDat[mt,1])
    x<-expDat[match(p,rownames(expDat)),]
    rownames(x)<-g
    z<-x[apply(is.na(x),1,sum)==0,]
    return(z)
  })
    
  output$expMap <- renderPlot({    
    
    kp<-keepSamples()
    kp[is.na(kp)]<-FALSE
    
    # Metagene based on genes in heatmap
    gdat<-t(apply(eDat(),1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    ss<-rank(svd(gdat)$v[,1])/ncol(gdat)
    if(cor(apply(gdat,2,mean),ss)<0) ss<-abs(1-ss)
    
    gdat[gdat< -3]<- -3
    gdat[gdat>3]<-3
    
    gdat<-gdat[,kp]
    
    # Proliferation Metagene
    gps<-c("OSTC","MCM6","RPA3","MCM7","PCNA","XRCC6","KPNA2","ANLN","RNASEH2A",
           "PBK","GMNN","RRM1","CDC45","MAD2L1","RAN","DUT","RRM2","CDK7",
           "MLH3","SMC4","SMC3","POLD2","POLE2","BCCIP","GINS2","TREX1",
           "BUB3","FEN1","DBF4B","MOB4","CCNE1","RPA1","POLE3","RFC4","MCM3",
           "CHEK1","CCND1","CDC37")
    gps<-unlist(lapply(tolower(gps),.simpleCap))
    gpsProbes<-as.vector(na.omit(annDat[match(gps,annDat[,3]),1]))
    gpsDat<-t(apply(expDat[as.vector(na.omit(match(gpsProbes,rownames(expDat)))),],1,
                    function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    gpsMG<-rank(svd(gpsDat)$v[,1])/ncol(gpsDat)
    if(cor(apply(gpsDat,2,mean),gpsMG)<0) gpsMG<-abs(1-gpsMG)
    
    # EMT metagene
    emt<-c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2",
           "CDH1","CLDN4","CLDN7","TJP3","MUC1")
    emt<-unlist(lapply(tolower(emt),.simpleCap))
    emtProbes<-as.vector(na.omit(annDat[match(emt,annDat[,3]),1]))
    emtDat<-t(apply(expDat[as.vector(na.omit(match(emtProbes,rownames(expDat)))),],1,
                    function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    emtMG<-rank(svd(emtDat)$v[,1])/ncol(emtDat)
    if(cor(apply(emtDat,2,mean),emtMG)<0) emtMG<-abs(1-emtMG)

    ## TH1 metagene
    th1<-c("GZMB","GNLY","IFNG","IRF1","CD3Z","CD8A","TBX21","TNFRSF10A")
    th1<-unlist(lapply(tolower(th1),.simpleCap))
    th1Probes<-as.vector(na.omit(annDat[match(th1,annDat[,3]),1]))
    th1Dat<-t(apply(expDat[as.vector(na.omit(match(th1Probes,rownames(expDat)))),],1,
                    function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    th1MG<-rank(svd(th1Dat)$v[,1])/ncol(th1Dat)
    if(cor(apply(th1Dat,2,mean),th1MG)<0) th1MG<-abs(1-th1MG)
    
    ## Zmiz1 metagene
    zCor<-cor(expDat['Zmiz1',],t(expDat))
    zcorDat<-expDat[abs(zCor)>0.7,]
    zcorDat<-t(apply(zcorDat,1,function(x) 
      (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    zcorMG<-rank(svd(zcorDat)$v[,1])/ncol(zcorDat)
    if(cor(apply(zcorDat,2,mean),zcorMG)<0) zcorMG<-abs(1-zcorMG)

    if(input$colors=="greenred"){
      cc<-rbind(greenred(length(emtMG))[rank(emtMG)],
                greenred(length(th1MG))[rank(th1MG)],
                greenred(length(gpsMG))[rank(gpsMG)],
                greenred(length(gpsMG))[rank(zcorMG)],
                greenred(length(ss))[rank(ss)])
      ccx<-greenred(50)
    } 
    if(input$colors=="bluered"){
      cc<-rbind(bluered(length(emtMG))[rank(emtMG)],
                bluered(length(th1MG))[rank(th1MG)],
                bluered(length(gpsMG))[rank(gpsMG)],
                bluered(length(gpsMG))[rank(zcorMG)],
                bluered(length(ss))[rank(ss)])
      ccx<-bluered(50)
    }
    if(input$colors=="blueyellow"){
      cc<-rbind(blue2yellow(length(emtMG))[rank(emtMG)],
                blue2yellow(length(th1MG))[rank(th1MG)],
                blue2yellow(length(gpsMG))[rank(gpsMG)],
                blue2yellow(length(gpsMG))[rank(zcorMG)],
                blue2yellow(length(ss))[rank(ss)])
      ccx<-blue2yellow(50)
    }
    if(input$colors=="matlab"){
      cc<-rbind(matlab.like(length(emtMG))[rank(emtMG)],
                matlab.like(length(th1MG))[rank(th1MG)],
                matlab.like(length(gpsMG))[rank(gpsMG)],
                matlab.like(length(gpsMG))[rank(zcorMG)],
                matlab.like(length(ss))[rank(ss)])
      ccx<-matlab.like(50)
    }
#    if(input$slider){
#      ccy<-ifelse(ss>quantile(ss,input$survSlider),"red","green")
#    }
#    else{
#      ccy<-ifelse(ss>quantile(ss,0.667),"Red",ifelse(ss>quantile(ss,0.333),"Black","Green"))
#    }
    
#    cc<-rbind(cc,ccy)
    rownames(cc)<-c("EMT","Immune Response","Proliferation",
                    "Zmiz1 metagene","Heatmap metagene")
    oo<-order(ss[kp])
   
# Add Zmiz1 (+ other insertions) info
    mtc<-c(grey(0.7),"black")

    cc<-rbind(mtc[clinDat$Zmiz1.Status],
              t(apply(mutDat[-1,],1,function(x) mtc[x+1])),
              cc)[,kp]

    rownames(cc)[1]<-"Zmiz1 Status"
    
    if(input$order=="mt"){
      oo<-order(clinDat$Zmiz1.Status)
    }

    if(input$order!="cluster"){
        suppressWarnings(
        heatmap.mik(gdat[,oo], trace='none',
                    col=ccx, keysize=1.5, key=T,
                    ColSideColors=cc[,oo],mar=c(4,10),
                    labCol=colnames(expDat)[oo],Colv=F,cexCol=1.2)
        )
#        legend(0.18,1.05,nm,ncol=6,
#               fill=c(stc,grc,erc),bty='n')
      }
      if(input$order=="hc"){
        if(input$radioDist=="Euclidean"){
          suppressWarnings(
            heatmap.mik(gdat, trace='none',
                        col=ccx, keysize=1.5, key=T,
                        ColSideColors=cc,mar=c(4,10),
                        labCol=colnames(expDat),cexCol=1.2)
          )
        }
        if(input$radioDist=="Correlation"){
          suppressWarnings(
            heatmap.mik(gdat, trace='none',
                        col=ccx, keysize=1.5, key=T,
                        ColSideColors=cc,mar=c(4,10),
                        labCol=colnames(expDat),cexCol=1.2,
                        distfun=function(x) as.dist(1-cor(t(x))))
          )
        }
      }
  },height=600,width=1000)

  
  makeDataTable<-reactive({
    dt<-data.frame(Gene=rownames(eDat()),MGcor=mgCor(),CPHpv=geneCox())
    return(dt)
  })

  output$giDisplay <- renderDataTable({ 
    tt<-ttDat[match(rownames(eDat()),rownames(ttDat)),]
    data.frame(Gene=rownames(tt),
               FoldChange=round(2^tt$logFC,2),
               logFC=round(tt$logFC,2),
               AveLogExp=round(tt$AveExpr,2),
               RawP=round(tt$P.Value,5),
               FDRadjP=round(tt$adj.P.Val,5))
  }, options = list(orderClasses = TRUE,
                    autoWidth=FALSE,
                    pageLength = 200,
                    lengthChange=FALSE
                )
  )

  output$glDisplay <- renderTable({
    cbind(rownames(eDat()))
  },include.rownames=FALSE, include.colnames=FALSE)
  
})
