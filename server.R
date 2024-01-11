shinyServer(function(input, output,session) {

  output$slickr <- renderSlickR({
    imgs <- c("Home_protein.png","Home_pho.png")
    slickR(imgs,height = "350px")+ settings(slidesToShow=1,
                                            slidesToScroll= 1,
                                            autoplay = T,
                                            autoplaySpeed=5000)
  })
  
  observeEvent(input$dataset, {
    updateSelectizeInput(session, "symbol", choices = protein[[input$dataset]], server = TRUE)
  })
  observeEvent(input$pro_volcano_dataset, {
    updateSelectizeInput(session, "pro_volcano_gene", choices = protein[[input$pro_volcano_dataset]], server = TRUE)
  })
  observeEvent(input$corr_dataset, {
    updateSelectizeInput(session, "geneA", choices = protein[[input$corr_dataset]], server = TRUE)
    updateSelectizeInput(session, "geneB", choices = protein[[input$corr_dataset]], server = TRUE)
  })
  observeEvent(input$corr2_dataset, {
    updateSelectizeInput(session, "corr2_pro", choices = protein[[input$corr2_dataset]], server = TRUE)
  })
  observeEvent(input$survival_dataset, {
    updateSelectizeInput(session, "survival_gene", choices = protein[[input$survival_dataset]], server = TRUE)
  })
  observeEvent(input$pro_age_dataset, {
    updateSelectizeInput(session, "pro_age_protein", choices = protein[[input$pro_age_dataset]], server = TRUE)
  })
  observeEvent(input$pro_gender_dataset, {
    updateSelectizeInput(session, "pro_gender_protein", choices = protein[[input$pro_gender_dataset]], server = TRUE)
  })
  observeEvent(input$pro_stage_dataset, {
    updateSelectizeInput(session, "pro_stage_protein", choices = protein[[input$pro_stage_dataset]], server = TRUE)
  })
  
  observeEvent(input$mp_corr_dataset, {
    updateSelectizeInput(session, "mp_gene", choices = rna[[input$mp_corr_dataset]], server = TRUE)
  })
  
  observeEvent(input$pho_dataset, {
    updateSelectizeInput(session, "pho_site", choices = site[[input$pho_dataset]], server = TRUE)
  })
  
  observeEvent(input$pho_volcano_dataset, {
    updateSelectizeInput(session, "pho_volcano_site", choices = site[[input$pho_volcano_dataset]], server = TRUE)
  })
  
  observeEvent(input$pho_corr_dataset, {
    updateSelectizeInput(session, "pro_A", choices = protein[[input$pho_corr_dataset]], server = TRUE)
    updateSelectizeInput(session, "pho_A", choices = site[[input$pho_corr_dataset]], server = TRUE)
  })
  observeEvent(input$pho_survival_plot_dataset, {
    updateSelectizeInput(session, "pho_survival_plot_site", choices = site[[input$pho_survival_plot_dataset]], server = TRUE)
  })
  observeEvent(input$pho_survival_table_dataset, {
    updateSelectizeInput(session, "pho_survival_table_protein", choices = pho_survival_protein[[input$pho_survival_table_dataset]], server = TRUE)
  })
 
  observeEvent(input$pho_age_dataset, {
    updateSelectizeInput(session, "pho_age_site", choices = site[[input$pho_age_dataset]], server = TRUE)
  })
  observeEvent(input$pho_gender_dataset, {
    updateSelectizeInput(session, "pho_gender_site", choices = site[[input$pho_gender_dataset]], server = TRUE)
  })

  observeEvent(input$pho_stage_dataset, {
    updateSelectizeInput(session, "pho_stage_site", choices = site[[input$pho_stage_dataset]], server = TRUE)
  })

  observeEvent(input$pho_ks_corr_dataset, {
    updateSelectizeInput(session, "kinase", choices = Kinase[[input$pho_ks_corr_dataset]], server = TRUE)
  })

  plot <-  eventReactive(input$action, {
    protein <- read.table(paste0("protein_",input$dataset,".txt"),header = T,sep = "\t",check.names = F)
    class <- read.table(paste0("patient_",input$dataset,".txt"),header = T,sep = "\t")
    pro_input<- as.data.frame(t(protein[input$symbol,]))
    pro_input$tumor <- class$class
    pro_input$tumor <- factor(pro_input$tumor,levels = c("Tumor","Normal"))
    
    boxplot <- ggboxplot(pro_input, x = "tumor", y = input$symbol,color = "tumor",width = 0.4,
                         add = "jitter",size = 1.2,palette =c("Tumor"=input$pro_col_tumor,"Normal" =input$pro_col_normal),add.params=list(size=input$pro_df_size,jitter = 0.1))+
      stat_compare_means(method = input$pro_df_method,hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$dataset)+xlab(NULL)+ylab(paste0(input$symbol," Abundance\n(Log2 TMT ratio)"))+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
      
    boxplot
  })
  
  output$plot <- renderPlot({ 
    req(input$action)
    plot()
  })

  output$downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$dataset,"_",input$symbol,"_boxplot.pdf")
    },
    content = function (file) {
      ggsave(file, plot = plot(),he=8,wi=6)
    })
	

  pro_volcano_plot <- eventReactive(input$pro_volcano_action, {
    volcanoData <- read.csv(paste0("DEP.",input$pro_volcano_dataset,".csv"),row.names = 1)
    markerGene <- volcanoData[rownames(volcanoData)== input$pro_volcano_gene,]
    markerGene$labelGene <- rownames(markerGene)
    upGene <- rownames(volcanoData)[which(volcanoData$adj.P.Val< input$pro_volcano_FDR & volcanoData$logFC>input$pro_volcano_log2FC)]
    downGene <- rownames(volcanoData)[which(volcanoData$adj.P.Val< input$pro_volcano_FDR & volcanoData$logFC< -(input$pro_volcano_log2FC))]
    
    normdata <- volcanoData[-which(rownames(volcanoData) %in% c(upGene,downGene)),]
    updata <- volcanoData[which(rownames(volcanoData) %in% upGene),]
    downdata <- volcanoData[which(rownames(volcanoData) %in% downGene),]
    
    p1 <- ggplot()+
      geom_point(pch=20,cex=0.6,color="gray77",
                 data=normdata,aes(x=logFC,y=(-log10(adj.P.Val))))+
      geom_point(pch=20,cex=0.6,color="#FC4E07",
                 data=updata,aes(x=logFC,y=(-log10(adj.P.Val))))+
      geom_point(pch=20,cex=0.6,color="#00AFBB",
                 data=downdata,aes(x=logFC,y=(-log10(adj.P.Val))))+
      geom_point(pch=20,cex=3,color="red",
                 data=markerGene,aes(x=logFC,y=(-log10(adj.P.Val))))+
      theme_classic()+
      geom_vline(xintercept = c(-(input$pro_volcano_log2FC),input$pro_volcano_log2FC),linetype="dotted",size = 0.7)+
      geom_hline(yintercept = c(-log10(input$pro_volcano_FDR)),linetype="dotted",size = 0.7)+
      xlab("log2 (fold change)")+ylab("-log10 (FDR)")

    p1
  })
  
  output$pro_volcano_plot <- renderPlot({ 
    input$pro_volcano_action
    pro_volcano_plot()
  })
  
  output$pro_volcano_plot_download <- downloadHandler(
    filename = function() {
      paste0(input$pro_volcano_dataset,"_",input$pro_volcano_gene,"_volcano.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pro_volcano_plot(),he=6,wi=6)
    })

  pro_volcano_table <- eventReactive(input$pro_volcano_action, {
    DEP <- read.csv(paste0("DEP.",input$pro_volcano_dataset,".csv"),row.names = 1)
    DEP <- DEP[order(abs(DEP$logFC),decreasing = T),]
    DEP
  })
 
  output$pro_volcano_table <- DT::renderDataTable({
    input$pro_volcano_action
    formatSignif(datatable(pro_volcano_table()),columns=1:6,digits=3)
  })

  output$pro_volcano_csv_download <- downloadHandler(
    filename = function() {
      paste0(input$pro_volcano_dataset,"_volcano.csv")
    },
    content = function (file) {
      write.csv(pro_volcano_table(),file,row.names = T,quote = F)
    })
  
  corr_plot <- eventReactive(input$action2, {

    protein_tumor_nonImpute <- read.table(paste0(input$corr_dataset,"/","protein_Tumor_",input$corr_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    ab_nonImpute <- na.omit(as.data.frame(t(protein_tumor_nonImpute[c(input$geneA,input$geneB),])))
    
    corr_nonImpute <-ggscatter(ab_nonImpute,x = input$geneA, y = input$geneB,
                                conf.int = TRUE,size = input$pro_corr_size,color = input$pro_col_non_imputed_corr,title = "Filter missing value",
                                xlab = paste0(input$geneA," Abundance\n(Log2 TMT ratio)"),
                                ylab = paste0(input$geneB," Abundance\n(Log2 TMT ratio)"))+
      stat_cor(method = input$correlation,size=5)+geom_smooth(method = "lm",linewidth=1.2,colour="black")+
      theme(plot.title = element_text(hjust = 0.5))

    protein_tumor_impute <- read.table(paste0(input$corr_dataset,"/","protein_Tumor_",input$corr_dataset,"_Impute.txt"),header = T,sep = "\t",check.names=F)
    ab_impute <- as.data.frame(t(protein_tumor_impute[c(input$geneA,input$geneB),]))
    ab_impute$impute <- ifelse(rownames(ab_impute)%in%rownames(ab_nonImpute),"Non_imputed","Imputed")
    ab_impute$impute <- as.factor(ab_impute$impute) 
    corr_impute <- ggscatter(ab_impute, x = input$geneA, y = input$geneB,
                conf.int = TRUE,size = input$pro_corr_size,color  = "impute",palette = c(Imputed=input$pro_col_imputed_corr, Non_imputed=input$pro_col_non_imputed_corr),
                title = "Impute missing value",show.legend.text = T,
                xlab = paste0(input$geneA," Abundance\n(Log2 TMT ratio)"),
                ylab = paste0(input$geneB," Abundance\n(Log2 TMT ratio)"))+
      stat_cor(method = input$correlation,size=5)+geom_smooth(method = "lm",linewidth=1.2,colour="black")+
      theme(plot.title = element_text(hjust = 0.5))
    
   arrange <- ggarrange(corr_nonImpute,corr_impute,  nrow = 1, ncol = 2,legend = "right",common.legend = T)
   arrange
  })
  
  output$corr_plot <- renderPlot({
    input$action2
    corr_plot()
  })

  output$downloadPdf2 <- downloadHandler(
    filename = function() {
      paste0(input$corr_dataset,"_",input$geneA,"_",input$geneB,"_corr",".pdf")
    },
    content = function (file) {
      ggsave(file, plot = corr_plot(),he=6,wi=12)
    })
    
    
  corr_table_impute <- eventReactive(input$corr2_action, {
    protein_tumor_impute <- read.table(paste0(input$corr2_dataset,"/","protein_Tumor_",input$corr2_dataset,"_Impute.txt"),header = T,sep = "\t")
    pro_imputed_fraction <- read.csv(paste0(input$corr2_dataset,"/protein_imputed_fraction_",input$corr2_dataset,".csv"),row.names = 1,check.names=F)
    corr <- as.data.frame(t(protein_tumor_impute))

    corr_table_impute <- do.call(rbind, lapply(1:ncol(corr), function(x){
      dd = cor.test(corr[,input$corr2_pro], corr[,x], method = input$corr2_correlation )
      data.frame(Protein = input$corr2_pro, Correlated_proteins = colnames(corr)[x], Coefficient = signif(dd$estimate, digits = 4), p.value = signif(dd$p.value, digits = 4))
    }))
    corr_table_impute$Imputed_fraction <- pro_imputed_fraction$`Imputed_fraction(%)`
    corr_table_impute <- corr_table_impute[-c(which(corr_table_impute$Correlated_proteins==input$corr2_pro)),]
    corr_table_impute$adjP.value <- signif(p.adjust(corr_table_impute$p.value,method = "BH"),digits = 4)
    corr_table_impute <- corr_table_impute[order(corr_table_impute$Coefficient,decreasing = T),]
    corr_table_impute <- corr_table_impute[,c(1,2,5,3,4,6)]
    colnames(corr_table_impute) <- c("Protein","Correlated proteins","Imputed fraction(%)","Correlation coefficient","P.value","adjP.value")
    corr_table_impute
  })

  output$corr2_table_impute <- DT::renderDataTable({
    input$corr2_action
    DT::datatable(corr_table_impute(),rownames = F)
  })
  
  output$corr2_table_impute_download <- downloadHandler(
    filename = function() {
      paste0(input$corr2_dataset,"_",input$corr2_pro,"_impute_corr",".txt")
    },
    content = function (file) {
      write.table(corr_table_impute(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })

  corr_table_nonImpute <- eventReactive(input$corr2_action, {
    protein_tumor_nonImpute <- read.table(paste0(input$corr2_dataset,"/","protein_Tumor_",input$corr2_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    
    corr <- as.data.frame(t(protein_tumor_nonImpute))

    corr_table_nonImpute <- do.call(rbind, lapply(1:ncol(corr), function(x){
      dd = cor.test(corr[,input$corr2_pro], corr[,x], method = input$corr2_correlation,use="pairwise.complete.obs")
      data.frame(Protein = input$corr2_pro, Correlated_proteins = colnames(corr)[x], Coefficient = signif(dd$estimate, digits = 4), p.value = signif(dd$p.value, digits = 4))
    }))
    corr_table_nonImpute <- corr_table_nonImpute[-c(which(corr_table_nonImpute$Correlated_proteins==input$corr2_pro)),]
    corr_table_nonImpute$adjP.value <- signif(p.adjust(corr_table_nonImpute$p.value,method = "BH"),digits = 4)
    corr_table_nonImpute <- corr_table_nonImpute[order(corr_table_nonImpute$Coefficient,decreasing = T),]
    colnames(corr_table_nonImpute) <- c("Protein","Correlated proteins","Correlation coefficient","P.value","adjP.value")
    corr_table_nonImpute
  })

  output$corr2_table_nonImpute <- DT::renderDataTable({
    input$corr2_action
    DT::datatable(corr_table_nonImpute(),rownames = F)
  })
  
  output$corr2_table_nonImpute_download <- downloadHandler(
    filename = function() {
      paste0(input$corr2_dataset,"_",input$corr2_pro,"_nonImpute_corr",".txt")
    },
    content = function (file) {
      write.table(corr_table_nonImpute(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })
  

  survival_plot <- eventReactive(input$action3, {
    protein_tumor <- read.table(paste0(input$survival_dataset,"/","protein_Tumor_",input$survival_dataset,".txt"),header = T,sep = "\t",check.names=F)
    clinical <- read.table(paste0(input$survival_dataset,"/","clinical_",input$survival_dataset,".txt"),header = T,sep = "\t")
    clinical$gene <- as.numeric(protein_tumor[input$survival_gene,])
    clinical <- clinical[!is.na(clinical$OS_days),]
    if (input$cutoff == "Median") {
      clinical$Group <- ifelse(clinical$gene >= median(clinical$gene,na.rm = T),"High","Low")
      lab1 <- paste0(names(table(clinical$Group)[1])," (N=",table(clinical$Group)[1],")")
      lab2 <- paste0(names(table(clinical$Group)[2])," (N=",table(clinical$Group)[2],")")
      fit <- survfit(Surv(OS_days,Status) ~ Group,  
                     data = clinical)
      survival <- ggsurvplot(fit, data = clinical,palette = c(input$pro_col_high,input$pro_col_low),pval = TRUE, conf.int = T,conf.int.style = "step",
                             xlab = "Follow up time (days)",ylab = "Overall survival rate",legend.title = "",title = input$survival_gene,
                             legend = c(0.28,0.18),pval.size=5.2,legend.labs= c(lab1,lab2))%++%
        theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
              legend.background = element_blank())+
        guides(color=guide_legend(override.aes = list(shape = NA)))
      
    } else if (input$cutoff == "Optimal value") {
      res.cut <- surv_cutpoint(clinical, 
                               time = "OS_days", 
                               event = "Status", 
                               variables = c("gene"), 
                               minprop = 0.3
      )

        res.cat <- surv_categorize(res.cut,labels = c("Low", "High"))
        lab_high <- paste0(names(table(res.cat$gene)[1])," (N=",table(res.cat$gene)[1],")")
        lab_low <- paste0(names(table(res.cat$gene)[2])," (N=",table(res.cat$gene)[2],")")
        fit2 <- survfit(Surv(OS_days, Status) ~gene, data = res.cat)
        data.optimal <- res.cat
      

      survival <- ggsurvplot(fit2, data = data.optimal,palette = c(input$pro_col_high,input$pro_col_low),pval = TRUE, conf.int = T,conf.int.style = "step",
                             xlab = "Follow up time (days)",ylab = "Overall survival rate",legend.title = "",title = input$survival_gene,
                             legend = c(0.28,0.18),pval.size=5.2,
                             legend.labs= c(lab_high,lab_low))%++%
        theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
              legend.background = element_blank())+
        guides(color=guide_legend(override.aes = list(shape = NA)))
    }
    survival
  })
  
  output$survival_plot <- renderPlot({
    input$action3
    survival_plot()
  })

  output$downloadPdf3 <- downloadHandler(
    filename = function() {
      paste0(input$survival_dataset,"_",input$survival_gene,"_survival",".pdf")
    },
    content = function (file) {
      survival_plot()
      ggsave(file,he=6,wi=6)
    })

  mp_corr_plot <- eventReactive(input$action4, {
    protein_tumor_nonImpute <- read.table(paste0(input$mp_corr_dataset,"/","protein_Tumor_",input$mp_corr_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    protein_tumor_nonImpute <- as.data.frame(t(protein_tumor_nonImpute[input$mp_gene,]))
    rna <- read.table(paste0(input$mp_corr_dataset,"/","rna_",input$mp_corr_dataset,".txt"),header = T,sep = "\t",check.names=F)
    rna <- as.data.frame(t(rna[input$mp_gene,]))
    rna_pro_nonImpute <- merge(rna,protein_tumor_nonImpute,by = "row.names")
    colnames(rna_pro_nonImpute) <- c("patient","mRNA","Protein")
    rna_pro_nonImpute <- na.omit(rna_pro_nonImpute)
    
    mp_corr_nonImpute <-ggscatter(rna_pro_nonImpute,x = "mRNA", y = "Protein",
                               conf.int = TRUE,size = input$mp_corr_size,color = input$mp_col_non_imputed_corr,title = "Filter missing value",
                               subtitle= input$mp_gene,
                               xlab = "mRNA abundance",
                               ylab = "Protein abundance")+
      stat_cor(method = input$mp_correlation,size=5)+geom_smooth(method = "lm",linewidth=1.2,colour="black")+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
 
    protein_tumor_Impute <- read.table(paste0(input$mp_corr_dataset,"/","protein_Tumor_",input$mp_corr_dataset,"_Impute.txt"),header = T,sep = "\t",check.names=F)
    protein_tumor_Impute <- as.data.frame(t(protein_tumor_Impute[input$mp_gene,]))
    rna_pro_Impute <- merge(rna,protein_tumor_Impute,by = "row.names")
    colnames(rna_pro_Impute) <- c("patient","mRNA","Protein")
    rna_pro_Impute$impute <- ifelse(rna_pro_Impute$patient%in%rna_pro_nonImpute$patient,"Non_imputed","Imputed")
    rna_pro_Impute$impute <- as.factor(rna_pro_Impute$impute) 

    mp_corr_Impute <- ggscatter(rna_pro_Impute, x = "mRNA", y = "Protein",
                                conf.int = TRUE,size = input$mp_corr_size,color  = "impute",palette = c(Imputed=input$mp_col_imputed_corr, Non_imputed=input$mp_col_non_imputed_corr),
                                title = "Impute missing value",show.legend.text = T,
                                subtitle= input$mp_gene,
                                xlab = "mRNA abundance",
                                ylab = "Protein abundance")+
      stat_cor(method = input$mp_correlation,size=5)+geom_smooth(method = "lm",linewidth=1.2,colour="black")+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    mp_corr_Impute
    arrange <- ggarrange(mp_corr_nonImpute,mp_corr_Impute,  nrow = 1, ncol = 2,legend = "right",common.legend = T)
    arrange
  })
  
  output$mp_corr_plot <- renderPlot({
    input$action4
    mp_corr_plot()
  })

  output$downloadPdf4 <- downloadHandler(
    filename = function() {
      paste0(input$mp_corr_dataset,"_",input$mp_gene,"_","mRNA_Protein_Corr",".pdf")
    },
    content = function (file) {
      ggsave(file, plot = mp_corr_plot(),he=6,wi=12)
    })
	

  pro_age_plot <-  eventReactive(input$pro_age_action, {
    protein <- read.table(paste0("protein_Tumor_",input$pro_age_dataset,".txt"),header = T,sep = "\t",check.names = F)
    clinical <- read.table(paste0("clinical_",input$pro_age_dataset,".txt"),header = T,sep = "\t")
    gene_input<- as.data.frame(t(protein[input$pro_age_protein,]))
    gene_input$age <- clinical$Age
    
    if (input$pro_age_cutoff == "median") {
      cutoff <- median(gene_input$age)
    }else if (input$pro_age_cutoff == "custom") {
      cutoff <- input$pro_age_custom
    }
    gene_input$Group <- ifelse(gene_input$age >= cutoff,"Older","Younger")
    gene_input$Group <- factor(gene_input$Group,levels = c("Younger","Older"))
    
    violinplot <- ggviolin(gene_input, x = "Group", y = input$pro_age_protein,fill = "Group",width = 0.6,ylab = paste0(input$pro_age_protein," Abundance\n(Log2 TMT ratio)"),
                           add = "boxplot",size = 0.8,palette =c("Younger"=input$pro_age_col_young,"Older" =input$pro_age_col_old),add.params = list(width=0.1,size=0.7,fill="white")) +
      stat_compare_means(method = input$pro_age_method,hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pro_age_dataset)+xlab(NULL)+scale_x_discrete(labels = c(paste0("Younger\n(<",cutoff,")"),paste0("Older\n(>=",cutoff,")")))+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    violinplot
  })
  
  output$pro_age_plot <- renderPlot({
    input$pro_age_action
    pro_age_plot()
  })

  output$pro_age_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_age_dataset,"_",input$pro_age_protein,"_Age.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pro_age_plot(),he=8,wi=6)
    })

  pro_gender_plot <-  eventReactive(input$pro_gender_action, {
    protein <- read.table(paste0("protein_Tumor_",input$pro_gender_dataset,".txt"),header = T,sep = "\t",check.names = F)
    clinical <- read.table(paste0("clinical_",input$pro_gender_dataset,".txt"),header = T,sep = "\t")
    gene_input<- as.data.frame(t(protein[input$pro_gender_protein,]))

    gene_input$Gender <- clinical$Gender
    gene_input$Gender <- factor(gene_input$Gender)
    
    violinplot <- ggviolin(gene_input, x = "Gender", y = input$pro_gender_protein,fill = "Gender",width = 0.6,ylab = paste0(input$pro_gender_protein," Abundance\n(Log2 TMT ratio)"),
                           add = "boxplot",size = 0.8,palette =c("Female"=input$pro_gender_col_female,"Male" =input$pro_gender_col_male),add.params = list(width=0.1,size=0.7,fill="white")) +
      stat_compare_means(method = input$pro_gender_method,hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pro_gender_dataset)+xlab(NULL)+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    violinplot
  })
  
  output$pro_gender_plot <- renderPlot({
    input$pro_gender_action
    pro_gender_plot()
  })

  output$pro_gender_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_gender_dataset,"_",input$pro_gender_protein,"_Gender.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pro_gender_plot(),he=8,wi=6)
    })

  pro_stage_plot <-  eventReactive(input$pro_stage_action, {
    protein <- read.table(paste0("protein_Tumor_",input$pro_stage_dataset,".txt"),header = T,sep = "\t",check.names = F)
    clinical <- read.table(paste0("clinical_",input$pro_stage_dataset,".txt"),header = T,sep = "\t")
    gene_input<- as.data.frame(t(protein[input$pro_stage_protein,]))

    gene_input$Stage <- clinical$Tumor_stage
    gene_input$Stage <- factor(gene_input$Stage)
    
    violinplot <- ggviolin(gene_input, x = "Stage", y = input$pro_stage_protein,fill = "Stage",width = 0.7,palette = "RdYIGn",
                           ylab = paste0(input$pro_stage_protein," Abundance\n(Log2 TMT ratio)"),
                           add = "boxplot",size = 0.8,add.params = list(width=0.1,size=0.7,fill="white")) +
      stat_compare_means(method = "anova",hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pro_stage_dataset)+xlab(NULL)+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    violinplot
  })
  
  output$pro_stage_plot <- renderPlot({
    input$pro_stage_action
    pro_stage_plot()
  })

  output$pro_stage_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_stage_dataset,"_",input$pro_stage_protein,"_Stage.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pro_stage_plot(),he=8,wi=12)
    })

  pro_go_table <-  eventReactive(input$pro_go_action, {
    DEP <- read.csv(paste0("DEP.",input$pro_go_dataset,".csv"),row.names = 1)
    log2FC_cutoff = log2(input$pro_go_FC)
    pvalue_cutoff = input$pro_go_adjP
    if(input$pro_go_regulation=="Up"){
      DEP <- DEP[DEP$logFC>log2FC_cutoff&DEP$adj.P.Val<pvalue_cutoff,]
    }else{
      DEP <- DEP[DEP$logFC< c(-log2FC_cutoff)&DEP$adj.P.Val<pvalue_cutoff,]
    }
    DEP
  })
 
  output$pro_go_table <- DT::renderDataTable({
    input$pro_go_action
    formatSignif(datatable(pro_go_table()),columns=1:6,digits=3)
  })
  
  go <- eventReactive(input$pro_go_plot, {
    DEP_symbol <- row.names(pro_go_table())
    gene_ID <- bitr(DEP_symbol, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
    
    go.bp <- enrichGO(gene_ID$ENTREZID, OrgDb = "org.Hs.eg.db",ont="BP") 
    if(!is.null(go.bp)){
      go.bp <- clusterProfiler::simplify(go.bp)
      go.bp <- setReadable(go.bp, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
      go.bp <- data.frame(go.bp) 
      go.bp$ONTOLOGY <- "BP"
      #order by p-value value
      go.bp <- go.bp[order(go.bp$pvalue), ]
      go.bp$GeneRatio <- as.numeric(sub("(\\d+)/\\d+","\\1",go.bp$GeneRatio))/as.numeric(sub("\\d+/(\\d+)","\\1",go.bp$GeneRatio))
    }
    
    go.cc <- enrichGO(gene_ID$ENTREZID, OrgDb = "org.Hs.eg.db",ont="CC") 
    if(!is.null(go.cc)){
      go.cc <- clusterProfiler::simplify(go.cc)
      go.cc <- setReadable(go.cc, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
      go.cc <- data.frame(go.cc)
      go.cc$ONTOLOGY <- "CC"
      go.cc <- go.cc[order(go.cc$pvalue), ]
      go.cc$GeneRatio <- as.numeric(sub("(\\d+)/\\d+","\\1",go.cc$GeneRatio))/as.numeric(sub("\\d+/(\\d+)","\\1",go.cc$GeneRatio))
      
    }
    
    go.mf <- enrichGO(gene_ID$ENTREZID, OrgDb = "org.Hs.eg.db",ont="MF") 
    if(!is.null(go.mf)){
      go.mf <- clusterProfiler::simplify(go.mf)
      go.mf <- setReadable(go.mf, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
      go.mf <- data.frame(go.mf) 
      go.mf$ONTOLOGY <- "MF"
      go.mf <- go.mf[order(go.mf$pvalue), ]
      go.mf$GeneRatio <- as.numeric(sub("(\\d+)/\\d+","\\1",go.mf$GeneRatio))/as.numeric(sub("\\d+/(\\d+)","\\1",go.mf$GeneRatio))
    }
    
    go <- rbind(go.bp, go.cc, go.mf)
    if(is.null(go)){ 
      shinyWidgets::sendSweetAlert(session = session ,type = "error",title = "Error!" ,text = "No GO terms are significantly enriched.")
    }
    go
  })
  
  go.top <- eventReactive(input$pro_go_plot, {
    go.bp <- go()[go()$ONTOLOGY == "BP", ]
    go.cc <- go()[go()$ONTOLOGY == "CC", ]
    go.mf <- go()[go()$ONTOLOGY == "MF", ]
    go.bp.top <- go.bp[1:min(nrow(go.bp), 10),]
    go.cc.top <- go.cc[1:min(nrow(go.cc), 10),]
    go.mf.top <- go.mf[1:min(nrow(go.mf), 10),]
    go.top <- rbind(go.bp.top,go.cc.top,go.mf.top)
    go.top$Description <- factor(go.top$Description,levels = rev(go.top$Description))
    go.top
  })

  output$pro_go_list_download <- downloadHandler(
    filename = function() {
      paste0(input$pro_go_dataset,"_GO_Result.txt")
    },
    content = function (file) {
      write.table(go(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })

  go.bar <- eventReactive(input$pro_go_plot, {
    go.bar <- ggplot(data = go.top(), aes(x = Description, y = -log10(pvalue), fill = ONTOLOGY)) +
      geom_col(width = 0.9) + 
      coord_flip() + 
      theme_classic() + 
      scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
      labs(x = "", y = "-log10(pvalue)", title = "GO Enrichment") + 
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11))+ 
      scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a"), 
                        breaks = c("BP", "CC", "MF"),
                        labels = c("BP", "CC", "MF")) 
    
    go.bar
  })
  
  output$pro_go_barplot <- renderPlot({
    input$pro_go_plot
    go.bar()
  })
  
  output$pro_go_bar_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_go_dataset,"_GO_Barplot.pdf")
    },
    content = function (file) {
      ggsave(file,plot = go.bar(),he=8,wi=9)
    })

  go.dotplot <- eventReactive(input$pro_go_plot, {
    go.top.dot <- go.top()[order(go.top()$ONTOLOGY,-go.top()$GeneRatio),]
    go.top.dot$Description <- factor(go.top.dot$Description,levels = rev(go.top.dot$Description))
    go_bubble <- ggplot(go.top.dot, aes(x = GeneRatio, y = Description, color = -log10(pvalue))) +
      geom_point(aes(size = Count)) + 
      theme_bw() + 
      scale_y_discrete(labels = function(y) str_wrap(y, width = 50)) + 
      labs(size = "Counts", x = "GeneRatio", y = "", title = "GO Enrichment") +
      scale_color_gradient(low="blue",high ="red") + 
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + 
      facet_grid(ONTOLOGY~., scales = "free", space = "free") 
    go_bubble
  })
  
  output$pro_go_dotplot <- renderPlot({
    input$pro_go_plot
    go.dotplot()
  })
  
  output$pro_go_dot_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_go_dataset,"_GO_Dotplot.pdf")
    },
    content = function (file) {
      ggsave(file,plot = go.dotplot(),he=8,wi=8)
    })

  pro_kegg_table <-  eventReactive(input$pro_kegg_action, {
    DEP <- read.csv(paste0("DEP.",input$pro_kegg_dataset,".csv"),row.names = 1)
    log2FC_cutoff = log2(input$pro_kegg_FC)
    pvalue_cutoff = input$pro_kegg_adjP
    
    if(input$pro_kegg_regulation=="Up"){
      DEP <- DEP[DEP$logFC>log2FC_cutoff&DEP$adj.P.Val<pvalue_cutoff,]
    }else{
      DEP <- DEP[DEP$logFC< c(-log2FC_cutoff)&DEP$adj.P.Val<pvalue_cutoff,]
    }
    DEP
  })

  output$pro_kegg_table <- DT::renderDataTable({
    input$pro_kegg_action
    formatSignif(datatable(pro_kegg_table()),columns=1:6,digits=3)
  })
  
  
  kegg <- eventReactive(input$pro_kegg_plot, {
    
    DEP_symbol <- row.names(pro_kegg_table())
    gene_ID <- bitr(DEP_symbol, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
    
    kegg <- enrichKEGG(gene_ID$ENTREZID,kk)
    kegg <- setReadable(kegg, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    
    if(dim(kegg)[1]==0){ 
      shinyWidgets::sendSweetAlert(session = session ,type = "error",title = "Error!" ,text = "No KEGG terms are significantly enriched.")
    }
    
    kegg <- data.frame(kegg)
    kegg <- kegg[order(kegg$pvalue), ]
    kegg$GeneRatio <- as.numeric(sub("(\\d+)/\\d+","\\1",kegg$GeneRatio))/as.numeric(sub("\\d+/(\\d+)","\\1",kegg$GeneRatio))
    kegg
  })
  
  kegg.top <- eventReactive(input$pro_kegg_plot, {
    kegg.top <- head(kegg(),min(nrow(kegg()),10))
    kegg.top$Description <- factor(kegg.top$Description,levels = rev(kegg.top$Description))
    kegg.top
  })

  output$pro_kegg_list_download <- downloadHandler(
    filename = function() {
      paste0(input$pro_kegg_dataset,"_KEGG_Result.txt")
    },
    content = function (file) {
      write.table(kegg(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })

  kegg.bar <- eventReactive(input$pro_kegg_plot, {
    kegg.bar <- ggplot(data = kegg.top(), aes(x = Description, y = -log10(pvalue), fill = "Description")) +
      geom_col(width = 0.9) + 
      coord_flip() + 
      theme_classic() + 
      scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) + 
      labs(x = "", y = "-log10(pvalue)", title = "KEGG Enrichment") + 
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
            legend.position="none")+
      scale_fill_manual(values = c("#e41a1c"))
    
    kegg.bar
  })
  
  output$pro_kegg_barplot <- renderPlot({
    input$pro_kegg_plot
    kegg.bar()
  })
  
  output$pro_kegg_bar_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_kegg_dataset,"_KEGG_Barplot.pdf")
    },
    content = function (file) {
      ggsave(file,plot = kegg.bar(),he=7,wi=7)
    })
	
  kegg.dotplot <- eventReactive(input$pro_kegg_plot, {
    kegg.top.dot <- kegg.top()[order(-kegg.top()$GeneRatio),]
    kegg.top.dot$Description <- factor(kegg.top.dot$Description,levels = rev(kegg.top.dot$Description))
    
    kegg_bubble <- ggplot(kegg.top.dot, aes(x = GeneRatio, y = Description, color = -log10(pvalue))) +
      geom_point(aes(size = Count)) + 
      theme_bw() +
      scale_y_discrete(labels = function(y) str_wrap(y, width = 40)) + 
      labs(size = "Counts", x = "GeneRatio", y = "", title = "KEGG Enrichment") +
      scale_color_gradient(low="blue",high ="red") + 
      scale_size_continuous(range=c(4,10))+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth =1.5, linetype="solid"),
            panel.grid = element_line(linewidth =1, linetype="solid"),
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))  
    kegg_bubble
  })
  
  output$pro_kegg_dotplot <- renderPlot({
    input$pro_kegg_plot
    kegg.dotplot()
  })
  
  output$pro_kegg_dot_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_kegg_dataset,"_KEGG_Dotplot.pdf")
    },
    content = function (file) {
      ggsave(file,plot = kegg.dotplot(),he=8,wi=8)
    })
	
  gsea <- eventReactive(input$pro_gsea_action, {
    set.seed(123456789)
    gsea_gene <- strsplit(input$pro_gsea_gene, "\\s*,\\s*")[[1]]
    gsea_gene <- gsub(" ","",gsea_gene)
    DEP <- read.csv(paste0("DEP.",input$pro_gsea_dataset,".csv"),row.names = 1)
    DEP <- DEP %>% arrange(desc(logFC))
    gene_list <- DEP$logFC
    names(gene_list) <- rownames(DEP)
    gene_set <- data.frame(term = rep("Targer_GeneSet", length(gsea_gene)), gene = gsea_gene)
    res <- GSEA(gene_list, TERM2GENE = gene_set, pvalueCutoff = 1,minGSSize = 5)
    return(res)
  })
  
  output$pro_gsea_plot <- renderPlot({
    input$pro_gsea_action
    if (dim(gsea()@result)[1] == 0) {
      showModal(modalDialog(
        title = "GSEA Analysis",
        "The input mapped proteins are less than 5, and the GSEA result is empty.",
        footer = modalButton("Close"),
        easyClose = TRUE
      ))
    }else{
      pro_gsea_plot <- gseaplot2(gsea(), title = "Target_GeneSet", geneSetID = 1, pvalue_table = T)
      pro_gsea_plot
    }
  })
  
  output$pro_gsea_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pro_gsea_dataset,"_GSEA.pdf")
    },
    content = function (file) {
      ggsave(file,he=6,wi=8)
    })
	
  pro_ppi_table <-  eventReactive(input$pro_ppi_action, {
    DEP <- read.csv(paste0("PPI.",input$pro_ppi_dataset,".csv"))
    
    log2FC_cutoff = log2(input$pro_ppi_FC)
    pvalue_cutoff = input$pro_ppi_adjP
    
    if(input$pro_ppi_regulation=="Up"){
      DEP <- DEP[DEP$logFC>log2FC_cutoff&DEP$adj.P.Val<pvalue_cutoff,]
    }else{
      DEP <- DEP[DEP$logFC< c(-log2FC_cutoff)&DEP$adj.P.Val<pvalue_cutoff,]
    }
    
    dat <- DEP[1:min(nrow(DEP),200),] 
    dat
  })

  output$pro_ppi_table <- DT::renderDataTable({
    input$pro_ppi_action
    formatSignif(datatable(pro_ppi_table(),rownames = FALSE),columns=2:7,digits=3)
  })
  
  ppi <- eventReactive(input$pro_ppi_plot, {
    info <- string_db$get_interactions(pro_ppi_table()$STRING_id)
    info <- info %>% distinct(from, to, .keep_all = T)
    ##igraph
    #Convert stringID to Symbol
    links <- info %>%
      mutate(from = pro_ppi_table()[match(from, pro_ppi_table()$STRING_id), "Protein"]) %>% 
      mutate(to = pro_ppi_table()[match(to, pro_ppi_table()$STRING_id), "Protein"]) %>%  
      dplyr::select(from, to , last_col()) %>% 
      dplyr::rename(weight = combined_score)
    
    # Remove free interactions
    # If a link from of the links data frame appears only once, and the to appears only once, it is removed
    links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1,2,3)
    links_2
  })
  
  output$pro_ppi_link_download <- downloadHandler(
    filename = function() {
      paste0(input$pro_ppi_dataset,"_PPI_Links.txt")
    },
    content = function (file) {
      write.table(ppi(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })

  pho_plot <-  eventReactive(input$pho_action, {
    phospho <- read.table(paste0("pho_",input$pho_dataset,".txt"),header = T,sep = "\t",check.names = F)
    pho_class <- read.table(paste0("patient_",input$pho_dataset,".txt"),header = T,sep = "\t")
    pho_input<- as.data.frame(t(phospho[input$pho_site,]))
    pho_input$tumor <- pho_class$class
    pho_input$tumor <- factor(pho_input$tumor,levels = c("Tumor","Normal"))
    boxplot <- ggboxplot(pho_input, x = "tumor", y = input$pho_site, color = "tumor",width = 0.4,
                         add = "jitter",size = 1.2,palette =c("Tumor"=input$pho_col_tumor,"Normal" =input$pho_col_normal),add.params=list(size=input$pho_df_size,jitter = 0.1))+
      stat_compare_means(method = input$pho_df_method,hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pho_dataset)+xlab(NULL)+ ylab(paste0(input$pho_site," Abundance\n(Log2 TMT ratio)"))+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    boxplot
    
  })
  
  output$pho_plot <- renderPlot({ 
    input$pho_action
    pho_plot()
  })

  output$pho_downloadPdf <- downloadHandler(
    filename = function() {
      paste0("phospho_",input$pho_dataset,"_",input$pho_site,".pdf")
    },
    content = function (file) {
      ggsave(file, plot = pho_plot(),he=8,wi=6)
    })
	
  pho_volcano_plot <- eventReactive(input$pho_volcano_action, {
    volcanoData <- read.csv(paste0("DPS.",input$pho_volcano_dataset,".csv"),row.names = 1)
    markerSite <- volcanoData[rownames(volcanoData)== input$pho_volcano_site,]
    markerSite$labelSite <- rownames(markerSite)
    upSite <- rownames(volcanoData)[which(volcanoData$adj.P.Val< input$pho_volcano_FDR & volcanoData$logFC>input$pho_volcano_log2FC)]
    downSite <- rownames(volcanoData)[which(volcanoData$adj.P.Val< input$pho_volcano_FDR & volcanoData$logFC< -(input$pho_volcano_log2FC))]
    
    normdata <- volcanoData[-which(rownames(volcanoData) %in% c(upSite,downSite)),]
    updata <- volcanoData[which(rownames(volcanoData) %in% upSite),]
    downdata <- volcanoData[which(rownames(volcanoData) %in% downSite),]
    
    p1 <- ggplot()+
      geom_point(pch=20,cex=0.6,color="gray77",
                 data=normdata,aes(x=logFC,y=(-log10(adj.P.Val))))+
      geom_point(pch=20,cex=0.6,color="#FC4E07",
                 data=updata,aes(x=logFC,y=(-log10(adj.P.Val))))+
      geom_point(pch=20,cex=0.6,color="#00AFBB",
                 data=downdata,aes(x=logFC,y=(-log10(adj.P.Val))))+
      geom_point(pch=20,cex=3,color="red",
                 data=markerSite,aes(x=logFC,y=(-log10(adj.P.Val))))+
      theme_classic()+
      geom_vline(xintercept = c(-(input$pho_volcano_log2FC),input$pho_volcano_log2FC),linetype="dotted",size = 0.7)+
      geom_hline(yintercept = c(-log10(input$pho_volcano_FDR)),linetype="dotted",size = 0.7)+
      xlab("log2 (fold change)")+ylab("-log10 (FDR)")

    p1
  })
  
  output$pho_volcano_plot <- renderPlot({ 
    input$pho_volcano_action
    pho_volcano_plot()
  })
  
  output$pho_volcano_plot_download <- downloadHandler(
    filename = function() {
      paste0(input$pho_volcano_dataset,"_",input$pho_volcano_site,"_volcano.pdf")
    },
	
    content = function (file) {
      ggsave(file, plot = pho_volcano_plot(),he=6,wi=6)
    })
	
  pho_volcano_table <- eventReactive(input$pho_volcano_action, {
    DPS <- read.csv(paste0("DPS.",input$pho_volcano_dataset,".csv"),row.names = 1)
    DPS <- DPS[order(abs(DPS$logFC),decreasing = T),]
    DPS
  })

  output$pho_volcano_table <- DT::renderDataTable({
    input$pho_volcano_action
    formatSignif(datatable(pho_volcano_table()),columns=1:6,digits=3)
  })

  output$pho_volcano_csv_download <- downloadHandler(
    filename = function() {
      paste0(input$pho_volcano_dataset,"_phosphosite_volcano.csv")
    },
    content = function (file) {
      write.csv(pho_volcano_table(),file,row.names = T,quote = F)
    })
	
  
  pho_corr_plot <- eventReactive(input$pho_action2, {
    protein_tumor_nonImpute <- read.table(paste0(input$pho_corr_dataset,"/","protein_Tumor_",input$pho_corr_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    protein_tumor_nonImpute <- as.data.frame(t(protein_tumor_nonImpute[input$pro_A,]))
    pho_tumor_nonImpute <- read.table(paste0(input$pho_corr_dataset,"/","pho_Tumor_",input$pho_corr_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    pho_tumor_nonImpute <- as.data.frame(t(pho_tumor_nonImpute[input$pho_A,]))

    pro_pho_nonImupute <- merge(protein_tumor_nonImpute,pho_tumor_nonImpute,by="row.names")
    pro_pho_nonImupute <- na.omit(pro_pho_nonImupute)
    
    pho_corr_nonImpute <-ggscatter(pro_pho_nonImupute,x = input$pro_A, y = input$pho_A,
                               conf.int = TRUE,size = input$pho_corr_size,color = input$pho_col_non_imputed_corr,title = "Filter missing value",
                               xlab = paste0(input$pro_A," Abundance\n(Log2 TMT ratio)"),
                               ylab = paste0(input$pho_A," Abundance\n(Log2 TMT ratio)"))+
      stat_cor(method = input$pho_correlation,size=5)+geom_smooth(method = "lm",linewidth=1.2,colour="black")+
      theme(plot.title = element_text(hjust = 0.5))
    
    protein_tumor_impute <- read.table(paste0(input$pho_corr_dataset,"/","protein_Tumor_",input$pho_corr_dataset,"_Impute.txt"),header = T,sep = "\t",check.names=F)
    protein_tumor_impute <-  as.data.frame(t(protein_tumor_impute[input$pro_A,]))
    pho_tumor_impute <- read.table(paste0(input$pho_corr_dataset,"/","pho_Tumor_",input$pho_corr_dataset,"_Impute.txt"),header = T,sep = "\t",check.names=F)
    pho_tumor_impute <- as.data.frame(t(pho_tumor_impute[input$pho_A,]))
    
    pro_pho_Imupute <- merge(protein_tumor_impute,pho_tumor_impute,by="row.names")
    rownames(pro_pho_Imupute) <- pro_pho_Imupute$Row.names
    pro_pho_Imupute <- pro_pho_Imupute[,-1]
    pro_pho_Imupute$impute <- ifelse(rownames(pro_pho_Imupute)%in%pro_pho_nonImupute$Row.names,"Non_imputed","Imputed")
    pro_pho_Imupute$impute <- as.factor(pro_pho_Imupute$impute) 
    
    pho_corr_impute <- ggscatter(pro_pho_Imupute, x = input$pro_A, y = input$pho_A,
                             conf.int = TRUE,size = input$pho_corr_size,color  = "impute",palette = c(Imputed=input$pho_col_imputed_corr, Non_imputed=input$pho_col_non_imputed_corr),
                             title = "Impute missing value",show.legend.text = T,
                             xlab = paste0(input$pro_A," Abundance\n(Log2 TMT ratio)"),
                             ylab = paste0(input$pho_A," Abundance\n(Log2 TMT ratio)"))+
      stat_cor(method = input$pho_correlation,size=5)+geom_smooth(method = "lm",linewidth=1.2,colour="black")+
      theme(plot.title = element_text(hjust = 0.5))
    
    arrange <- ggarrange(pho_corr_nonImpute,pho_corr_impute,  nrow = 1, ncol = 2,legend = "right",common.legend = T)
    arrange
  })

  output$pho_corr_plot <- renderPlot({
    input$pho_action2
    pho_corr_plot()
  })

  output$pho_downloadPdf2 <- downloadHandler(
    filename = function() {
      paste0(input$pho_corr_dataset,"_",input$pro_A,"_",paste(input$pro_A,input$pho_A,sep = "_"),"_corr",".pdf")
    },
    content = function (file) {
      ggsave(file, plot = pho_corr_plot(),he=6,wi=12)
    })

  ks_cor_impute <- eventReactive(input$pho_action3, {
    protein_tumor_impute <- read.table(paste0(input$pho_ks_corr_dataset,"/","protein_Tumor_",input$pho_ks_corr_dataset,"_Impute.txt"),header = T,sep = "\t",check.names=F)
    protein_tumor_impute <- as.data.frame(t(protein_tumor_impute[c(input$kinase),]))
    pho_tumor_impute <- read.table(paste0(input$pho_ks_corr_dataset,"/","pho_Tumor_",input$pho_ks_corr_dataset,"_Impute.txt"),header = T,sep = "\t",check.names=F)
    pho_imputed_fraction <- read.csv(paste0(input$pho_ks_corr_dataset,"/phospho_imputed_fraction_",input$pho_ks_corr_dataset,".csv"),row.names = 1,check.names = F)
    ks <- merge(protein_tumor_impute,as.data.frame(t(pho_tumor_impute)),by = "row.names")

    ks_cor_impute <- do.call(rbind, lapply(3:ncol(ks), function(x){
      dd = cor.test(ks[,input$kinase], ks[,x], method = input$pho_ks_correlation )
      data.frame(kinase = input$kinase, substrate = colnames(ks)[x], R = signif(dd$estimate,4), p.value = signif(dd$p.value,4))
    }))
    ks_cor_impute$Imputed_fraction <- pho_imputed_fraction$`Imputed_fraction(%)`
    ks_cor_impute$adjP.value <- signif(p.adjust(ks_cor_impute$p.value,method = "BH"),digits = 4)
    ks_cor_impute <- ks_cor_impute[order(ks_cor_impute$R,decreasing = T),]
    ks_cor_impute <- ks_cor_impute[,c(1,2,5,3,4,6)]
    colnames(ks_cor_impute) <- c("Kinase","Substrates","Imputed fraction(%)","Correlation coefficient","P.value","adjP.value")
    ks_cor_impute
  })

  output$ks_table_impute <- DT::renderDataTable({
    input$pho_action3
    DT::datatable(ks_cor_impute(),rownames = F)
  })
  
  output$pho_ks_table_impute_download <- downloadHandler(
    filename = function() {
      paste0(input$pho_ks_corr_dataset,"_",input$kinase,"_impute_corr",".txt")
    },
    content = function (file) {
      write.table(ks_cor_impute(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })

  ks_cor_nonImpute <- eventReactive(input$pho_action3, {
    protein_tumor_nonImpute <- read.table(paste0(input$pho_ks_corr_dataset,"/","protein_Tumor_",input$pho_ks_corr_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    protein_tumor_nonImpute <- as.data.frame(t(protein_tumor_nonImpute[input$kinase,]))
    pho_tumor_nonImpute <- read.table(paste0(input$pho_ks_corr_dataset,"/","pho_Tumor_",input$pho_ks_corr_dataset,"_nonImpute.txt"),header = T,sep = "\t",check.names=F)
    pho_tumor_nonImpute <- as.data.frame(t(pho_tumor_nonImpute))
    ks <- merge(protein_tumor_nonImpute,pho_tumor_nonImpute,by = "row.names")

    ks_cor_nonImpute <- do.call(rbind, lapply(3:ncol(ks), function(x){
      dd = cor.test(ks[,input$kinase], ks[,x], method = input$pho_ks_correlation,use="pairwise.complete.obs" )
      data.frame(kinase = input$kinase, substrate = colnames(ks)[x], R = signif(dd$estimate,4), p.value = signif(dd$p.value,4))
    }))
    ks_cor_nonImpute$adjP.value <- signif(p.adjust(ks_cor_nonImpute$p.value,method = "BH"),digits = 4)
    ks_cor_nonImpute <- ks_cor_nonImpute[order(ks_cor_nonImpute$R,decreasing = T),]
    colnames(ks_cor_nonImpute) <- c("Kinase","Substrates","Correlation coefficient","P.value","adjP.value")
    ks_cor_nonImpute
  })

  output$ks_table_nonImpute <- DT::renderDataTable({
    input$pho_action3
    DT::datatable(ks_cor_nonImpute(),rownames = F)
  })
  output$pho_ks_table_nonImpute_download <- downloadHandler(
    filename = function() {
      paste0(input$pho_ks_corr_dataset,"_",input$kinase,"_nonImpute_corr",".txt")
    },
    content = function (file) {
      write.table(ks_cor_nonImpute(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })

  PX <- eventReactive(input$pho_action4, {
    PX <- read.csv(paste0(input$KSEA_dataset,"/","KSEA.",input$KSEA_dataset,".csv"),header = T)
    PX        
  })
  choice <- eventReactive(input$pho_action4, {
    switch (input$PhosphoSite_dataset,
            "no" = FALSE,
            "yes" = TRUE
    )
  })

  KSEA.Bar <- eventReactive(input$pho_action4, {
    KSEA.Bar <- KSEA.Barplot(KSData,
                 PX(),
                 NetworKIN = choice(),
                 NetworKIN.cutoff = input$NetworKIN.cutoff,
                 m.cutoff = input$m.cutoff,
                 p.cutoff = input$p.cutoff,
                 export = F)
    KSEA.Bar
   
  })
    
  output$ksea_barplot <- renderPlot({
    input$pho_action4
    KSEA.Bar()
  },
  height = eventReactive(input$pho_action4, {
    height <- KSEA.Bar.height()*100
    return(height)
  })
  )

  output$KSEA_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$KSEA_dataset,"_KSEA_barplot.pdf")
    },
    content = function (file) {
      pdf(file,width = 10,        
          height = 2*KSEA.Bar.height())
      par(mai=c(1,1,.4,.4))
      KSEA.Barplot(KSData,
                   PX(),
                   NetworKIN = choice(),
                   NetworKIN.cutoff = input$NetworKIN.cutoff,
                   m.cutoff = input$m.cutoff,
                   p.cutoff = input$p.cutoff,
                   export = F)
      dev.off()
    })
  
  KSEA_scores <- eventReactive(input$pho_action4, {
    KSEA_scores <- KSEA.Scores(KSData, PX(), NetworKIN = choice(), NetworKIN.cutoff = input$NetworKIN.cutoff)
    KSEA_scores
  })
  
  output$ksea_scores <- DT::renderDataTable({
    formatSignif(DT::datatable(KSEA_scores(),rownames = F),columns=c(2,3,5,6,7),digits=3)
  })

  output$KSEA_scores_download <- downloadHandler(
    filename = function() {
      paste0(input$KSEA_dataset," Kinase Scores.csv")
    },
    content = function (file) {
      write.csv(KSEA_scores(),file,row.names = F,quote = F)
    })
  
  KSEA_links <- eventReactive(input$pho_action4, {
    KSEA_links <- KSEA.KS_table(KSData, PX(), NetworKIN = choice(), NetworKIN.cutoff = input$NetworKIN.cutoff)
    row.names(KSEA_links) <- c(1:nrow(KSEA_links))
    KSEA_links
  })
  output$ksea_links <- DT::renderDataTable({
    formatSignif(DT::datatable(KSEA_links(),rownames = F),columns=c(5),digits=3)
  })

  output$KSEA_links_download <- downloadHandler(
    filename = function() {
      paste0(input$KSEA_dataset," Kinase-Substrate Links.csv")
    },
    content = function (file) {
      write.csv(KSEA_links(),file,row.names = F,quote = F)
    })


  pho_survival_plot <- eventReactive(input$pho_action5, {
    pho_tumor <- read.table(paste0(input$pho_survival_plot_dataset,"/","pho_Tumor_",input$pho_survival_plot_dataset,".txt"),header = T,sep = "\t",check.names=F)
    clinical <- read.table(paste0(input$pho_survival_plot_dataset,"/","clinical_",input$pho_survival_plot_dataset,".txt"),header = T,sep = "\t")
    clinical$site <- as.numeric(pho_tumor[input$pho_survival_plot_site,])
    clinical <- clinical[!is.na(clinical$OS_days),]
    if (input$pho_cutoff == "Median") {
      clinical$Group <- ifelse(clinical$site >= median(clinical$site,na.rm = T),"High","Low")
      lab1 <- paste0(names(table(clinical$Group)[1])," (N=",table(clinical$Group)[1],")")
      lab2 <- paste0(names(table(clinical$Group)[2])," (N=",table(clinical$Group)[2],")")
      fit <- survfit(Surv(OS_days,Status) ~ Group, 
                     data = clinical)
      survival <- ggsurvplot(fit, data = clinical,palette = c(input$pho_col_high,input$pho_col_low),pval = TRUE, conf.int = T,conf.int.style = "step",
                             xlab = "Follow up time (days)",ylab = "Overall survival rate",title = input$pho_survival_plot_site,
                             legend.title = "",
                             legend = c(0.3,0.18),pval.size=5.2,
                             legend.labs= c(lab1,lab2))%++%
        theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
          legend.background = element_blank())+
        guides(color=guide_legend(override.aes = list(shape = NA)))
      
    } else if (input$pho_cutoff == "Optimal value") {
      res.cut <- surv_cutpoint(clinical, 
                               time = "OS_days", 
                               event = "Status", 
                               variables = c("site"), 
                               minprop = 0.3
      )
        res.cat <- surv_categorize(res.cut,labels = c("Low", "High"))
        lab_high <- paste0(names(table(res.cat$site)[1])," (N=",table(res.cat$site)[1],")")
        lab_low <- paste0(names(table(res.cat$site)[2])," (N=",table(res.cat$site)[2],")")
        fit2 <- survfit(Surv(OS_days, Status) ~site, data = res.cat)
        data.optimal <- res.cat

      survival <- ggsurvplot(fit2, data = data.optimal,palette = c(input$pho_col_high,input$pho_col_low),pval = TRUE, conf.int = T,conf.int.style = "step",
                             xlab = "Follow up time (days)",ylab = "Overall survival rate",title = input$pho_survival_plot_site,
                             legend.title = "",
                             legend = c(0.3,0.18),pval.size=5.2,
                             legend.labs= c(lab_high,lab_low))%++%
        theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
              legend.background = element_blank())+
        guides(color=guide_legend(override.aes = list(shape = NA)))
    }
    survival
  })
  
  output$pho_survival_plot <- renderPlot({
    input$pho_action5
    pho_survival_plot()
  })

  output$pho_downloadPdf5 <- downloadHandler(
    filename = function() {
      paste0(input$pho_survival_plot_dataset,"_",input$pho_survival_plot_site,"_survival",".pdf")
    },
    content = function (file) {
      pho_survival_plot()
      ggsave(file,he=7,wi=7)
    })
	
  pho_survival_table <- eventReactive(input$pho_survival_table_action, {
    pho_tumor <- read.table(paste0(input$pho_survival_table_dataset,"/","pho_Tumor_",input$pho_survival_table_dataset,".txt"),header = T,sep = "\t",check.names=F)
    clinical <- read.table(paste0(input$pho_survival_table_dataset,"/","clinical_",input$pho_survival_table_dataset,".txt"),header = T,sep = "\t")
    clinical <- clinical[!is.na(clinical$OS_days),]
    
    gene_all <- str_split(rownames(pho_tumor),"_",simplify = T)[,1]
    gene_input <- as.data.frame(t(pho_tumor[which(gene_all%in%input$pho_survival_table_protein),]))
    survival.site <- merge(clinical[,c("Patient","OS_days","Status")],gene_input,by.x = "Patient",by.y = "row.names")
    
	##Calculate the median value of each phosphosite in all samples and assign it to a
    a <- apply(survival.site[,-c(1:3),drop=F], 2, median)
    a=data.frame(a)
    
	##According to the median value, the samples were divided into high and low phosphorylation level groups
    for (i in 4:ncol(survival.site)) {
      n=colnames(survival.site)[i]
      survival.site[,n] = ifelse(survival.site[,n]>= a[n,],"High","Low")
    }
    
    survival.p <- c()
    for (s in 4:ncol(survival.site)){
      surv_diff <- survdiff(Surv(OS_days, Status)~survival.site[,s],data = survival.site)
      p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
      survival.p <- c(survival.p,p.value)
    }
    p.adj <- signif(p.adjust(survival.p,method = "BH"),digits = 4)
    survival.result <- data.frame(Site=colnames(survival.site)[4:ncol(survival.site)],pvalue=signif(survival.p,digits = 4),p.adjust = p.adj)
    survival.result
  })

  output$pho_survival_table <- DT::renderDataTable({
    input$pho_survival_table_action
    DT::datatable(pho_survival_table(),rownames = F)
  })

  output$pho_survival_table_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pho_survival_table_dataset,"_",input$pho_survival_table_protein,"_survial",".txt")
    },
    content = function (file) {
      write.table(pho_survival_table(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    })
	
  pho_age_plot <-  eventReactive(input$pho_age_action, {
    phospho <- read.table(paste0("pho_Tumor_",input$pho_age_dataset,".txt"),header = T,sep = "\t",check.names = F)
    clinical <- read.table(paste0("clinical_",input$pho_age_dataset,".txt"),header = T,sep = "\t")
    site_input<- as.data.frame(t(phospho[input$pho_age_site,]))
    site_input$age <- clinical$Age
    site_input <- na.omit(site_input)
    
    if (input$pho_age_cutoff == "median") {
      cutoff <- median(site_input$age)
    }else if (input$pho_age_cutoff == "custom") {
      cutoff <- input$pho_age_custom
    }
    site_input$Group <- ifelse(site_input$age >= cutoff,"Older","Younger")
    site_input$Group <- factor(site_input$Group,levels = c("Younger","Older"))
    
    violinplot <- ggviolin(site_input, x = "Group", y = input$pho_age_site,fill = "Group",width = 0.6,
                           add = "boxplot",size = 0.8,palette =c("Younger"=input$pho_age_col_young,"Older" =input$pho_age_col_old),add.params = list(width=0.1,size=0.7,fill="white")) +
      stat_compare_means(method = input$pho_age_method,hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pho_age_dataset)+xlab(NULL)+scale_x_discrete(labels = c(paste0("Younger\n(<",cutoff,")"),paste0("Older\n(>=",cutoff,")")))+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    violinplot
  })
  
  output$pho_age_plot <- renderPlot({
    input$pho_age_action
    pho_age_plot()
  })

  output$pho_age_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pho_age_dataset,"_",input$pho_age_site,"_Age.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pho_age_plot(),he=8,wi=6)
    })
	
  pho_gender_plot <-  eventReactive(input$pho_gender_action, {
    phospho <- read.table(paste0("pho_Tumor_",input$pho_gender_dataset,".txt"),header = T,sep = "\t",check.names = F)
    clinical <- read.table(paste0("clinical_",input$pho_gender_dataset,".txt"),header = T,sep = "\t")
    site_input<- as.data.frame(t(phospho[input$pho_gender_site,]))
    site_input$Gender <- clinical$Gender
    site_input <- na.omit(site_input)
    site_input$Gender <- factor(site_input$Gender)
    
    violinplot <- ggviolin(site_input, x = "Gender", y = input$pho_gender_site,fill = "Gender",width = 0.6,
                           add = "boxplot",size = 0.8,palette =c("Female"=input$pho_gender_col_female,"Male" =input$pho_gender_col_male),add.params = list(width=0.1,size=0.7,fill="white")) +
      stat_compare_means(method = input$pho_gender_method,hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pho_gender_dataset)+xlab(NULL)+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    violinplot
  })
  
  output$pho_gender_plot <- renderPlot({
    input$pho_gender_action
    pho_gender_plot()
  })

  output$pho_gender_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pho_gender_dataset,"_",input$pho_gender_site,"_Gender.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pho_gender_plot(),he=8,wi=6)
    })

  pho_stage_plot <-  eventReactive(input$pho_stage_action, {
    phospho <- read.table(paste0("pho_Tumor_",input$pho_stage_dataset,".txt"),header = T,sep = "\t",check.names = F)
    clinical <- read.table(paste0("clinical_",input$pho_stage_dataset,".txt"),header = T,sep = "\t")
    site_input<- as.data.frame(t(phospho[input$pho_stage_site,]))
    
    site_input$Stage <- clinical$Tumor_stage
    site_input <- na.omit(site_input)
    site_input$Stage <- factor(site_input$Stage)
    
    violinplot <- ggviolin(site_input, x = "Stage", y = input$pho_stage_site,fill = "Stage",width = 0.7,palette = "RdYIGn",
                           add = "boxplot",size = 0.8,add.params = list(width=0.1,size=0.7,fill="white")) +
      stat_compare_means(method = "anova",hjust = 0.5,label.x.npc = "center",vjust = 0,size = 5)+
      labs(title = input$pho_stage_dataset)+xlab(NULL)+
      theme(legend.position="none",plot.title=element_text(hjust=0.5))
    violinplot
  })
  
  output$pho_stage_plot <- renderPlot({
    input$pho_stage_action
    pho_stage_plot()
  })

  output$pho_stage_downloadPdf <- downloadHandler(
    filename = function() {
      paste0(input$pho_stage_dataset,"_",input$pho_stage_site,"_Stage.pdf")
    },
    content = function (file) {
      ggsave(file, plot = pho_stage_plot(),he=8,wi=12)
    })
  
})