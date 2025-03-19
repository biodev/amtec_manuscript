get.metabric.exprs <- function(metabric.exprs.file){
    
    ma.dt <- fread(metabric.exprs.file)
    
    ma.dt <- ma.dt[is.na(Entrez_Gene_Id)==F]
    
    avg.exprs <- rowMeans(as.matrix(ma.dt[,3:ncol(ma.dt)]), na.rm=T)
    
    ma.dt[,avg_exprs:=avg.exprs]
    
    ma.dt <- ma.dt[order(-avg_exprs),]
    
    uniq.dt <- ma.dt[!duplicated(Entrez_Gene_Id),]
    
    uniq.dt <- uniq.dt[is.na(Hugo_Symbol)==F,]
    
    stopifnot(uniq.dt[,.N,by=Hugo_Symbol][,all(N == 1)])
    
    ma.mat <- as.matrix(uniq.dt[,-c("Hugo_Symbol", "Entrez_Gene_Id", "avg_exprs"),with=F])
    rownames(ma.mat) <- uniq.dt$Hugo_Symbol
    
    ma.mat <- ma.mat[complete.cases(ma.mat),]
    
    ma.mat <- ma.mat[matrixStats::rowMads(ma.mat) > 0,]
    
    ma.mat
}

get.metabric.clin <- function(metabric.pat.file, metabric.sample.file, ma.mat){
    
    pat.clin <- fread(metabric.pat.file, header=T, skip=4)
    
    pat.clin[,RFS:=as.numeric(sapply(strsplit(RFS_STATUS, "\\:"), function(x) ifelse(length(x) == 0, NA, x[1])))]
    
    pat.clin[,OS:=as.numeric(sapply(strsplit(OS_STATUS, "\\:"), function(x) ifelse(length(x) == 0, NA, x[1])))]
    
    samp.clin <- fread(metabric.sample.file,  skip=4)
    
    stopifnot(samp.clin[,.N,by=PATIENT_ID][,all(N==1)])
    
    pat.basal <- pat.clin[CLAUDIN_SUBTYPE=="Basal"]
    
    pat.basal <- merge(pat.basal, samp.clin[,.(PATIENT_ID, GRADE, TUMOR_SIZE)], by="PATIENT_ID")
    
    pat.basal <- pat.basal[PATIENT_ID %in% colnames(ma.mat),]
    
    pat.basal[,r_acc:=PATIENT_ID]
    
    pat.basal
}

subset.exprs <- function(abund, clin){
    
    ma.mat <- abund[,clin$PATIENT_ID]
    
    ma.mat <- ma.mat[matrixStats::rowMads(ma.mat) > 0,]
    
    ma.mat
}

.score.genes.ssgsea <- function(ma.mat, geneset, gene.universe){
    
    use.mm <- gene.universe[gene_id %in% rownames(ma.mat)]
    
    ma.mat <- ma.mat[use.mm$gene_id,]
    
    mod.list.symb <- lapply(split(geneset, by="signature"), function(x) x$gene_id)
    
    mod.enrich <- gsva(ma.mat,mod.list.symb , method="ssgsea", ssgsea.norm=F)
    
    mod.scores <- data.table(reshape2::melt(mod.enrich, as.is=T))
    names(mod.scores) <- c("signature", "r_acc", "value")
    
    mod.scores
}

score.m3s <- function(mm, m3s, ma.mat, signature){
    
    m3s$signature <- signature
    
    m3s.dt <- .score.genes.ssgsea(ma.mat, m3s, mm)
    
    m3s.dt
}

compute.zscore <- function(gl, abund){
    
    use.gl <- gl[gene_id %in% rownames(abund),]
    
    gl.list <- split(use.gl$gene_id, use.gl$signature)
    
    tmp.abund <- t(scale(t(abund[unique(use.gl$gene_id),,drop=F])))
    
    tmp.res <- rbindlist(lapply(gl.list, function(y){
        
        numer <- colSums(tmp.abund[y,,drop=F])
        as.data.table(numer / sqrt(length(y)), keep.rownames=T)
        
    }), idcol="signature")
    
    tmp.res[,.(signature, r_acc=V1, value=V2)]
    
}

combine.metabric.preds <- function(meta.clin, meta.ssg, meta.mzb1, ssg.preds, mzb1.preds){
    
    meta.ssg[,value:=(value-min(value))/(max(value)-min(value))]
    
    meta.ssg$preds <- predict(ssg.preds$`Bx1.Mod3 (brown)`$fit, newdata=meta.ssg)
    
    meta.mzb1$preds <- predict(mzb1.preds$Bx1$fit, newdata=meta.mzb1)
    
    comb.preds <- merge(dcast(r_acc~signature, value.var=c("value", "preds"),data=meta.ssg),
                        dcast(r_acc~signature, value.var=c("value", "preds"),data=meta.mzb1), 
                        by="r_acc")
    
    meta.clin <- merge(meta.clin, comb.preds, by="r_acc")
    
    meta.clin
}

fit.metabric.full.models <- function(meta.clin){
    
    base.model <- coxph(Surv(time=OS_MONTHS, event=OS, type="right")~1+pspline(AGE_AT_DIAGNOSIS) + pspline(LYMPH_NODES_EXAMINED_POSITIVE) + GRADE + pspline(TUMOR_SIZE) + strata(as.character(COHORT)), 
                        data=as.data.frame(meta.clin))
    
    use.coefs <- names(meta.clin)[grepl("value_", names(meta.clin)) | grepl("preds_", names(meta.clin))]
    
    rbindlist(lapply(use.coefs, function(x){
        
        tmp.form <- paste0('Surv(time=OS_MONTHS, event=OS, type="right")~', x, '+pspline(AGE_AT_DIAGNOSIS) + pspline(LYMPH_NODES_EXAMINED_POSITIVE) + GRADE + pspline(TUMOR_SIZE) + strata(as.character(COHORT))')
        
        tmp.mod <- coxph(as.formula(tmp.form), data=as.data.frame(meta.clin))
        
        lrt <- as.data.table(broom::tidy(anova(base.model, tmp.mod)))[2,.(term=ifelse(grepl("preds", x), paste0(x, "SDPR"), x), lrt_logLik=logLik, lrt_statistic=statistic, lrt_df=df, lrt_p.value=p.value)]
        
        coef.dt <- merge(as.data.table(broom::tidy(tmp.mod)), lrt, by="term", all.x=T, all.y=F)
        
        coef.dt$signature <- x
        
        coef.dt
    }))
    
}

summarize.metabric.univ <- function(meta.preds){
    
    preds.names <- names(meta.preds)[grepl("preds_", names(meta.preds))]
    
    rbindlist(lapply(setNames(preds.names, preds.names), function(x){
        
        use.form <- paste('Surv(time=OS_MONTHS, event=OS, type="right")~', x)
        
        as.data.table(broom::glance(survdiff(as.formula(use.form), data=as.data.frame(meta.preds))))
        
    }), idcol="term")
    
}

plot.metabric.km <- function(meta.preds){
    
    meta.preds[,mTNBC3s_top3:=factor(ifelse(preds_mTNBC3s_Top3 == "PD", "Low mTNBC3s_top3", "High mTNBC3s_top3"))]
    
    basal.p <- ggsurvfit(survfit2(Surv(time=OS_MONTHS/12, event=OS, type="right")~mTNBC3s_top3, data=as.data.frame(meta.preds)), linewidth=1) +
        add_censor_mark(size = 2, alpha = 0.75) +
        scale_color_manual(values=setNames(plasma(n=5)[c(1, 4)], c("Low mTNBC3s_top3", "High mTNBC3s_top3"))) +
        scale_fill_manual(values=setNames(plasma(n=5)[c(1, 4)], c("Low mTNBC3s_top3", "High mTNBC3s_top3"))) +
        #add_confidence_interval() +
        ggtitle("Overall Survival Basal Classified Metabric Cohort") + xlab("Years") +
        add_risktable() +
        guides(color=guide_legend(position = "bottom"),
               fill=guide_legend(position = "bottom")) +
        theme_bw(base_size=10)
    
    ggsave(basal.p, file="figures/brown_ssgsea_metabric-new.pdf", width=4, height=4)
    
    "figures/brown_ssgsea_metabric-new.pdf"
}

write.full.models <- function(meta.res, meta.univ, calgb.res=NULL){
    
    openxlsx::write.xlsx(list(metabric_full_models=meta.res, metabric_univ_results=meta.univ, calgb_univ_results=calgb.res), file="output/metabric_calgb_full_model_summary.xlsx")
    
    "output/metabric_calgb_full_model_summary.xlsx"
}