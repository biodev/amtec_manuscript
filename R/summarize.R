fit.blia.blis.tcga <- function(comb.mes, dtype.rna, batch.bl){
    
    comb.mes <- merge(comb.mes[,.(cur_label, r_acc=ptid, PC1, PC2)], dtype.rna, by="r_acc")
    
    batch.bl <- batch.bl[final_call %in% c("BLIA", "BLIS")]
    
    tcga.mes.bl <- merge(batch.bl, comb.mes, by="r_acc")
    
    tcga.mat <- dcast(ptid+final_call~make.names(cur_label), value.var="PC1",data=tcga.mes.bl)
    tcga.mat[,fc_fac:=factor(final_call)]
    
    tcga.mat.tbl <- as_tibble(tcga.mat)
    
    mod.bb.accs <- bind_rows(lapply(names(tcga.mat.tbl)[grepl("Mod\\d+", names(tcga.mat.tbl))], function(x){
        
        use.form <- as.formula(paste("fc_fac", "~", x))
        
        bb_fit <- 
            logistic_reg() %>% 
            set_engine("glm") %>%
            fit(use.form, data = tcga.mat.tbl)
        
        bb_preds <- augment(
            bb_fit, 
            new_data=tcga.mat.tbl)
        
        summarize(bb_preds, 
                  acc=accuracy_vec(truth=fc_fac, estimate=.pred_class),
                  mcc=mcc_vec(truth=fc_fac, estimate=.pred_class),
                  f1=f_meas_vec(truth=fc_fac, estimate=.pred_class)
        ) %>%
            mutate(Module=x) %>%
            select(Module, acc, mcc, f1)
        
    }))
    
    openxlsx::write.xlsx(list(Pred_Accuracy=arrange(mod.bb.accs, -acc)), file="output/module_burstein_preds.xlsx")
    
    "output/module_burstein_preds.xlsx"
    
}

enrich.modules <- function(kme.annot, geneset.path){
    
    gmt.df <- clusterProfiler::read.gmt(geneset.path)
    
    mod.enrich <- rbindlist(lapply(split(kme.annot, by="cur_label"), function(x){
        
        if(length(intersect(x$gene_name, gmt.df$gene)) == 0){
            
            NULL
            
        }else{
            
            tmp <- clusterProfiler::enricher(
                gene=x$gene_name,
                pvalueCutoff = 1,
                pAdjustMethod = "none",
                universe=kme.annot$gene_name,
                minGSSize = 5,
                maxGSSize = 1000,
                qvalueCutoff = 1,
                TERM2GENE=gmt.df,
                TERM2NAME = NA
            )
            
            tmp.dt <- as.data.table(tmp)
            
            tmp.dt
            
        }
        
        
    }), idcol="cur_label")
    
    mod.enrich <- mod.enrich[cur_label != "Mod0 (grey)"]
    
    mod.enrich[,padj:=p.adjust(pvalue, method="BY")]
    
    mod.enrich
    
}

.fit.ctree <- function(x, rescale.value=F){
    
    if(rescale.value == T){
        x[,value:=(value-min(value))/(max(value)-min(value))]
    }
    
    tmp <- ctree(bin_resp~., data=as.data.frame(x[,c("bin_resp", "value"),with=F]), control=ctree_control( testtype ="Univariate", mincriterion = 0.95 ,minbucket=3, minsplit=3, maxdepth=1))
    
    tmp.split <- split_node(node_party(tmp[[1]]))
    
    if (is.null(tmp.split$breaks)){
        tmp.dt <- data.table(split=NA, pval=node_party(tmp[[1]])$info$p.value)
    }else{
        tmp.dt <- data.table(split=tmp.split$breaks, pval=node_party(tmp[[1]])$info$p.value)
    }
    
    tmp.acc <- caret::confusionMatrix(data = predict(tmp), reference = x$bin_resp)
    
    tmp.mcc <- mcc_vec(truth=x$bin_resp, estimate=predict(tmp))
    
    tmp.dt <- cbind(tmp.dt, t(tmp.acc$overall), t(tmp.acc$byClass), mcc=tmp.mcc)
    
    list(summary=tmp.dt, fit=tmp, preds=cbind(x, preds=predict(tmp)))
}

.pred.ctree <- function(val, trn, keep.cols, rescale.value=F){
    
    if (rescale.value == T){
        val[,value:=(value-min(value))/(max(value)-min(value))]
    }
    
    tmp.preds <- predict(trn$fit, newdata=as.data.frame(val))
    
    tmp.acc <- caret::confusionMatrix(data = tmp.preds, reference = val$bin_resp)
    
    tmp.mcc <- mcc_vec(truth=val$bin_resp, estimate=tmp.preds)
    
    list(summary=data.table(val[1,keep.cols, with=F], t(tmp.acc$overall),
                            t(tmp.acc$byClass), mcc=tmp.mcc ),
         preds=cbind(val, preds=tmp.preds))
}

classify.valid.mes <- function(comb.mes, mod.res, valid.samps){
    
    comb.mes <- merge(valid.samps, comb.mes,by.x="r_acc" , by.y="ptid")
    
    comb.mes[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    comb.mes[,value:=PC1]
    
    split.mod <- split(comb.mes, by=c("Biopsy", "cur_label"))
    
    res.list <- mapply(function(val, trn){
        
        .pred.ctree(val, trn, keep.cols=c("Biopsy", "cur_label"))
        
    }, split.mod, mod.res[names(split.mod)], SIMPLIFY=F)
    
    res.list
}


ctree.mes <- function(comb.mes, all.rna){
    
    comb.mes <- merge(all.rna, comb.mes, by.x=c("r_acc", "cohort"), by.y=c("ptid",  "cohort"))
    
    comb.mes <- comb.mes[cohort == "AMTEC"]
    
    comb.mes[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    comb.mes[,value:=PC1]
    
    mod.res <- lapply(split(comb.mes, by=c("biopsy", "cur_label")), .fit.ctree)
    
    mod.res
}

ctree.mzb1 <- function(mm, kal.genes, dtype.rna){
    
    mzb1.z <- t(scale(t(log2(kal.genes$kal.tpm + 1))))[mm[gene_name == "MZB1"]$gene_id,]
    
    mzb1.dt <- as.data.table(mzb1.z, keep.rownames=T)
    mzb1.dt$exprs <- log2(kal.genes$kal.tpm + 1)[mm[gene_name == "MZB1"]$gene_id, mzb1.dt$rn]
    names(mzb1.dt) <- c("r_acc", "value", "exprs")
    
    amtec.rna <- merge(dtype.rna, mzb1.dt, by="r_acc")
    
    amtec.rna[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    lapply(split(amtec.rna, by="biopsy"), function(x){
        
        .fit.ctree(x)
    })
}

ctree.ssgsea <- function(mm, kal.genes, dtype.rna, m3s){
    
    mod.list <- lapply(split(m3s, by="cur_label"), function(x) x$gene_id)
    
    kal.abund <- kal.genes$kal.tpm[mm$gene_id,]
    
    amtec.enrich <- suppressWarnings(gsva(kal.abund,mod.list, method="ssgsea", ssgsea.norm=F))
    
    enrich.dt <- data.table(reshape2::melt(amtec.enrich, as.is=T))
    
    enrich.dt <- merge(enrich.dt[,.(cur_label=Var1, r_acc=Var2, value)], dtype.rna, by="r_acc")
    
    enrich.dt[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    lapply(split(enrich.dt, by=c("biopsy", "cur_label")), function(x){
        
        .fit.ctree(x, rescale.value = T)
    })
    
}

classify.valid.mzb1 <- function(mm, val.samps, val.exprs.file, mzb1.res){
    
    load(val.exprs.file)
    
    #valid.mzb1 <- .rank.single.gene(kal.genes$abundance, mm[gene_name == "MZB1"]$gene_id)
    
    valid.mzb1 <- as.data.table(log2(kal.genes$abundance+1)[mm[gene_name == "MZB1"]$gene_id,], keep.rownames=T)
    names(valid.mzb1) <- c("r_acc", "exprs")
    
    valid.mzb1 <- merge(valid.mzb1, val.samps, by="r_acc")
    
    valid.mzb1[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    split.mzb1 <- split(valid.mzb1, by="Biopsy")
    
    res.list <- mapply(function(val, trn){
        
        tmp.val <- copy(val)
        
        tmp.val[,value:=(exprs - mean(trn$preds$exprs) / sd(trn$preds$exprs))]
        
        .pred.ctree(tmp.val, trn, keep.cols=c("Biopsy"))
        
    }, split.mzb1, mzb1.res[names(split.mzb1)], SIMPLIFY=F)
    
    res.list
    
}

classify.valid.ssgsea <- function(mm, val.samps, val.exprs.file, ssg.res, m3s){
    
    load(val.exprs.file)
    
    mod.list <- lapply(split(m3s, by="cur_label"), function(x) x$gene_id)
    
    kal.abund <- kal.genes$abundance[mm$gene_id,val.samps$r_acc]
    
    valid.enrich <- suppressWarnings(gsva(kal.abund,mod.list, method="ssgsea", ssgsea.norm=F))
    
    enrich.dt <- data.table(reshape2::melt(valid.enrich, as.is=T))
    
    enrich.dt <- merge(enrich.dt[,.(cur_labels=Var1, r_acc=Var2, value)], val.samps, by="r_acc")
    
    enrich.dt[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    split.enrich <- split(enrich.dt, by=c("Biopsy", "cur_labels"))
    
    res.list <- mapply(function(val, trn){
        
        .pred.ctree(val, trn, keep.cols=c("Biopsy", "cur_labels"), rescale.value=T)
        
    }, split.enrich, ssg.res[names(split.enrich)], SIMPLIFY=F)
    
    res.list
    
}

.score.genes.mean <- function(ma.mat, geneset){
    
    common.genes <- intersect(geneset$gene_name, rownames(ma.mat))
    
    melt.mat <- data.table(reshape2::melt(ma.mat[common.genes,], as.is=T))
    names(melt.mat) <- c("gene_name","PATIENT_ID", "value")
    
    melt.mat <- merge(melt.mat, geneset, by="gene_name", allow.cartesian=T)
    
    mean.scores <- melt.mat[,.(value=mean(value)),by=.(PATIENT_ID, signature)]
    
    mean.scores
}

classify.zhang.ctree <- function(zhang, gene.info, all.samps, kal.genes, valid.samps, val.exprs.file){
    
    kal.val <- local({get(load(val.exprs.file))})
    
    #uniqify genes
    gene.info[,symb_count:=.N,by=gene_name]
    gene.info <- gene.info[symb_count == 1]
    
    stopifnot(length(setdiff(zhang$gene_name, gene.info$gene_name))==0)
    
    all.samps[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    all.samps <- all.samps[cohort == "AMTEC"]
    
    logTPM <- log2(kal.genes$kal.tpm[gene.info$gene_id,all.samps$r_acc]+1)
    rownames(logTPM) <- gene.info$gene_name
    
    sig.exprs <- .score.genes.mean(logTPM, zhang)
    
    sig.exprs <- merge(all.samps, sig.exprs, by.x="r_acc", by.y="PATIENT_ID")
    
    sig.ctrees <- lapply(split(sig.exprs, by=c("biopsy", "signature")), function(x){
        
        x.med <- median(x$value)
        
        x$preds <- factor(ifelse(x$value > x.med, "SDPR", "PD"))
        
        tmp.dt <- data.table(split=x.med, pval=NA_real_)
        
        tmp.acc <- caret::confusionMatrix(data = x$preds, reference = x$bin_resp)
        
        tmp.mcc <- mcc_vec(truth=x$bin_resp, estimate=x$preds)
        
        tmp.dt <- cbind(tmp.dt, t(tmp.acc$overall), t(tmp.acc$byClass), mcc=tmp.mcc)
        
        list(preds=x, summary=tmp.dt)
    })
    
    valid.samps[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    validTPM <- log2(kal.val$abundance[gene.info$gene_id,valid.samps$r_acc]+1)
    rownames(validTPM) <- gene.info$gene_name
    
    sig.valid.exprs <-.score.genes.mean(validTPM, zhang)
    
    valid.sigs <- merge(valid.samps, sig.valid.exprs, by.x="r_acc", by.y="PATIENT_ID")
    
    split.valid <- split(valid.sigs, by=c("Biopsy", "signature"))
    
    valid.ctrees <- mapply(function(val, trn){
        
        val$preds <- factor(ifelse(val$value > trn$summary$split, "SDPR", "PD"))
        
        tmp.acc <- caret::confusionMatrix(data = val$preds, reference = val$bin_resp)
        
        tmp.mcc <- mcc_vec(truth=val$bin_resp, estimate=val$preds)
        
        names(val) <- tolower(names(val))
        
        list(summary=data.table(val[1,.(signature, biopsy)], t(tmp.acc$overall), t(tmp.acc$byClass), mcc=tmp.mcc),
             preds=val)
        
    }, split.valid, sig.ctrees[names(split.valid)], SIMPLIFY=F)
    
    list(amtec=sig.ctrees, valid=valid.ctrees)
}

merge.clin.ssg.preds <- function(ssg, valid.ssg, clin, amtec.samps, valid.samps){
    
    ssg.prob <- rbindlist(lapply(ssg[grep("brown", names(ssg))], function(x) x$preds), idcol="desc")
    ssg.valid.prob <- rbindlist(lapply(valid.ssg[grep("brown", names(valid.ssg))], function(x) x$preds), idcol="desc")
    
    comb.ssg <- rbind(ssg.prob[,.(r_acc, ssg_preds=preds)], ssg.valid.prob[,.(r_acc, ssg_preds=preds)])
    
    comb.samps <- rbind(amtec.samps, valid.samps[,names(amtec.samps),with=F])
    
    clin <- merge(clin, comb.samps, by.x=c("De-Identified.Id", "Biopsy"), by.y=c("ptid", "Biopsy"), all.x=T)
    
    clin <- merge(clin, comb.ssg, by="r_acc", all.x=T)
    
    clin[,Consensus_Call_ssGSEA:=Consensus_Call ]
    clin[is.na(Consensus_Call) | Consensus_Call=="INDETERMINATE",Consensus_Call_ssGSEA:=ifelse(as.character(ssg_preds) == "PD", "BLIS", "BLIA")]
    
    clin[,ssg_preds:=NULL]
    
    clin[,Collapsed_Call:=Final_Call]
    
    clin[,Final_Call:=NULL]
    
    #recode collapsed
    
    clin[Collapsed_Call=="BLIA",Collapsed_Call:="Non-BLIS/LAR"]
    clin[Collapsed_Call=="BLIS",Collapsed_Call:="BLIS/LAR"]
    
    clin
}

.summarize.ctree.fit <- function(mod, valid, var.name, acc.name="Accuracy", valid.acc.name="Accuracy"){
    
    mod.summary <- cbind(rbindlist(lapply(mod, function(x) x$summary), idcol="desc"), cohort="AMTEC")
    valid.summary <- cbind(rbindlist(lapply(valid, function(x) x$summary), idcol="desc"), cohort="Validation")
    
    tmp.mod <- mod.summary[,c("cohort", "desc", acc.name),with=F]
    names(tmp.mod)[3] <- "Accuracy"
    
    tmp.valid <- valid.summary[,c("cohort", "desc", valid.acc.name),with=F]
    names(tmp.valid)[3] <- "Accuracy"
    
    comb.summary <- rbind(tmp.mod, tmp.valid)
    
    comb.summary[,biopsy:=str_match(desc, "(Bx\\d)")[,2]]
    comb.summary[,variable:=var.name]
    comb.summary
}

summarize.models <- function(bb.calls, mods, valid.mods, mzb, valid.mzb, ssg, valid.ssg, ssg.meta, valid.ssg.meta,zhang.preds){
    
    #accuracy of blia blis
    
    bb.calls[,bin_resp:=factor(ifelse(Best_Response == "PD", "PD", "SDPR"))]
    
    bb.res <- rbindlist(lapply(c(BLIA_BLIS_IHC="Consensus_Call", BLIA_BLIS_IHC_Collapsed="Collapsed_Call", BLIA_BLIS_IHC_ssGSEA="Consensus_Call_ssGSEA"), function(x){
        
        tmp <- bb.calls[,c("Biopsy", "cohort","bin_resp", x),with=F]
        tmp$call_fac <- factor(bb.calls[[x]], levels=c("LAR","BLIS", "BLIA"), labels=c("PD","PD", "SDPR"))
        
        acc.res <- tmp[,data.table(t(confusionMatrix(data=call_fac, reference=bin_resp)$overall), ssize=sum(is.na(call_fac)==F)),by=.(biopsy=Biopsy, cohort)]
        
        acc.res
    }), idcol="variable")
    
    bb.res[,desc:=paste(biopsy, variable, sep=".")]
    
    brown.summary <- .summarize.ctree.fit(mods[grep("brown", names(mods))], valid.mods[grep("brown", names(valid.mods))], "Brown_ME")
    
    mzb1.summary <- .summarize.ctree.fit(mzb, valid.mzb, "MZB1")
    
    ssg.summary <- .summarize.ctree.fit(ssg[grep("brown", names(ssg))], valid.ssg[grep("brown", names(valid.ssg))], "ssGSEA")
    
    ssg.meta.summary <- .summarize.ctree.fit(ssg.meta[grep("brown", names(ssg.meta))], valid.ssg.meta[grep("brown", names(valid.ssg.meta))], "ssGSEA_Top3")
    
    pB.summary <- .summarize.ctree.fit(zhang.preds$amtec[grepl("pB", names(zhang.preds$amtec))], zhang.preds$valid[grepl("pB", names(zhang.preds$valid))], "Zhang_pB")
    
    all.summary <- rbind(bb.res[,names(brown.summary),with=F], brown.summary, mzb1.summary, ssg.summary, ssg.meta.summary, pB.summary)
    
    all.summary
}

generate.summary.matrix <- function(bb.calls, mods, mzb, ssg, ssg.3, zhang.preds, perf.summary){
    
    #brown module
    
    brown.preds <- rbindlist(lapply(mods[grep("brown", names(mods))], function(x) x$preds))[,.(r_acc, Brown_ME=paste0("p",preds))]
    
    bb.calls.m <- merge(bb.calls, brown.preds, all.x=T, by="r_acc")
    
    #mzb1
    
    mzb.preds <- rbindlist(lapply(mzb, function(x) x$preds))[,.(r_acc, MZB1=paste0("p",preds))]
    
    bb.calls.m <- merge(bb.calls.m, mzb.preds, all.x=T, by="r_acc")
    
    #ssgsea
    
    ssg.preds <- rbindlist(lapply(ssg[grep("brown", names(ssg))], function(x) x$preds))[,.(r_acc, ssGSEA=paste0("p",preds))]
    
    bb.calls.m <- merge(bb.calls.m, ssg.preds, all.x=T, by="r_acc")
    
    #ssgsea top3
    
    ssg.3.preds <- rbindlist(lapply(ssg.3[grep("brown", names(ssg.3))], function(x) x$preds))[,.(r_acc, ssGSEA_Top3=paste0("p",preds))]
    
    bb.calls.m <- merge(bb.calls.m, ssg.3.preds, all.x=T, by="r_acc")
    
    #zhang pb
    
    pb.preds <- rbindlist(lapply(zhang.preds$amtec[grepl("pB", names(zhang.preds$amtec))], function(x) x$preds))[,.(r_acc, Zhang_pB=paste0("p",preds))]
    
    bb.calls.m <- merge(bb.calls.m, pb.preds, all.x=T, by="r_acc")
    
    bb.calls.m[,BLIA_BLIS_IHC:=Consensus_Call]
    
    bb.calls.m[,BLIA_BLIS_IHC_Collapsed:=Collapsed_Call]
    
    bb.calls.m[,BLIA_BLIS_IHC_ssGSEA:=Consensus_Call_ssGSEA ]
    
    melt.bb <- melt(id.vars=c("De-Identified.Id","Biopsy", "Short_Location", "cohort"), 
                    measure.vars=c("BLIA_BLIS_IHC",  "BLIA_BLIS_IHC_ssGSEA", "Brown_ME", "Zhang_pB", "MZB1", "ssGSEA", "ssGSEA_Top3", "Best_Response"), 
                    data=bb.calls.m, as.is=T)
    
    all.combs <- data.table(expand.grid(list(`De-Identified.Id`=unique(melt.bb$`De-Identified.Id`), Biopsy=unique(melt.bb$Biopsy), 
                                             variable=unique(melt.bb$variable)), stringsAsFactors = F))
    all.combs <- merge(all.combs, unique(melt.bb[,.(`De-Identified.Id`, cohort)]), by=c("De-Identified.Id"))
    
    melt.bb <- merge(all.combs, melt.bb, by=c("De-Identified.Id", "Biopsy", "variable", "cohort"), all=T)
    
    resp.ord <- unique(bb.calls.m[,.(`De-Identified.Id`, Best_Response)])[order(Best_Response, `De-Identified.Id`)]
    
    var.levs <- c("Best_Response", "BLIA_BLIS_IHC_ssGSEA" , "BLIA_BLIS_IHC","Zhang_pB","MZB1", "ssGSEA_Top3","ssGSEA",  "Brown_ME")
    var.labels <- c("Best Response", "Burstein+mTNBC3s", "Burstein", "Zhang pB", "MZB1", "mTNBC3s_top3","mTNBC3s",  "mTNBC3e")
    
    melt.bb[,`:=`(pat_fac=factor(`De-Identified.Id`, levels=resp.ord$`De-Identified.Id`, ordered=T), 
                  var_fac=factor(variable, levels=var.levs, labels=var.labels, ordered=T))]
    
    melt.bb[,br_cats:=""]
    melt.bb[variable == "Best_Response", br_cats:="br"]
    melt.bb[variable %in% c( "BLIA_BLIS_IHC", "BLIA_BLIS_IHC_ssGSEA"), br_cats:="bb"]
    melt.bb[,value:=ifelse(is.na(value), "NC", value)]
    
    perf.summary <- perf.summary[variable %in% var.levs]
    
    perf.summary[,var_fac:=factor(variable, levels=var.levs, labels=var.labels,ordered=T)]
    perf.summary[,Biopsy:=biopsy]
    perf.summary[,br_cats:=""]
    perf.summary[variable == "Best_Response", br_cats:="br"]
    perf.summary[variable %in% c("BLIA_BLIS_IHC", "BLIA_BLIS_IHC_ssGSEA"), br_cats:="bb"]
    perf.summary[,value:="NC"]
    
    #add in no calls
    
    nc.sum <- melt.bb[,.(NC=sum(value == "INDETERMINATE"), total=sum(value != "NC")),by=.(var_fac, Biopsy, cohort)]
    
    perf.summary <- merge(perf.summary, nc.sum, by=c("var_fac", "Biopsy", "cohort"))
    
    perf.summary[,nc_rate:=ifelse(NC > 0, paste0("(NC: ", scales::percent(NC/total, accuracy=.1), ")"), "")]
    
    list(data=melt.bb, summary=perf.summary)
    
}

amtec.exprs.de.by.biopsy <- function(amtec.abunds, amtec.rna, gene.info){
    
    amtec.mat <- amtec.abunds$kal.counts
    amtec.dge <- DGEList(amtec.mat)
    
    stopifnot((max(amtec.dge$samples$lib.size) / min(amtec.dge$samples$lib.size)) < 3)
    
    #if above is true, can use limma trend for this
    
    #from this: https://support.bioconductor.org/p/95818/
    
    amtec.rna[,`:=`(bin_resp=ifelse(short_resp == "PD", "PD", "SDPR"))]
    amtec.rna[,biop_resp:=paste(Biopsy, bin_resp, sep="_")]
    
    design <- model.matrix(~0+biop_resp, data=amtec.rna)
    rownames(design) <- amtec.rna$r_acc
    colnames(design) <- gsub("biop_resp", "", colnames(design))
    
    con <- makeContrasts(Bx1=Bx1_PD-Bx1_SDPR, Bx2=Bx2_PD-Bx2_SDPR, levels=design)
    
    keep <- filterByExpr(amtec.dge[,rownames(design)], design=design)
    
    amtec.dge <- amtec.dge[keep,,keep.lib.sizes=FALSE]
    
    amtec.dge <- calcNormFactors(amtec.dge)
    
    logCPM <- cpm(amtec.dge, log=T, prior.count=3)
    
    fit <- lmFit(logCPM[,rownames(design)], design)
    
    fit <- contrasts.fit(fit, con)
    
    fit <- eBayes(fit, trend=TRUE)
    
    res.bx1 <- topTable(fit, coef="Bx1", n=Inf)
    
    res.bx2 <- topTable(fit, coef="Bx2", n=Inf)
    
    res.bx1.dt <- merge(gene.info, data.table(gene_id=rownames(res.bx1), res.bx1), by="gene_id")
    
    res.bx2.dt <- merge(gene.info, data.table(gene_id=rownames(res.bx2), res.bx2), by="gene_id")
    
    comb.exprs.de <- rbind(
        cbind(biopsy="Bx1", res.bx1.dt),
        cbind(biopsy="Bx2", res.bx2.dt)
    )
    
    comb.exprs.de
    
}

amtec.dnet.by.biopsy <- function(l.conc, all.samps, mm){
    
    all.samps[,ptid_biop:=paste(ptid, biopsy, sep=".")]
    
    all.samps <- all.samps[cohort == "AMTEC"]
    
    mod.stats.dt <- rbindlist(lapply(l.conc$module.stats, function(x) rbindlist(lapply(x, function(y) y$per_gene), idcol="ptid_biop")), idcol="cur_labels")
    
    melt.lion <- melt(measure.vars=c("ClusterCoef", "Connectivity", "MAR"), id.vars=c("gene_id", "ptid_biop"), 
                      data=mod.stats.dt[cur_labels %in% c("Mod1 (turquoise)", "Mod2 (blue)", "Mod3 (brown)")], variable.factor=F)
    
    melt.lion[,value:=log2(value+1)]
    
    all.samps[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    all.samps[,biop_resp:=paste(biopsy, bin_resp, sep="_")]
    
    design <- model.matrix(~0+biop_resp, data=all.samps)
    rownames(design) <- all.samps$ptid_biop
    colnames(design) <- gsub("biop_resp", "", colnames(design))
    
    con <- makeContrasts(Bx1=Bx1_PD-Bx1_SDPR, Bx2=Bx2_PD-Bx2_SDPR, levels=design)
    
    net.conc.de <- rbindlist(lapply(split(melt.lion, by=c("variable")), function(x){
        
        tmp.mat <- reshape2::acast(gene_id~ptid_biop, value.var="value", data=x)
        
        fit <- lmFit(tmp.mat[,rownames(design)], design)
        
        fit <- contrasts.fit(fit, con)
        
        fit <- eBayes(fit)
        
        res.bx1 <- topTable(fit, coef="Bx1", n=Inf)
        
        res.bx2 <- topTable(fit, coef="Bx2", n=Inf)
        
        rbind(
            data.table(biopsy="Bx1",gene_id=rownames(res.bx1), res.bx1),
            data.table(biopsy="Bx2",gene_id=rownames(res.bx2), res.bx2)
        )
        
    }), idcol="network_concept")
    
    net.conc.de <- merge(mm, net.conc.de, by="gene_id")
    
    net.conc.de
    
}

amtec.exprs.de.delta <- function(amtec.abunds, amtec.rna, gene.info){
    
    amtec.mat <- amtec.abunds$kal.counts
    amtec.dge <- DGEList(amtec.mat)
    
    stopifnot((max(amtec.dge$samples$lib.size) / min(amtec.dge$samples$lib.size)) < 3)
    
    #can use limma trend for this then
    
    #from this: https://support.bioconductor.org/p/95818/
    
    amtec.rna[,`:=`(biop_fac=factor(Biopsy), pt_fac=factor(ptid), bin_resp=factor(ifelse(short_resp == "PD", "PD", "SDPR")))]
    
    design <- model.matrix(~0 + pt_fac + biop_fac:bin_resp, data=amtec.rna)
    design <- design[,-grep("Bx1", colnames(design))] # Get to full rank
    colnames(design) <- sub(":", "_", colnames(design)) # Rename for convenience
    rownames(design) <- amtec.rna$r_acc
    
    con <- makeContrasts(biop_facBx2_bin_respPD - biop_facBx2_bin_respSDPR, levels=design)
    
    keep <- filterByExpr(amtec.dge[,rownames(design)], design=design)
    
    amtec.dge <- amtec.dge[keep,,keep.lib.sizes=FALSE]
    
    amtec.dge <- calcNormFactors(amtec.dge)
    
    logCPM <- cpm(amtec.dge, log=T, prior.count=3)
    
    fit <- lmFit(logCPM[,rownames(design)], design)
    
    fit <- contrasts.fit(fit, con)
    
    fit <- eBayes(fit, trend=TRUE)
    
    res <- topTable(fit, n=Inf)
    
    res.dt <- data.table(gene_id=rownames(res), res)
    
    res.dt <- merge(gene.info, res.dt, by="gene_id")
    
    #test
    
    melt.cpm <- data.table(reshape2::melt(logCPM, as.is=T))
    names(melt.cpm) <- c("gene_id", "r_acc", "value")
    
    melt.cpm <- merge(amtec.rna, melt.cpm, by="r_acc")
    
    wide.cpm <- dcast(ptid+bin_resp~biop_fac, value.var="value", data=melt.cpm[gene_id == "ENSG00000122852.14"])
    
    wide.cpm[,diff:=Bx2-Bx1]
    
    stopifnot(isTRUE(all.equal(wide.cpm[,mean(diff),by=bin_resp][,diff(V1)],res["ENSG00000122852.14","logFC"] )))
    
    #end test
    
    cbind(biopsy="Delta", res.dt)
    
}

amtec.dnet.delta <- function(all.samps, l.conc, mm){
    
    all.samps[,ptid_biop:=paste(ptid, biopsy, sep=".")]
    
    all.samps <- all.samps[cohort == "AMTEC"]
    
    mod.stats.dt <- rbindlist(lapply(l.conc$module.stats, function(x) rbindlist(lapply(x, function(y) y$per_gene), idcol="ptid_biop")), idcol="cur_labels")
    
    melt.lion <- melt(measure.vars=c("ClusterCoef", "Connectivity", "MAR"), id.vars=c("gene_id", "ptid_biop"), 
                      data=mod.stats.dt[cur_labels %in% c("Mod1 (turquoise)", "Mod2 (blue)", "Mod3 (brown)")], variable.factor=F)
    
    melt.lion[,value:=log2(value+1)]
    
    all.samps[,`:=`(biop_fac=factor(biopsy), pt_fac=factor(ptid), bin_resp=factor(ifelse(short_resp == "PD", "PD", "SDPR")))]
    
    design <- model.matrix(~0 + pt_fac + biop_fac:bin_resp, data=all.samps)
    design <- design[,-grep("Bx1", colnames(design))] # Get to full rank
    colnames(design) <- sub(":", "_", colnames(design)) # Rename for convenience
    rownames(design) <- all.samps$ptid_biop
    
    con <- makeContrasts(biop_facBx2_bin_respPD - biop_facBx2_bin_respSDPR, levels=design)
    
    net.conc.de <- rbindlist(lapply(split(melt.lion, by=c("variable")), function(x){
        
        tmp.mat <- reshape2::acast(gene_id~ptid_biop, value.var="value", data=x)
        
        fit <- lmFit(tmp.mat[,rownames(design)], design)
        
        fit <- contrasts.fit(fit, con)
        
        fit <- eBayes(fit)
        
        tmp.res <- topTable(fit, n=Inf)
        
        data.table(gene_id=rownames(tmp.res), tmp.res)
        
    }), idcol="network_concept")
    
    net.conc.de <- merge(mm, net.conc.de, by="gene_id")
    
    #test
    
    melt.cpm <- merge(all.samps, melt.lion, by="ptid_biop")
    
    wide.cpm <- dcast(ptid+variable+bin_resp~biop_fac, value.var="value", data=melt.cpm[gene_id == "ENSG00000002933.7"])
    
    wide.cpm[,diff:=Bx2-Bx1]
    
    wide.cpm[,mean(diff),by=.(variable, bin_resp)][,diff(V1), by=.(network_concept=variable)][order(network_concept)]
    
    stopifnot(
        isTRUE(
            all.equal(
                wide.cpm[,mean(diff),by=.(variable, bin_resp)][,diff(V1), by=.(network_concept=variable)][order(network_concept)],
                net.conc.de[gene_id == "ENSG00000002933.7",.(network_concept, V1=logFC)][order(network_concept)] 
            )))
    
    #end test
    
    cbind(biopsy="Delta", net.conc.de)
    
}

.fit.multi.ctree <- function(x){
    
    tmp <- ctree(bin_resp~., data=as.data.frame(x[,c("bin_resp", "value"),with=F]), control=ctree_control( testtype ="Univariate", mincriterion = 0.95 ,minbucket=3, minsplit=3, maxdepth=1))
    
    tmp.split <- split_node(node_party(tmp[[1]]))
    
    if (is.null(tmp.split$breaks)){
        tmp.dt <- data.table(split=NA)
    }else{
        tmp.dt <- data.table(split=tmp.split$breaks)
    }
    
    tmp.acc <- caret::confusionMatrix(data = predict(tmp), reference = x$bin_resp)
    
    tmp.dt <- cbind(tmp.dt, t(tmp.acc$overall))
    
    tmp.2 <- ctree(bin_resp~., data=as.data.frame(x[,c("bin_resp", "value", "z"),with=F]), control=ctree_control( testtype ="Univariate", mincriterion = 0.95 ,minbucket=1, minsplit=1, maxdepth=2))
    
    tmp.2.acc <- caret::confusionMatrix(data = predict(tmp.2), reference = x$bin_resp)
    
    tmp.3 <- ctree(bin_resp~., data=as.data.frame(x[,c("bin_resp", "z"),with=F]), control=ctree_control( testtype ="Univariate", mincriterion = 0.95 ,minbucket=1, minsplit=1, maxdepth=2))
    
    tmp.3.acc <- caret::confusionMatrix(data = predict(tmp.3), reference = x$bin_resp)
    
    tmp.dt <- cbind(tmp.dt, gene_only_acc=tmp.3.acc$overall["Accuracy"], both_acc=tmp.2.acc$overall["Accuracy"])
    
    list(summary=tmp.dt, fit=tmp.2, preds=cbind(x, preds=predict(tmp.2)))
}

ctree.gene.net.exprs <- function(all.samps, l.conc, exprs){
    
    all.samps[,ptid_biop:=paste(ptid, biopsy, sep=".")]
    
    all.samps <- all.samps[cohort == "AMTEC"]
    
    mod.stats.dt <- rbindlist(lapply(l.conc$module.stats, function(x) rbindlist(lapply(x, function(y) y$per_gene), idcol="ptid_biop")), idcol="cur_labels")
    
    lion.k <- merge(all.samps[,.(ptid_biop, r_acc, biopsy, short_resp)],  mod.stats.dt, by="ptid_biop")
    
    lion.k[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    melt.lion <- melt(measure.vars=c("ClusterCoef", "Connectivity", "MAR"), id.vars=c("gene_id","cur_labels", "ptid_biop", "r_acc", "biopsy", "short_resp" ,"bin_resp"), data=lion.k[cur_labels != "Mod0 (grey)"], variable.factor=F)
    melt.lion[,value:=log2(value+1)]
    
    melt.lion <- merge(melt.lion, exprs[,.(gene_id, r_acc, z)], by=c("gene_id", "r_acc"))
    
    gene.res <- lapply(split(melt.lion, by=c("biopsy", "cur_labels", "gene_id", "variable")), function(x){
        
        .fit.multi.ctree(x)
        
    })
    
    gene.res
    
}

classify.gene.net.delta.ctree <- function(mm, all.rna, lion.c, scaled.exprs, valid.rna, valid.lion.c){
    
    gene.spec <- rbindlist(lapply(unlist(lion.c$module.stats, recursive=F), function(x) x$per_gene), idcol="desc")
    gene.spec[,c("cur_label", "ptid", "biopsy"):=tstrsplit(desc, "\\.")]
    melt.gs<- melt(id.vars=c("cur_label", "gene_id", "ptid", "biopsy"), measure.vars=c("Connectivity", "ClusterCoef", "MAR"), data=gene.spec)
    melt.gs[,value:=log2(value+1)]
    gene.spec.m <- merge(all.rna, melt.gs, by=c("ptid", "biopsy"))
    wide.gene.spec <- dcast(ptid+short_resp+gene_id+cur_label+variable~biopsy, value.var="value", data=gene.spec.m)
    wide.gene.spec[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    wide.gene.spec[,diff:=Bx2-Bx1]
    
    #note expression here relates to scaled as scaled differences ends up being (bx2-bx1)/s with same test stat and p-value
    wide.se <- dcast(cohort+ptid+short_resp+gene_id~biopsy, value.var="z", data=scaled.exprs[cohort == "AMTEC" & gene_id %in% mm$gene_id ])
    wide.se[,diff:=Bx2-Bx1]
    
    #fit ctree
    
    comb.gene <- merge(wide.gene.spec[,.(ptid, gene_id, cur_label, variable=as.character(variable), bin_resp, value=diff)], wide.se[,.(ptid, gene_id, z=diff)], by=c("ptid","gene_id"), allow.cartesian=T)
    
    cg.res <- lapply(split(comb.gene, by=c("gene_id", "variable", "cur_label")), function(x){
        .fit.multi.ctree(x)
    })
    
    valid.gene.spec <- rbindlist(lapply(unlist(valid.lion.c$module.stats, recursive=F), function(x) x$per_gene), idcol="desc")
    valid.gene.spec[,c("cur_label", "ptid", "biopsy"):=tstrsplit(desc, "\\.")]
    valid.melt.gs<- melt(id.vars=c("cur_label", "gene_id", "ptid", "biopsy"), measure.vars=c("Connectivity", "ClusterCoef", "MAR"), data=valid.gene.spec)
    valid.melt.gs[,value:=log2(value+1)]
    valid.gene.spec.m <- merge(valid.rna, valid.melt.gs, by=c("ptid", "biopsy"))
    valid.wide.gene.spec <- dcast(ptid+short_resp+gene_id+cur_label+variable~biopsy, value.var="value", data=valid.gene.spec.m)
    valid.wide.gene.spec[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    valid.wide.gene.spec[,diff:=Bx2-Bx1]
    
    valid.wide.se <- dcast(cohort+ptid+short_resp+gene_id~biopsy, value.var="z", data=scaled.exprs[cohort == "Validation" & gene_id %in% mm$gene_id ])
    valid.wide.se[,diff:=Bx2-Bx1]
    
    valid.comb.gene <- merge(valid.wide.gene.spec[,.(ptid, gene_id, cur_label, variable=as.character(variable), bin_resp, value=diff)], valid.wide.se[,.(ptid, gene_id, z=diff)], by=c("ptid","gene_id"), allow.cartesian=T)
    
    #limit to the paired validation samples
    valid.comb.gene <- valid.comb.gene[is.na(value)==F]
    
    split.valid <- split(valid.comb.gene, by=c("gene_id", "variable", "cur_label"))
    
    valid.list <- mapply(function(val, trn){
        
        .pred.ctree(val, trn, keep.cols=c("gene_id", "variable", "cur_label"))
        
    }, split.valid, cg.res[names(split.valid)], SIMPLIFY=F)
    
    list(amtec=cg.res, valid=valid.list)
    
}

write.de <- function(de.table, filename){
    
    de.table[,`:=`(AveExpr=NULL, B=NULL)]
    
    openxlsx::write.xlsx(list(de.table), file=file.path("output", filename))
    
    file.path("output", filename)
}
