parse.trans.info <- function(trans.file){
    trans.info <- fread(trans.file, skip=1, header=F)
    trans.info <- trans.info[V3 == "transcript"]
    trans.info[,gene_id:=sub("\"", "", sub("\\s*gene_id\\s+\"", "", sapply(strsplit(V9, ";"), function(x) x[grepl("gene_id", x)])))]
    trans.info[,transcript_id:=sub("\"", "", sub("\\s*transcript_id\\s+\"", "", sapply(strsplit(V9, ";"), function(x) x[grepl("transcript_id", x)])))]
    trans.info[,gene_name:=sub("\"", "", sub("\\s*gene_name\\s+\"", "", sapply(strsplit(V9, ";"), function(x) x[grepl("gene_name", x)])))]
    trans.info[,gene_type:=sub("\"", "", sub("\\s*gene_type\\s+\"", "", sapply(strsplit(V9, ";"), function(x) x[grepl("gene_type", x)])))]
    
    trans.info[,V9:=NULL]
    
    trans.info
}

gene.info.from.trans <- function(trans.info){
    
    gene.info <- unique(trans.info[,.(gene_id, gene_name, gene_type)])
    
    gene.info
    
}

read.matrix.from.file <- function(inp.file){
    
    inp.dt <- fread(inp.file)
    
    inp.mat <- as.matrix(inp.dt[,-1,with=F])
    rownames(inp.mat) <- inp.dt$V1
    
    inp.mat
}


process.subtypes <- function(id.map.file, pancan.clin.file){
    
    id.map <- fread(id.map.file)
    brca.ids <- id.map[Disease=="BRCA"]
    brca.ids[,c("project","tss", "participant", "sample_vial", "portion_analyte", "plate", "center"):=tstrsplit(AliquotBarcode, "\\-")]
    
    brca.ids[,`:=`(sample=substr(sample_vial, 1, 2), vial=substr(sample_vial, 3,3))]
    brca.ids[,sample_vial:=NULL]
    
    tumor.samps <- brca.ids[sample == "01"]
    tumor.samps[,Sample.ID:=paste(project, tss, participant, sample, sep="-")]
    
    #tumor.samps[,.N]#1135
    
    pancan <- fread(pancan.clin.file)
    names(pancan) <- make.names(names(pancan))
    
    tumor.samps <- merge(tumor.samps, pancan, by="Sample.ID", all.x=T, all.y=F)
    
    #remove duplicate samples
    
    tumor.samps[,num_samps:=.N,by=Patient.ID]
    
    tumor.annots <- tumor.samps[num_samps == 1 & Subtype != ""]
    
    tumor.annots
}

filter.subtypes <- function(tcga.annots){
    
    filt.annots <- tcga.annots[Subtype == "BRCA_Basal"]
    
    filt.annots
}

tcga.kallisto.abund <- function(tpm.mat, count.mat, tumor.annots){
    
    tpm.mat <- tpm.mat[,tumor.annots$CGHubAnalysisID]
    
    count.mat <- count.mat[,tumor.annots$CGHubAnalysisID]
    
    stopifnot(all(colnames(count.mat) == colnames(tpm.mat)))
    
    #as we only want to sum the tpms to gene-level using `tximport`, create a placeholder for length
    ph.length <- matrix(1, ncol=ncol(count.mat), nrow=nrow(count.mat), dimnames=list(rownames(count.mat), colnames(count.mat)))
    
    t2g <- setNames(as.data.frame(tstrsplit(rownames(tpm.mat), "\\|", keep=c(1:2, 6))), c("transcript_id", "gene_id", "symbol"))
    
    #in case we want the counts too
    
    osf.abund <- summarizeToGene(
        object=list(abundance=tpm.mat, counts=count.mat, length=ph.length),
        tx2gene=t2g,
        countsFromAbundance = c("scaledTPM"),
        ignoreAfterBar=T
    )
    
    dge <- DGEList(osf.abund$counts, remove.zeros=F)
    
    stopifnot(all(colnames(dge) == tumor.annots$CGHubAnalysisID))
    
    list(counts=dge, tpm=osf.abund$abundance)
    
}

filter.tcga.kallisto.tpm <- function(tcga, filt.annots){
    
    filt.dge <- tcga$tpm[,filt.annots$CGHubAnalysisID]
    
    #initial filtering, originally for counts, but seems to work for TPM
    wgcna.keep <- (rowSums(filt.dge < 10) >= (.9*ncol(filt.dge))) == F
    
    rawTPM <- log2(filt.dge[wgcna.keep,] + 1)
    
    #determing low expression outliers based on medians
    
    samp.meds <- matrixStats::colMedians(rawTPM)
    
    med.zs <- scale(samp.meds, scale=T)
    
    rm.samps <- med.zs < -1.5
    
    stopifnot(all(colnames(rawTPM) == colnames(filt.dge)))
    
    filt.dge <- filt.dge[,rm.samps == F]
    
    wgcna.keep <- (rowSums(filt.dge < 10) >= (.9*ncol(filt.dge))) == F
    
    list(dge=filt.dge, wgcna.keep=wgcna.keep, num_removed=sum(rm.samps))
    
}

form.dge.kallisto <- function(kal.counts){
    
    log2(kal.counts$kal.tpm + 1)
}

make.combined.annot <- function(dtype.rna, tcga, filt.tcga, comb.cohort){
    
    cur.tcga <- filt.tcga[CGHubAnalysisID %in% colnames(tcga$dge)]
    
    comb.annot <- merge(dtype.rna[,.(ptid, cohort=comb.cohort, biopsy=Biopsy, short_resp, r_acc)], cur.tcga[,.(ptid=Patient.ID, cohort="TCGA", r_acc=CGHubAnalysisID )], by=c("ptid", "cohort", "r_acc"), all=T)
    
    stopifnot(comb.annot[,.N] == (dtype.rna[,.N] + cur.tcga[,.N]))
    
    comb.annot
    
}

#batch correct using TCGA as the reference
combine.abundance <- function(tcga, amtec, comb.annot){
    
    tcga.ltpm <- log2(tcga$dge[tcga$wgcna.keep,]+1)
    
    init.norm <- cbind(tcga.ltpm, amtec[rownames(tcga.ltpm),])[,comb.annot$r_acc]
    
    use.batch <- ifelse(comb.annot$cohort == "TCGA", 2, 1)
    
    combat.norm <- ComBat(dat=init.norm, batch=use.batch, mod=NULL, par.prior=TRUE, prior.plots=F, mean.only=F, ref.batch=2)
    
    combat.norm
}

get.validation.exprs <- function(kal.count.file){
    
    load(kal.count.file)
    
    log2(kal.genes$abundance + 1)
}

#orginal paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6349443/#SD1
call.blia.blis <- function(centroids, exprs, gene.info){
    
    centr.gene.info <- merge(centroids, gene.info, by.x="Gene", by.y="gene_name")
    
    stopifnot(centr.gene.info[,.N,by=Gene][,all(N==1)])
    
    centr.gene.info <- centr.gene.info[gene_id %in% rownames(exprs)]
    
    centr.mat <- as.matrix(centr.gene.info[,.(LAR, MES, BLIA, BLIS)])
    rownames(centr.mat) <- centr.gene.info$gene_id
    
    full.cor <- cor(exprs[rownames(centr.mat),], centr.mat, method="spearman")
    full.cor.dt <- data.table(reshape2::melt(full.cor, as.is=T))
    
    full.top <- full.cor.dt[,.SD[order(-value)][1],by=Var1]
    full.second <- full.cor.dt[,.SD[order(-value)][2],by=Var1]
    
    full.toptwo <- merge(full.top[,.(r_acc=Var1, call=Var2, call_cor=value)], full.second[,.(r_acc=Var1, second_call=Var2, second_call_cor=value)], by="r_acc")
    
    # < .1 for mixed subtype is based off of Ding et al paper linked to above.
    full.toptwo[,comp_call:=ifelse((call_cor - second_call_cor) >= .1, call, paste(call, second_call, sep="/"))]
    
    full.toptwo[,final_call:=sapply(strsplit(comp_call, "\\/"), function(x) paste(sort(x), collapse="/"))]
    
    full.toptwo[,comp_call:=NULL]
    
    list(call=full.toptwo, matrix=full.cor)
    
}

batch.call.blia.blis <- function(centroids, comb.ltpm,dtype.rna, gene.info, use.cohort){
    
    use.cohort <- dtype.rna[cohort == use.cohort]
    
    logTPM <- comb.ltpm[,use.cohort$r_acc]
    
    call.blia.blis(centroids, logTPM, gene.info)
    
}

get.module.membership <- function(gene.univ, mod.file, exprs.file, blia.blis.cents){
    
    load(mod.file)
    load(exprs.file)#datExpr.2, datExpr.5
    
    exprs.size <- sub("k", "", strsplit(basename(mod.file), "_")[[1]][2])
    
    print(exprs.size)
    
    if (exprs.size == "two"){
        use.exprs <- datExpr.2
    }else if (exprs.size == "five"){
        use.exprs <- datExpr.5
    }else{
        stop("Unexpected WGCNA size")
    }
    
    #kme 
    
    kmes <- signedKME(
        datExpr=use.exprs, 
        datME=MEs, 
        outputColumnName = "kME", 
        corFnc = "bicor", 
        corOptions = "maxPOutliers = .1")
    
    kme.dt <- data.table(reshape2::melt(as.matrix(kmes), as.is=T))
    
    mod.map <- data.table(gene_id=colnames(datExpr.2), label=paste0("kME", moduleLabels), color=moduleColors)
    
    kme.dt <- merge(mod.map, kme.dt[,.(gene_id=Var1, label=Var2, kme=value)],by=c("gene_id","label"))
    
    kme.dt[,cur_label:=paste0(sub("kME", "Mod", label), " (", color, ")")]
    
    col.map <- unique(kme.dt[,.(cur_label, color)])
    
    kme.annot <- merge(gene.univ[,.(gene_id, gene_name, gene_type)], kme.dt, by="gene_id")
    
    kme.annot[,label:=sub("kME", "", label)]
    
    kme.annot[,in_blia_blis:=gene_name %in% blia.blis.cents$Gene]
    
    kme.annot
    
}

proc.mes <- function(comb.ltpm, dtype.rna, wgcna.map, use.cohort){
    
    
    mod.map <- unique(wgcna.map[,.(color, label, cur_label)])
    
    tcga.exprs <- t(comb.ltpm[,dtype.rna[cohort=="TCGA"]$r_acc])
    
    tcga.prs <- lapply(split(wgcna.map, by="cur_label"), function(x){
        
        z.exprs <- scale(tcga.exprs[,x$gene_id])
        
        av.expr <- rowMeans(z.exprs)
        
        tmp.prs <- prcomp(z.exprs, center=F, scale.=F)
        
        tmp.me <- data.table(ptid=rownames(tmp.prs$x), tmp.prs$x[,paste0("PC", 1:2)])
        
        if(sign(cor(tmp.prs$x[names(av.expr),"PC1"], av.expr)) < 0){
            
            tmp.me$PC1 <- -tmp.me$PC1
            
        }
        
        list(mes=tmp.me, pca=tmp.prs, centers=attr(z.exprs,"scaled:center"), scales=attr(z.exprs,"scaled:scale"))
        
    })
    
    tcga.mes <- rbindlist(lapply(tcga.prs, function(x){
        
        x$mes
        
    }), idcol="cur_label")
    
    
    #predict the amtec MEs based on eigenvectors from TCGA
    
    amtec.exprs <- t(comb.ltpm[,dtype.rna[cohort == use.cohort]$r_acc])
    
    amtec.mes <- rbindlist(lapply(tcga.prs, function(x){
        
        z.exprs <- scale(amtec.exprs[,names(x$centers)], center=x$centers, scale=x$scales)
        
        av.expr <- rowMeans(z.exprs)
        
        tmp.preds <- z.exprs[,rownames(x$pca$rotation)] %*% x$pca$rotation
        
        tmp.me <- data.table(ptid=rownames(tmp.preds), tmp.preds[,paste0("PC", 1:2)])
        
        if(sign(cor(tmp.preds[names(av.expr),"PC1"], av.expr)) < 0){
            
            tmp.me$PC1 <- -tmp.me$PC1
            
        }
        
        tmp.me
        
    }), idcol="cur_label")
    
    comb.mes <- rbind(cbind(tcga.mes, cohort="TCGA"),
                      cbind(amtec.mes, cohort=use.cohort))
    
    
    comb.mes[cur_label != "Mod0 (grey)"]
    
}

write.mm <- function(mm){
    
    openxlsx::write.xlsx(list(mm[,.(gene_id, gene_name, gene_type, module_name=cur_label, kme)]), file="output/module_membership.xlsx")
    
    "output/module_membership.xlsx"
}

#parse Table S1 of Tamborero et al 2018: https://aacrjournals.org/clincancerres/article/24/15/3717/80876/A-Pan-cancer-Landscape-of-Interactions-between
compile.tamborero.sets <- function(file){
    
    sets <- data.table(openxlsx::read.xlsx(file, sheet="TableS1A"))
    
    cast.sets <- melt(data=sets, id.vars=c("Cell.type", "Cell.type.abbreviation"), measure.vars=paste0("X", 4:42), variable.factor=F, value.name="gene_name")
    
    cast.sets <- cast.sets[is.na(gene_name) == F]
    
    cast.sets[,variable:=NULL]
    
    cast.sets
    
}

get.scaled.exprs <- function(gene.info, comb.lcpm, main.samps, valid.lcpm, val.samps){
    
    #main
    
    main.exprs <- data.table(reshape2::melt(comb.lcpm, as.is=T))
    colnames(main.exprs) <- c("gene_id", "r_acc", "exprs")
    
    main.exprs <- merge(main.samps, main.exprs, by=c("r_acc"))
    
    tcga.exprs <- main.exprs[cohort == "TCGA", .(tcga_mean=mean(exprs), tcga_sd=sd(exprs)),by=gene_id]
    
    amtec.exprs <- merge(main.exprs[cohort == "AMTEC"], tcga.exprs, by="gene_id")
    
    #validation
    
    
    valid.exprs <- data.table(reshape2::melt(valid.lcpm, as.is=T))
    colnames(valid.exprs) <- c("gene_id", "r_acc", "exprs")
    
    valid.exprs <- merge(val.samps, valid.exprs, by="r_acc")
    
    valid.tcga.exprs <- valid.exprs[cohort == "TCGA", .(tcga_mean=mean(exprs), tcga_sd=sd(exprs)),by=gene_id]
    
    valid.main.exprs <- merge(valid.exprs[cohort == "Validation"], valid.tcga.exprs, by="gene_id")
    
    comb.exprs <- rbind(amtec.exprs[,.(gene_id, r_acc, biopsy, ptid, short_resp, cohort, exprs, tcga_mean, tcga_sd)], valid.main.exprs[,.(gene_id, r_acc, biopsy, ptid, short_resp, cohort, exprs, tcga_mean, tcga_sd)])
    
    comb.exprs[,z:=(exprs-tcga_mean)/tcga_sd]
    
    comb.exprs <- merge(gene.info, comb.exprs, by="gene_id")
    
    comb.exprs
}


