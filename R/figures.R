.add.p.dots <- function(dt){
    dt[,dots:="N.S."]
    dt[p.value < .05, dots:="*"]
    dt[p.value < .01, dots:="**"]
    dt[p.value < .001, dots:="***"]
    
    dt
}

mm.hallmark.plot <- function(hall.enrich){
    
    #limit to only modules with significance
    
    sig.mods <- hall.enrich[padj < .05,.N,by=cur_label]
    
    hall.enrich <- hall.enrich[cur_label %in% sig.mods$cur_label]
    
    hall.enrich[,num_genes:=as.numeric(sapply(strsplit(GeneRatio, "\\/"), function(x) x[1]))]
    
    hall.sum <- hall.enrich[,.(minp=min(pvalue)),by=ID][order(-minp)]
    
    #mod.sum <- hall.enrich[,.(minp=min(pvalue)),by=cur_label][order(minp)]
    
    hall.enrich[,id_fac:=factor(ID, levels=hall.sum$ID, ordered=T, labels=sub("HALLMARK_", "", hall.sum$ID))]
    hall.enrich[,mod_fac:=factor(cur_label, ordered=T)]#levels=mod.sum$cur_label, 
    
    mod.cols <- hall.enrich[,.N,by=cur_label]
    mod.cols[,cols:=str_match(cur_label, "\\((\\w+)\\)")[,2]]
    
    hall.enrich[,mod_fill:=ifelse(padj < .05, cur_label, "none")]
    
    p1 <- ggplot(data=hall.enrich[pvalue < .05], mapping=aes(x=mod_fac, y=id_fac, size=num_genes, fill=mod_fill)) + 
        geom_point(shape=21) + 
        scale_fill_manual(values=c(setNames(mod.cols$cols, mod.cols$cur_label), "none"="white"), guide="none") +
        scale_size_continuous(name = "Number of Genes") +
        guides(size=guide_legend(position = "bottom")) +
        theme_bw(base_size = 10) + ylab("") + xlab("") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    ggsave(p1, file="figures/hallmark_fig.pdf", width=4.2, height=3.5)
    
    "figures/hallmark_fig.pdf"
}

tcga.bb.plot <- function(comb.mes, batch.bl, dtype.rna){
    
    comb.mes <- merge(comb.mes[,.(cur_label, r_acc=ptid, PC1, PC2)], dtype.rna, by="r_acc")
    
    batch.bl <- batch.bl[final_call %in% c("BLIA", "BLIS")==T]
    
    tcga.mes.bl <- merge(batch.bl, comb.mes, by="r_acc")
    
    bl.diffs <- tcga.mes.bl[,broom::tidy(t.test(PC1~final_call)), by=cur_label]
    
    bl.diffs[,sig:=p.value < .05]
    
    bl.diffs <- .add.p.dots(bl.diffs)
    
    tcga.mes.bl[,mod_fac:=factor(cur_label, levels=bl.diffs[order(p.value),cur_label], ordered=T)]
    bl.diffs[,mod_fac:=factor(cur_label, levels=bl.diffs[order(p.value),cur_label], ordered=T)]
    
    bl.plot <- ggplot(data=tcga.mes.bl[cur_label %in% bl.diffs[sig==T]$cur_label], mapping=aes(x=final_call, y=PC1)) + 
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(width=.25) + 
        facet_wrap(~mod_fac, ncol=3) + 
        annotate("segment", x=1, xend=2, y=31, yend=31) +
        geom_text(data=bl.diffs[sig == T], mapping=aes(x=1.5, y=32, label=dots), size=10) +
        theme_bw(base_size = 10) + xlab("") + ylab("Module Eigengene") + ylim(c(-30, 35))
    
    ggsave(bl.plot, file="figures/blia_s_by_me_v2.pdf", width=7, height=3)
    
    openxlsx::write.xlsx(list(stats=bl.diffs, ns=tcga.mes.bl[,.N,by=cur_label]), file="output/blia_s_by_me_stats.xlsx")
    
    "figures/blia_s_by_me_v2.pdf"
    
}

brown.kme.figure <- function(mm, all.rna, cell.marker.file){
    
    cellmarker <- data.table(openxlsx::read.xlsx(cell.marker.file))
    
    bc.cellmarker <- cellmarker[is.na(Symbol)==F]
    
    bc.cell.sum <- bc.cellmarker[,.(.N,num_in_bc=sum(grepl("[Bb]reast", cancer_type))),by=.(Symbol, cell_name)]
    
    bc.cell.sum <- bc.cell.sum[,.SD[N == max(N)],by=.(Symbol)]
    
    top.cell.assoc <- bc.cell.sum[,.(cell_type=paste(sort(cell_name), collapse=";"), N=max(N), num_in_bc=max(num_in_bc)),by=.(Symbol)]
    
    cell.mm <- merge(top.cell.assoc[N >= 5], mm[cur_label %in% c("Mod3 (brown)")], by.x="Symbol", by.y="gene_name", all.x=F, all.y=T)
    
    cell.mm[,ct_count:=.N,by=cell_type]
    
    cell.mm[ct_count < 5, cell_type:=NA_character_]
    
    cell.mm[,disp_genes:=ifelse(is.na(cell_type)==F, Symbol, "")]
    
    pos <- position_jitter(seed = 122, height=.25, width=0)
    
    p1 <- ggplot(data=cell.mm, mapping=aes(x=kme, y=cur_label)) + geom_boxplot(outlier.shape=NA) + 
        geom_point(position = pos, size=3, shape=21,mapping=aes(fill=cell_type, alpha=is.na(cell_type))) + 
        geom_label_repel(mapping=aes(label=disp_genes, fill=cell_type), 
                         position = pos, min.segment.length = 0, box.padding = 0.5, max.overlaps=Inf, size=2.88, force=10) +
        scale_fill_discrete(name="Cell Type") +
        scale_alpha_manual(values=c(`TRUE`=.5, `FALSE`=1), guide="none") +
        guides(fill=guide_legend(position = "bottom")) +
        xlab("kME") + ylab("") + coord_cartesian(xlim=c(.3, 1.1)) + theme_bw(base_size=10)
    
    ggsave(p1, file="figures/mod_kme_plot.pdf", width=4, height=4)
    
    "figures/mod_kme_plot.pdf"
    
}

brown.bcell.cors <- function(gsets, amtec.exprs, mes, gene.info, all.rna, resp.cols){
    
    amtec.exprs <- amtec.exprs[matrixStats::rowMads(amtec.exprs) > 0,]
    
    mes[,r_acc:=ptid]
    mes[,ptid:=NULL]
    
    gene.info[,count:=.N,by=gene_name]
    
    stopifnot(gene.info[gene_name %in% gsets$gene_name,all(count == 1)])
    
    gene.info <- gene.info[count == 1]
    
    gene.info <- gene.info[gene_id %in% rownames(amtec.exprs)]
    
    amtec.exprs <- amtec.exprs[gene.info$gene_id,]
    rownames(amtec.exprs) <- gene.info$gene_name
    
    gset.list <- split(gsets$gene_name, gsets$Cell.type)
    
    immune.sets.res <- gsva(expr=amtec.exprs, gset.idx.list=gset.list)
    
    melt.res <- data.table(reshape2::melt(immune.sets.res, as.is=T))
    
    melt.res <- merge(all.rna, melt.res, by.x="r_acc", by.y="Var2")
    
    melt.res <- merge(melt.res, mes, by=c("r_acc","cohort"), allow.cartesian=T)
    
    stopifnot(melt.res[,.N,by=.(cur_label, Var1, biopsy)][,all(N==13)])
    
    melt.cors <- melt.res[,broom::tidy(cor.test(x=PC1, y=value)),by=.(cur_label, Var1, biopsy)]
    
    melt.cors <- .add.p.dots(melt.cors)
    
    melt.cors[,fdr:=p.adjust(p.value)]
    
    p1 <- ggplot(data=melt.res[cur_label == "Mod3 (brown)" & Var1 == "B cells"], mapping=aes(x=PC1, y=value)) +
        geom_point() +
        geom_smooth(method="lm", formula=y~x, se=F) +
        geom_text(data=melt.cors[cur_label == "Mod3 (brown)" & Var1 == "B cells"], mapping=aes(x=5, y=.55, label=dots), size=10)+
        facet_grid(biopsy~.) +
        theme_bw(base_size = 10) + xlab("Mod3 (brown)") + ylab("B cells (GSVA Score)")
    
    ggsave(p1, file="figures/tamb_bcells_vs_brown.pdf", width=2, height=3)
    
    openxlsx::write.xlsx(list(cors=melt.cors, ns=melt.res[,.N,by=.(cur_label, Var1, biopsy)]), file="output/tamb_bcells_vs_brown_stats.xlsx")
    
    "figures/tamb_bcells_vs_brown.pdf"
}

brown.mihc.cors <- function(comb.mes, dtype.rna, filt.mihc.file){
    
    comb.mes <- merge(comb.mes[,.(cur_label, r_acc=ptid, PC1, PC2)], dtype.rna, by="r_acc")
    
    amtec.mes <- comb.mes[cohort == "AMTEC"]
    
    amtec.mihc <- data.table(openxlsx::read.xlsx(filt.mihc.file))
    
    amtec.mihc <- merge(amtec.mes, amtec.mihc[,.(ptid=PatientID, biopsy, cell_type, mean_val=`cells/mm2`)],by=c("ptid", "biopsy"), allow.cartesian=T)
    
    #brown and cd20
    
    bx.cors <- amtec.mihc[,as.data.table(broom::tidy(cor.test(x=PC1, y=mean_val))),by=.(cur_label, cell_type ,biopsy)]
    
    bx.cors <- .add.p.dots(bx.cors)
    
    bcd20.p <- ggplot(data=amtec.mihc[cur_label == "Mod3 (brown)" & cell_type == "CD20"], mapping=aes(x=PC1, y=log10(mean_val))) + 
        geom_point() +
        facet_grid(biopsy~.) +
        geom_smooth(method="lm", se=F) +
        geom_text(data=bx.cors[cur_label == "Mod3 (brown)" & cell_type == "CD20"], mapping=aes(x=5, y=2.5, label=dots, size=dots)) +
        scale_size_manual(values=c(`*`=10, `**`=10, `***`=10, `N.S.`=5), guide="none") +
        theme_bw(base_size = 10) + ylab(bquote(mIHC~CD20~log10(cells/mm^{2}))) + xlab("Mod3 (Brown)")
    
    ggsave(bcd20.p, file="figures/mihc_bcells_vs_brown.pdf", width=2, height=3)
    
    openxlsx::write.xlsx(list(cors=bx.cors, ns=amtec.mihc[,.N,by=.(cur_label, cell_type ,biopsy)]), file="output/mihc_bcells_vs_brown_stats.xlsx")
    
    "figures/mihc_bcells_vs_brown.pdf"
    
}

plot.modules.by.response <- function(all.rna, comb.mes, resp.cols, me.classif, use.biopsy){
    
    comb.mes <- merge(all.rna, comb.mes, by.x=c("r_acc", "cohort"), by.y=c("ptid",  "cohort"))
    
    comb.mes[,short_resp:=ifelse(is.na(short_resp), "N/A", short_resp)]
    comb.mes[,biopsy:=ifelse(is.na(biopsy), "N/A", biopsy)]
    
    resp.cols <- c(resp.cols, "N/A"="grey")
    
    pred.loc <- me.classif[[paste(use.biopsy,"Mod3 (brown)", sep=".")]]$summary$split
    
    mod.pca.amtec.1 <- ggplot(data=comb.mes[cur_label=="Mod3 (brown)" & ((cohort == "TCGA") | (cohort == "AMTEC" & biopsy==use.biopsy))], mapping=aes(x=PC1, y=PC2, fill=short_resp,  alpha=cohort)) + 
        geom_point(shape=21, size=3) +
        scale_fill_manual(values=resp.cols, name="Response") +
        geom_vline(xintercept=pred.loc, linetype="dashed")+
        scale_alpha_manual(values=c(TCGA=.25, AMTEC=1), name="Cohort") +
        theme_bw() + ggtitle(paste("Mod3 (brown) AMTEC", use.biopsy))
    
    ggsave(mod.pca.amtec.1, file=paste0("figures/tcga_amtec_brown_pca_delta_",use.biopsy,".pdf"), width=4, height=2.5)
    
    paste0("figures/tcga_amtec_brown_pca_delta_",use.biopsy,".pdf")
    
}

plot.all.modules.by.response <- function(all.rna, valid.rna, comb.mes, val.mes, resp.cols, me.classif, plot.line=T, file.suf=""){
    
    comb.mes <- merge(all.rna, comb.mes, by.x=c("r_acc", "cohort"), by.y=c("ptid",  "cohort"))
    comb.mes[,biopsy:=ifelse(is.na(biopsy), "N/A", biopsy)]
    
    #also validation data
    
    val.mes <- merge(valid.rna, val.mes, by.x=c("r_acc", "cohort"), by.y=c("ptid",  "cohort"))
    val.mes[,biopsy:=ifelse(is.na(biopsy), "N/A", biopsy)]
    
    pred.loc <- rbindlist(lapply(me.classif[grepl("Mod[123]",names(me.classif))], function(x) x$summary), idcol="desc")
    pred.loc[,c("biopsy", "cur_label"):=tstrsplit(desc, "\\.")]
    
    bx.list <- lapply(c("Bx1", "Bx2"), function(use.biopsy){
        
        mod.pca.amtec.1 <- ggplot(data=comb.mes[(cur_label %in% pred.loc$cur_label) & ((cohort == "TCGA") | (cohort == "AMTEC" & biopsy==use.biopsy))], mapping=aes(x=PC1, y=PC2, fill=short_resp,  alpha=cohort)) + 
            geom_point(shape=21, size=3) +
            scale_fill_manual(values=resp.cols, name="AMTEC Response", na.value="grey") +
            scale_alpha_manual(values=c(TCGA=.25, AMTEC=1), name="Cohort") +
            facet_grid(~cur_label) +
            theme_bw() + ggtitle(paste("AMTEC", use.biopsy))
        
        if (plot.line){
            mod.pca.amtec.1 <- mod.pca.amtec.1 + geom_vline(data=pred.loc[biopsy == use.biopsy], mapping=aes(xintercept=split), linetype="dashed")
        }
        
        mod.pca.valid.1 <- ggplot(data=val.mes[(cur_label %in% pred.loc$cur_label) & ((cohort == "TCGA") | (cohort == "Validation" & biopsy==use.biopsy))], mapping=aes(x=PC1, y=PC2, fill=short_resp,  alpha=cohort)) + 
            geom_point(shape=21, size=3) +
            scale_fill_manual(values=resp.cols, name="AMTEC Response", na.value="grey") +
            scale_alpha_manual(values=c(TCGA=.25, Validation=1), name="Cohort") +
            facet_wrap(~cur_label) +
            theme_bw() + ggtitle(paste("Validation", use.biopsy))
        
        if (plot.line){
            mod.pca.valid.1 <- mod.pca.valid.1 + geom_vline(data=pred.loc[biopsy == use.biopsy], mapping=aes(xintercept=split), linetype="dashed")
        }
        
        mod.pca.amtec.1+mod.pca.valid.1 
        
    })
    
    out.file <- paste0("figures/tcga_amtec_all_pca_delta", file.suf ,".pdf")
    
    ggsave(wrap_plots(bx.list, nrow=2, guides="collect"), file=out.file, width=14, height=6)
    
    out.file
    
}

immune.rep.by.me <- function(mes, btcr.res, resp.cols, rna.annot){
    
    btcr.res <- merge(btcr.res, rna.annot[,.(ptid, Biopsy=biopsy, cohort, r_acc, short_resp)], by=c("ptid", "Biopsy"))
    
    btcr.res.me <- merge(btcr.res, mes[cur_label == "Mod3 (brown)",.(r_acc=ptid, PC1)], by="r_acc")
    
    #Correlate read abundance with brown ME
    
    me3.cors <- btcr.res.me[cohort == "AMTEC" & type %in% c("BCR"),broom::tidy(cor.test(x=PC1, y=log2(Abundance+1))),by=.(Biopsy, `#chain`)]
    
    me3.cors <- .add.p.dots(me3.cors)
    
    me.brc.plot <- ggplot(data=btcr.res.me[cohort == "AMTEC" & type %in% c("BCR")], 
                          mapping=aes(x=PC1, y=log2(Abundance+1))) +
        geom_point(size=3, shape=21, mapping=aes(fill=short_resp), color="black") +
        scale_fill_manual(values=resp.cols, name="Best Response") +
        geom_smooth(method="lm", formula=y~x, se=F, inherit.aes = F, mapping=aes(x=PC1, y=log2(Abundance+1)), alpha=.75) +
        geom_text(data=me3.cors, mapping=aes(x=0, y=20, label=paste0("r: ", round(estimate, digits=3), dots)), inherit.aes = F) +
        facet_grid(Biopsy~`#chain`) +
        theme_bw(base_size=10) + xlab("Mod3 (brown) Eigengene") + ylab("Abundance")
    
    ggsave(me.brc.plot, file="figures/me_by_bcr_v1.pdf",  width=7, height=3)
    
    openxlsx::write.xlsx(list(res=me3.cors, ns=btcr.res.me[,.N,by=.(Biopsy, `#chain`)]), file="output/me_by_bcr_v1_stats.xlsx")
    
    "figures/me_by_bcr_v1.pdf"
}

ighg.diversity.plot <- function(btcr.res, rna.annot, resp.cols){
    
    btcr.res <- merge(btcr.res, rna.annot[,.(ptid, Biopsy=biopsy, cohort, short_resp)], by=c("ptid", "Biopsy"))
    
    #warnings are type conversions
    melt.btcr <- suppressWarnings(melt(data=btcr.res, id.vars=c("ptid", "Biopsy", "short_resp", "cohort", "#chain", "Abundance"), measure.vars=c("Entropy", "D50", "Gini", "Evenness")))
    
    div.p <- ggplot(
        data=melt.btcr[`#chain` %in% c("IGHG") & cohort == "AMTEC"],
        mapping = aes(x = ifelse(is.na(value), -Inf, value) , y=log2(Abundance+1))) + 
        geom_point(size=3, shape=21, mapping=aes(fill=short_resp, alpha=is.na(value))) +
        #geom_smooth(se=F, formula=y~x, method="loess", alpha=.25) +
        scale_fill_manual(values=resp.cols, name="Best Response") +
        scale_alpha_manual(values=c(`TRUE`=.5, `FALSE`=1), name="Undefined") +
        facet_grid(Biopsy~variable, scales="free") +
        theme_bw(base_size=10) + xlab("") + ylab("Abundance") +
        ggtitle("IGHG")
    
    ggsave(div.p, file="figures/ighg_div_v1.pdf", width=7.5, height=3)
    
    "figures/ighg_div_v1.pdf"
    
}

ighg.eveness.plot <- function(btcr.res, rna.annot, resp.cols, suf="", title=""){
    
    btcr.res <- merge(btcr.res, rna.annot[,.(ptid, Biopsy=biopsy, cohort, short_resp)], by=c("ptid", "Biopsy"))
    
    ighg.subset <- btcr.res[`#chain` == "IGHG" & cohort == "AMTEC"]
    
    ighg.subset.med <- ighg.subset[,.(med_evenness=median(Evenness)),by=.(ptid, Biopsy, short_resp)]
    
    ighg.subset.med[,bin_resp:=factor(ifelse(short_resp == "PD", "PD", "SDPR"))]
    
    test.sum <- ighg.subset.med[,broom::tidy(wilcox.test(med_evenness~bin_resp)), by=.(Biopsy)]
    
    test.sum <- .add.p.dots(test.sum)
    
    ss.plot <- ggplot(data=ighg.subset.med[is.na(med_evenness)==F], mapping=aes(x=short_resp, y=med_evenness)) +
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(height=0, width=.1, mapping=aes(fill=short_resp), shape=21, size=3) +
        scale_fill_manual(values=resp.cols, name="Best Response") +
        annotate("segment", x=c(1, 1, 2, 2, 3, 2.5), xend=c(1, 2.5, 3, 2, 3, 2.5), y=c(.975, 1, .96875, .96875, .96875, 1),  yend=c(1, 1, .96875, .96875-.025, .96875-.025, .96875)) +
        geom_text(data=test.sum, mapping=aes(x=2, y=1.1, label=dots, size=dots), inherit.aes=F) +
        scale_size_manual(values=c(`*`=10, `**`=10, `***`=10, `N.S.`=5), guide="none") +
        facet_wrap(~Biopsy) +
        ylim(c(0, 1.15)) +
        theme_bw() + xlab("") + ylab("Evenness") + ggtitle(title)
    
    ggsave(ss.plot, file=paste0("figures/ighg_evenness", suf,".pdf"), width=6, height=2.5)
    
    openxlsx::write.xlsx(list(res=test.sum, ns=ighg.subset.med[,.N,by=Biopsy]), file=paste0("output/ighg_evenness", suf,"_stats.xlsx"))
    
    paste0("figures/ighg_evenness", suf,".pdf")
    
}

tumor.purity.by.response <- function(bb.calls, resp.cols){
    
    bb.calls[,bin_resp:=factor(ifelse(Best_Response == "PD", "PD", "SDPR"))]
    
    use.calls <- bb.calls[cohort == "AMTEC" & Consensus_Call != "LAR"]
    
    .get.preds <- function(dt){
        lfit <- glm(bin_resp~as.numeric(Percent_Tumor), data=dt, family=binomial(link="logit"))
        cbind(dt, 
              response=factor(ifelse(predict(lfit, type="response") >= .5, "pSDPR", "pPD"))
        )
        
    }
    
    pred.calls <- use.calls[,.get.preds(.SD), by=Biopsy]
    
    mod.res <- use.calls[,broom::tidy(glm(bin_resp~as.numeric(Percent_Tumor), family=binomial(link="logit"))), by=Biopsy]
    
    #one dimensional decision boundary
    dec.bound <- mod.res[,.(xint=-estimate[1]/estimate[2]),by=Biopsy]
    
    mod.res <- mod.res[term != "(Intercept)"]
    
    mod.res <- .add.p.dots(mod.res)
    
    p1 <- ggplot(data=pred.calls, mapping=aes(x=as.numeric(Percent_Tumor), y=bin_resp)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size=3, mapping=aes(color=response, shape=Short_Location), height=.1, width=0) +
        facet_grid(Biopsy~.) +
        scale_shape_discrete(name="") +
        scale_color_manual(values=c(pPD="blue", pSDPR="red"), name="Predicted Response") +
        geom_text(data=mod.res, mapping=aes(y=2.35, x=85, label=dots)) +
        geom_vline(data=dec.bound, mapping=aes(xintercept=xint), linetype="dashed") +
        theme_bw(base_size=10) + 
        ylab("") + xlab("Percent Tumor")
    
    ggsave(p1, file="figures/amtec_tumor_purity_by_resp.pdf", width=4, height=3.5)
    
    openxlsx::write.xlsx(list(res=mod.res, ns=use.calls[,.N,by=Biopsy]), file="output/amtec_tumor_purity_by_resp.xlsx")
    
    "figures/amtec_tumor_purity_by_resp.pdf"
    
}

tumor.purity.tissue.by.mod.plot <- function(mes, bb.calls, amtec.samps, resp.cols){
    
    mes <- merge(amtec.samps, mes, by.x="r_acc", by.y="ptid", all.x=T)
    
    bb.mes <- merge(bb.calls, mes, by.x=c("De-Identified.Id", "Biopsy", "cohort"), by.y=c("ptid", "Biopsy", "cohort"))
    
    bb.mes <- bb.mes[cur_label %in% c("Mod1 (turquoise)", "Mod2 (blue)", "Mod3 (brown)")]
    
    me.cors <- bb.mes[cohort == "AMTEC" & Consensus_Call != "LAR",broom::tidy(cor.test(y=PC1, x=as.numeric(Percent_Tumor))),by=.(Biopsy, cur_label)]
    
    me.cors <- .add.p.dots(me.cors)
    
    p1 <- ggplot(data=bb.mes[cohort == "AMTEC" & Consensus_Call != "LAR"], mapping=aes(y=PC1, x=as.numeric(Percent_Tumor))) + 
        geom_point(size=3, mapping=aes(color=Best_Response, shape=Short_Location)) +
        facet_grid(Biopsy~cur_label) +
        geom_text(data=me.cors, inherit.aes=F, mapping=aes(x=60, y=25, label=dots, size=dots)) +
        geom_smooth(formula=y~x, method="lm", se=F, alpha=.5) +
        scale_color_manual(values=resp.cols, name="") +
        scale_size_manual(values=c(`*`=10, `**`=10, `***`=10, `N.S.`=5), guide="none") +
        scale_shape_discrete(name="") +
        ylim(c(-20, 30)) +
        theme_bw(base_size=10) + ylab("Module Eigengene") + xlab("Percent Tumor")
    
    ggsave(p1, file="figures/amtec_tumor_purity_by_me.pdf", width=8, height=3.5)
    
    openxlsx::write.xlsx(list(res=me.cors, ns=bb.mes[cohort == "AMTEC" & Consensus_Call != "LAR",.N,by=.(Biopsy, cur_label)]), file="output/amtec_tumor_purity_by_me.xlsx")
    
    "figures/amtec_tumor_purity_by_me.pdf"
    
}

tissue.by.mod.plot <- function(mes, bb.calls, amtec.samps, resp.cols){
    
    mes <- merge(amtec.samps, mes, by.x="r_acc", by.y="ptid", all.x=T)
    
    bb.mes <- merge(bb.calls, mes, by.x=c("De-Identified.Id", "Biopsy", "cohort"), by.y=c("ptid", "Biopsy", "cohort"))
    
    bb.mes <- bb.mes[cur_label %in% c("Mod1 (turquoise)", "Mod2 (blue)", "Mod3 (brown)")]
    
    bb.mes[, location_fac:=factor(Short_Location, levels=c("Breast", "Liver", "Lung", "Soft Tissue", "Lymph node"), ordered=T)]
    
    bb.mes[,is_ln:=factor(ifelse(Short_Location == "Lymph node", "LN", "Other"))]
    
    ln.assoc <- bb.mes[cohort == "AMTEC" & Consensus_Call != "LAR",broom::tidy(t.test(PC1~is_ln)),by=.(Biopsy, cur_label)]
    
    ln.assoc <- .add.p.dots(ln.assoc)
    
    p1 <- ggplot(data=bb.mes[cohort == "AMTEC" & Consensus_Call != "LAR"], mapping=aes(y=PC1, x=location_fac)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size=3, mapping=aes(color=short_resp, shape=Short_Location), height=0, width=.1) +
        facet_grid(Biopsy~cur_label) +
        geom_text(data=ln.assoc, inherit.aes=F, mapping=aes(x=3, y=30, label=dots), size=5) +
        scale_color_manual(values=resp.cols, name="") +
        scale_shape_discrete(name="") +
        annotate("segment", x=c(1, 4, 1, 3, 3, 5), xend=c(1, 4, 4, 3, 5, 5), y=c(20, 20, 22, 22, 24, 20),  yend=c(22, 22, 22, 24, 24, 24)) +
        theme_bw() + 
        ylim(c(-20, 35)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        xlab("") + ylab("Module Eigengene")
    
    ggsave(p1, file="figures/amtec_tissue_by_me.pdf", width=8, height=3.5)
    
    openxlsx::write.xlsx(list(res=ln.assoc, ns=bb.mes[cohort == "AMTEC" & Consensus_Call != "LAR",.N,by=.(Biopsy, cur_label)]), file="output/amtec_tissue_by_me.xlsx")
    
    "figures/amtec_tissue_by_me.pdf"
    
}

me.bb.heatmap <- function(mes, blia.blis, amtec.samps, bb.cols, resp.cols){
    
    bb.cols <- append(bb.cols, c(INDETERMINATE="white"))
    
    blia.blis <- blia.blis[(is.na(Consensus_Call)) | (Consensus_Call != "LAR")]
    
    blia.blis <- blia.blis[cohort == "AMTEC"]
    
    mes <- merge(amtec.samps, mes, by.x="r_acc", by.y="ptid", all.x=T)
    
    use.mes <- merge(blia.blis, mes, by.x=c("De-Identified.Id", "Biopsy", "cohort"), by.y=c("ptid", "Biopsy", "cohort"))
    
    use.mes <- use.mes[cohort == "AMTEC" & cur_label %in% c("Mod1 (turquoise)", "Mod2 (blue)", "Mod3 (brown)")]
    
    biopsy.list <- lapply(c("Bx1", "Bx2"), function(tmp.biopsy){
        
        pts.ord <- blia.blis[Biopsy == tmp.biopsy][order(factor(Consensus_Call, levels=c("BLIS", "INDETERMINATE", "BLIA"), ordered=T), 
                                                         factor(Best_Response, levels=c("PD", "SD", "PR"), ordered=T))]
        
        use.mes[,pt_fac:=factor(`De-Identified.Id`, levels=pts.ord$`De-Identified.Id`, ordered=T)]
        
        p1 <- ggplot(data=use.mes[Biopsy == tmp.biopsy], mapping=aes(x=pt_fac, y=cur_label, fill=PC1)) +
            geom_tile(color="black") +
            theme_classic() +
            geom_hline(yintercept=0, linetype="dashed") +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.x=element_blank()) +
            ggtitle(tmp.biopsy) +
            ylab("")
        
        if (tmp.biopsy == "Bx1"){
            p1 <- p1 + scale_fill_gradient2(low=bb.cols["BLIS"], mid=bb.cols["INDETERMINATE"], high=bb.cols["BLIA"])
        }else{
            p1 <- p1 + scale_fill_gradient2(low=bb.cols["BLIS"], mid=bb.cols["INDETERMINATE"], high=bb.cols["BLIA"], guide="none")
        }
        
        blia.blis[,pt_fac:=factor(`De-Identified.Id`, levels=pts.ord$`De-Identified.Id`, ordered=T)]
        blia.blis[,cur_label:="Burstein"]
        
        p2 <- ggplot(data=blia.blis[Biopsy == tmp.biopsy], mapping=aes(x=pt_fac, y=cur_label, fill=Consensus_Call)) +
            geom_tile(color="black") +
            scale_fill_manual(values=bb.cols) +
            theme_classic() +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y = element_blank())
        
        blia.blis[,cur_label:="Best Response"]
        
        p3 <- ggplot(data=blia.blis[Biopsy == tmp.biopsy], mapping=aes(x=pt_fac, y=cur_label, fill=Best_Response)) +
            geom_tile(color="black") +
            scale_fill_manual(values=resp.cols) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  axis.ticks.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())
        
        (p1 / p2 / p3) + plot_layout(heights=c(3/5, 1/5, 1/5))
        
    })
    
    ggsave(wrap_plots(biopsy.list) + plot_layout(guides="collect"), file="figures/burstein_me_comp.pdf", width=15, height=4.5)
    
    "figures/burstein_me_comp.pdf"
}

m3s.exprs.plot <- function(comb.mes,comb.lcpm, main.samps, sc.exp, mm, m3s, m3s.meta){
    
    #TCGA
    
    main.exprs <- data.table(reshape2::melt(comb.lcpm, as.is=T))
    colnames(main.exprs) <- c("gene_id", "r_acc", "exprs")
    
    main.exprs <- merge(main.samps, main.exprs, by=c("r_acc"))
    
    tcga.exprs <- main.exprs[cohort == "TCGA"]
    
    tcga.exprs[,`:=`(tcga_mean=mean(exprs), tcga_sd=sd(exprs)),by=gene_id]
    
    tcga.exprs[,z:=(exprs-tcga_mean)/(tcga_sd)]
    
    #Combine with AMTEC
    
    all.exprs <- rbind(
        tcga.exprs,
        sc.exp[biopsy=="Bx1",names(tcga.exprs),with=F]
    )
    
    brown.mes <- comb.mes[cur_label == "Mod3 (brown)"]
    
    base.col <- circlize::colorRamp2(seq(-6, 6, length = 50), plasma(50, end=.9))
    
    colfun <- circlize::colorRamp2(seq(-2, 2, length = 50), plasma(50))
    
    brown.mes <- brown.mes[cohort == "TCGA"][order(PC1),]
    
    z.mat <- reshape2::acast(gene_id~r_acc, value.var="z",data=all.exprs[cohort == "TCGA"])
    
    m3s.meta <- m3s.meta[order(-kme),]
    
    m3s.meta.df <- as.data.frame(t(z.mat[m3s.meta$gene_id,]))
    names(m3s.meta.df) <- m3s.meta$gene_name
    
    sub.ha <- HeatmapAnnotation(df=m3s.meta.df[brown.mes$ptid,],
                                PC1 = anno_barplot(brown.mes$PC1, gp=gpar(col=base.col(brown.mes$PC1), fill=base.col(brown.mes$PC1))),
                                show_legend = c(F, F, F, F), height=unit(1, "in"), col=lapply(m3s.meta.df, function(x) colfun),
                                simple_anno_size = unit(.1, "in"), annotation_name_gp = gpar(cex=.75, fontface=c("plain", "bold", "plain", "plain")))
    
    gene.kme <- mm[cur_label == "Mod3 (brown)"][order(-kme)]
    
    gene.kme[,geneset:=""]
    
    gene.kme[gene_id %in% m3s$gene_id,geneset:="mTNBC3s"]
    
    gene.kme[,geneset:=factor(geneset, levels=c("mTNBC3s", ""), ordered=T)]
    
    hr <- rowAnnotation(kME=anno_barplot(gene.kme$kme, gp=gpar(col="grey")))
    
    ht <- Heatmap(z.mat[gene.kme$gene_id,brown.mes$ptid], col=colfun, top_annotation=sub.ha, right_annotation=hr,
                  show_column_names=F, cluster_columns=F, cluster_rows=F, show_row_names=F, row_names_side="left",
                  row_title_rot=0, show_row_dend=F, split=gene.kme$geneset,
                  heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
                                              legend_width = unit(5, "cm"), title_position = "leftcenter",
                                              title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                              labels_gp = gpar(fontsize = 8)), 
                  name="Scaled Exprs", border=T, column_title = "Mod3 (Brown)")
    
    
    pdf(file="figures/mzb1_brown_eigengene.pdf", width=6, height=3)
    
    draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T, padding = unit(c(2, 11, 2, 2), "mm"))
    
    decorate_annotation("IGKC", { 
        grid.text("mTNBC3s_top3", x=unit(0, "npc")-unit(.33, "npc"), y = unit(.5, "npc") , just = "left", gp=gpar(fontsize = 13.2)) 
    })
    
    dev.off()
    
    "figures/mzb1_brown_eigengene.pdf"
    
}

plot.summary.matrix <- function(perf.list, resp.cols){
    
    melt.bb <- perf.list$data
    perf.summary <- perf.list$summary
    
    perf.summary <- perf.summary[cohort == "AMTEC" & Biopsy == "Bx1" & variable %in% c("Best_Response", "BLIA_BLIS_IHC_ssGSEA" , "BLIA_BLIS_IHC", "Zhang_pB","MZB1", "ssGSEA_Top3","ssGSEA",  "Brown_ME")]
    melt.bb <- melt.bb[cohort == "AMTEC" & Biopsy == "Bx1" & variable %in% c("Best_Response", "BLIA_BLIS_IHC_ssGSEA" , "BLIA_BLIS_IHC", "Zhang_pB","MZB1", "ssGSEA_Top3","ssGSEA",  "Brown_ME")]
    
    br.lets <- melt.bb[variable == "Best_Response"]
    br.lets[,code:=substr(Short_Location, 1, 2)]
    br.lets[is.na(code), code:=""]
    
    p1 <- ggplot(data=melt.bb, mapping=aes(y=var_fac, x=pat_fac, fill=value, alpha=value)) + geom_tile(color="black") +
        geom_text(data=br.lets, mapping=aes(label=code)) +
        geom_text(data=perf.summary, mapping=aes(x=Inf, label=paste(scales::percent(Accuracy, accuracy=.1), nc_rate)), hjust = 0) +
        coord_cartesian(clip="off") +
        facet_grid(Biopsy+br_cats~cohort, scales="free", space="free") + xlab("") + ylab("") +
        scale_fill_manual(values=c(BLIA="#F8766D", BLIS="#00BFC4", INDETERMINATE="white", LAR="purple", pPD="blue", pSDPR="red", NC="grey", resp.cols), na.value="grey", guide="none") +
        scale_alpha_manual(values=c(BLIA=1, BLIS=1, INDETERMINATE=1, LAR=1, pPD=.75, pSDPR=.75, PD=1, SD=1, PR=1, NC=1), guide="none",na.value=1) +
        theme_classic() + 
        theme(strip.background.y=element_blank(), strip.text.y=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(panel.spacing.x= unit(8, "lines"), plot.margin = unit(c(0,8,0,0), "lines"))
    
    #margin: top, right, bottom, and left
    
    ggsave(p1, file="figures/base_summary_fig.pdf", width=8.5, height=4)
    
    "figures/base_summary_fig.pdf"
}

de.by.net.conc.biopsy <- function(comb.exprs.de, net.conc.de){
    
    comb.net.exprs <- merge(net.conc.de, comb.exprs.de[,.(gene_id, biopsy, exprs_t=t, exprs_fdr=adj.P.Val)], by=c("gene_id", "biopsy"))
    
    comb.net.exprs[,is_sig:=(abs(t) >= 2 & abs(exprs_t) >= 2) | (abs(t) >= 2 & abs(exprs_t) < 2)]
    
    text.dt <- data.table(x=c(-6, -6, 0, 0, 3.25, 3.25),y=c(5, -4, 5, -4, 5, -4), label=c(1, 4, 2, 5, 3, 6))
    
    p1 <- ggplot(data=comb.net.exprs, mapping=aes(x=exprs_t, y=t, color=cur_label)) +
        geom_point(mapping=aes(alpha=is_sig)) +
        scale_color_manual(values=c("Mod1 (turquoise)"="turquoise", "Mod2 (blue)"="blue", "Mod3 (brown)"="brown"), name="Module") +
        scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.1), name="", guide="none") +
        geom_hline(yintercept=c(-2, 2), linetype="dashed") +
        geom_vline(xintercept=c(-2, 2), linetype="dashed") +
        geom_text(data=text.dt, mapping=aes(x=x, y=y, label=label), inherit.aes=F, size=5) +
        geom_text_repel(data=comb.net.exprs[network_concept=="Connectivity" & biopsy == "Bx1" & gene_name == "IFI27"],
                        mapping=aes(label=gene_name), min.segment.length=0, box.padding = 0.5, color="black", size=4,
                        ylim=c(3, Inf), xlim=c(-Inf, -2)) +
        facet_grid(biopsy~network_concept) +
        theme_bw() + xlab("Expression T-statistic") + ylab("Network Concept T-statistic")
    
    ggsave(p1, file="figures/network_concepts_by_exprs.png", width=10, height=6)
    
    "figures/network_concepts_by_exprs.png"
}

ifi27.de.by.net.conc.biopsy  <- function(comb.exprs.de, net.conc.de){
    
    comb.net.exprs <- merge(net.conc.de, comb.exprs.de[,.(gene_id, biopsy, exprs_t=t, exprs_fdr=adj.P.Val)], by=c("gene_id", "biopsy"))
    
    comb.net.exprs[,is_sig:=(abs(t) >= 2 & abs(exprs_t) >= 2) | (abs(t) >= 2 & abs(exprs_t) < 2)]
    
    text.dt <- data.table(x=c(-4, -4, 0, 0, 3.25, 3.25),y=c(5, -4, 5, -4, 5, -4), label=c(1, 4, 2, 5, 3, 6))
    
    p1 <- ggplot(data=comb.net.exprs[network_concept=="Connectivity" & biopsy == "Bx1"], mapping=aes(x=exprs_t, y=t, color=cur_label)) +
        geom_point(mapping=aes(alpha=is_sig)) +
        scale_color_manual(values=c("Mod1 (turquoise)"="turquoise", "Mod2 (blue)"="blue", "Mod3 (brown)"="brown"), name="Module") +
        scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.1), name="", guide="none") +
        geom_hline(yintercept=c(-2, 2), linetype="dashed") +
        geom_vline(xintercept=c(-2, 2), linetype="dashed") +
        geom_text(data=text.dt, mapping=aes(x=x, y=y, label=label), inherit.aes=F, size=5) +
        geom_text_repel(data=comb.net.exprs[network_concept=="Connectivity" & biopsy == "Bx1" & gene_name == "IFI27"],
                        mapping=aes(label=gene_name), min.segment.length=0, box.padding = 0.5, color="black", size=4,
                        ylim=c(3, Inf), xlim=c(-Inf, -2)) +
        facet_grid(biopsy~network_concept) +
        theme_bw() + xlab("Expression T-statistic") + ylab("Network Concept T-statistic")
    
    ggsave(p1, file="figures/ifi27_network_concepts_by_exprs.pdf", width=5, height=3.5)
    
    "figures/ifi27_network_concepts_by_exprs.pdf"
    
}

de.by.net.conc.delta <- function(res.dt, net.conc.de){
    
    comb.net.exprs <- merge(net.conc.de, res.dt[,.(gene_id, exprs_t=t)], by=c("gene_id"))
    
    comb.net.exprs[,is_sig:=(abs(t) >= 2 & abs(exprs_t) >= 2) | (abs(t) >= 2 & abs(exprs_t) < 2)]
    
    text.dt <- data.table(x=c(-4, -4, 0, 0, 3.25, 3.25),y=c(4.5, -4, 4.5, -4, 4.5, -4), label=c(1, 4, 2, 5, 3, 6))
    
    p1 <- ggplot(data=comb.net.exprs, mapping=aes(x=exprs_t, y=t, color=cur_label)) +
        geom_point(mapping=aes(alpha=is_sig)) +
        scale_color_manual(values=c("Mod1 (turquoise)"="turquoise", "Mod2 (blue)"="blue", "Mod3 (brown)"="brown"), name="Module") +
        scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.1), name="", guide="none") +
        geom_hline(yintercept=c(-2, 2), linetype="dashed") +
        geom_vline(xintercept=c(-2, 2), linetype="dashed") +
        geom_text(data=text.dt, mapping=aes(x=x, y=y, label=label), inherit.aes=F, size=5) +
        geom_text_repel(data=comb.net.exprs[network_concept=="MAR" & gene_name == "KRT23"],
                        mapping=aes(label=gene_name), min.segment.length=0, box.padding = 0.5, color="black", size=4,
                        ylim=c(3, Inf), xlim=c(-Inf, -2)) +
        facet_grid(biopsy~network_concept) +
        theme_bw() + xlab("Expression T-statistic") + ylab("Network Concept T-statistic")
    
    ggsave(p1, file="figures/network_concepts_by_exprs_paired_diff.png", width=10, height=3)
    
    "figures/network_concepts_by_exprs_paired_diff.png"
}

krt23.de.by.net.conc.delta <- function(res.dt, net.conc.de){
    
    comb.net.exprs <- merge(net.conc.de, res.dt[,.(gene_id, exprs_t=t)], by=c("gene_id"))
    
    comb.net.exprs[,is_sig:=(abs(t) >= 2 & abs(exprs_t) >= 2) | (abs(t) >= 2 & abs(exprs_t) < 2)]
    
    text.dt <- data.table(x=c(-4, -4, 0, 0, 3.25, 3.25),y=c(4.5, -4, 4.5, -4, 4.5, -4), label=c(1, 4, 2, 5, 3, 6))
    
    p1 <- ggplot(data=comb.net.exprs[network_concept=="MAR"], mapping=aes(x=exprs_t, y=t, color=cur_label)) +
        geom_point(mapping=aes(alpha=is_sig)) +
        scale_color_manual(values=c("Mod1 (turquoise)"="turquoise", "Mod2 (blue)"="blue", "Mod3 (brown)"="brown"), name="Module") +
        scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.1), name="", guide="none") +
        geom_hline(yintercept=c(-2, 2), linetype="dashed") +
        geom_vline(xintercept=c(-2, 2), linetype="dashed") +
        geom_text(data=text.dt, mapping=aes(x=x, y=y, label=label), inherit.aes=F, size=5) +
        geom_text_repel(data=comb.net.exprs[network_concept=="MAR" & gene_name == "KRT23"],
                        mapping=aes(label=gene_name), min.segment.length=0, box.padding = 0.5, color="black", size=4,
                        ylim=c(3, Inf), xlim=c(-Inf, -2)) +
        facet_grid(biopsy~network_concept) +
        theme_bw() + xlab("Expression T-statistic") + ylab("Network Concept T-statistic")
    
    ggsave(p1, file="figures/krt23_network_concepts_by_exprs_paired_diff.pdf", width=5, height=3.5)
    
    "figures/krt23_network_concepts_by_exprs_paired_diff.pdf"
}

make.lioness.mats <- function(mm, lion.conc, lion.ar, all.rna){
    
    use.map <- mm[cur_label == "Mod2 (blue)"]
    
    bx1.pd.ar <- lion.ar[use.map$gene_id,use.map$gene_id,
                         all.rna[cohort == "AMTEC" & biopsy == "Bx1" & short_resp == "PD",paste(ptid, biopsy, sep=".")]]
    
    dimnames(bx1.pd.ar)[1:2] <- list(use.map$gene_name, use.map$gene_name)
    
    .trans.mat <- function(mat){
        scaled.mat <- mat / max(mat)
        diag(scaled.mat) <- 1
        scaled.mat
    }
    
    bx1.pd.sum <- .trans.mat(bx1.pd.ar[,,1])
    
    for(i in 2:dim(bx1.pd.ar)[3]){
        bx1.pd.sum <- bx1.pd.sum + .trans.mat(bx1.pd.ar[,,i])
    }
    
    bx1.pd.sum <- bx1.pd.sum / dim(bx1.pd.ar)[3]
    
    use.cols <- circlize::colorRamp2(breaks=c(0, .25, .5), colors=viridis::mako(n=3, alpha = 1, begin = 0, end = 1, direction = -1))
    
    hc = columnAnnotation(foo = anno_mark(at = which(colnames(bx1.pd.sum) == "IFI27"), side="bottom",labels = "IFI27"))
    
    hr <- rowAnnotation(IFI27=anno_simple(bx1.pd.sum[,"IFI27"], col=use.cols))
    
    ht1 <- Heatmap(bx1.pd.sum, col=use.cols, show_column_names=F, show_row_names=F,
                   clustering_method_rows="average", clustering_method_columns="average",
                   name="Adjacency", column_title="PD Patients Bx1", #right_annotation=h1r,
                   heatmap_legend_param = list(direction = "horizontal"), show_row_dend=F, show_column_dend=F,
                   bottom_annotation=hc, right_annotation=hr)
    
    bx1.sdpr.ar <- lion.ar[use.map$gene_id,use.map$gene_id,
                           all.rna[cohort == "AMTEC" & biopsy == "Bx1" & short_resp != "PD",paste(ptid, biopsy, sep=".")]]
    
    dimnames(bx1.sdpr.ar )[1:2] <- list(use.map$gene_name, use.map$gene_name)
    
    bx1.sdpr.sum <- .trans.mat(bx1.sdpr.ar[,,1])
    
    for(i in 2:dim(bx1.sdpr.ar)[3]){
        bx1.sdpr.sum <- bx1.sdpr.sum + .trans.mat(bx1.sdpr.ar[,,i])
    }
    
    bx1.sdpr.sum <- bx1.sdpr.sum / dim(bx1.sdpr.ar)[3]
    
    hc = columnAnnotation(foo = anno_mark(at = which(colnames(bx1.sdpr.sum) == "IFI27"), side="bottom",labels = "IFI27"))
    
    hr <- rowAnnotation(IFI27=anno_simple(bx1.sdpr.sum[,"IFI27"], col=use.cols))
    
    ht2 <- Heatmap(bx1.sdpr.sum, col=use.cols, show_column_names=F, show_row_names=F,
                   clustering_method_rows="average", clustering_method_columns="average",
                   name="Adjacency", column_title="SD/PR Patients Bx1", #right_annotation=h1r,
                   heatmap_legend_param = list(direction = "horizontal"), show_row_dend=F, show_column_dend=F,
                   bottom_annotation=hc, right_annotation=hr)
    
    png(file="figures/pd_lioness_diff.png", width=3, height=4, units="in", res=300)
    
    draw(ht1, heatmap_legend_side = "bottom")
    
    dev.off()
    
    png(file="figures/sdpr_lioness_diff.png", width=3, height=4, units="in", res=300)
    
    draw(ht2, heatmap_legend_side = "bottom")
    
    dev.off()
    
    "figures/pd_lioness_diff.png"
}

plot.gene.net.results <- function(mm, best.gn, all.samps, valid.samps, lion, valid.lion, preds, exprs, resp.cols){
    
    #unfortunately amtec is defined by ptid_biop not r_acc...
    all.samps[,ptid_biop:=paste(ptid, biopsy, sep=".")]
    all.samps <- all.samps[cohort == "AMTEC"]
    
    valid.samps[,ptid_biop:=paste(ptid, biopsy, sep=".")]
    valid.samps <- valid.samps[cohort == "Validation"]
    
    lion.dt <- rbindlist(lapply(lion$module.stats[[best.gn$cur_label]], function(x) x$per_gene), idcol="ptid_biop")
    lion.dt <- merge(all.samps[,.(r_acc, ptid_biop, cohort, biopsy, short_resp)], lion.dt, by="ptid_biop")
    
    valid.lion.dt <- rbindlist(lapply(valid.lion$module.stats[[best.gn$cur_label]], function(x) x$per_gene), idcol="ptid_biop")
    valid.lion.dt <- merge(valid.samps[,.(r_acc, ptid_biop, cohort, biopsy, short_resp)], valid.lion.dt, by=c("ptid_biop"))
    
    comb.lion <- rbind(lion.dt[,names(valid.lion.dt),with=F], valid.lion.dt)
    
    comb.lion <- merge(comb.lion, mm[,.(gene_id, gene_name)], by="gene_id")
    
    comb.lion <- merge(comb.lion, exprs[,.(gene_id, r_acc, z)], by=c("gene_id", "r_acc"))
    
    comb.lion <- comb.lion[gene_id == best.gn$gene_id]
    
    #plot the cutpoint as well
    
    best.fit <- preds[[best.gn$desc]]$fit
    
    y.break <- split_node(node_party(best.fit[[1]]))$breaks
    
    x.break <- split_node(node_party(best.fit[[3]]))$breaks
    
    comb.lion$value <- comb.lion[[best.gn$net_conc]]
    
    #label the extreme samples
    
    highlight.samps <- comb.lion[cohort == "AMTEC" & biopsy == best.gn$biopsy]
    highlight.samps <- highlight.samps[order(highlight.samps[[best.gn$net_conc]])][c(1, ceiling(.N/2), .N)]
    highlight.samps[,c("ptid","biopsy"):=tstrsplit(ptid_biop, "\\.")]
    
    p1 <- ggplot(data=comb.lion[biopsy == best.gn$biopsy], mapping=aes(x=z, y=log2(value+1), fill=short_resp)) + geom_point(shape=21, size=3) +
        scale_fill_manual(values=resp.cols, name="Best Response") + 
        facet_wrap(~biopsy+cohort) + 
        geom_hline(yintercept=y.break, linetype="dashed") +
        geom_vline(xintercept=x.break, linetype="dashed") +
        #geom_text_repel(data=highlight.samps, mapping=aes(label=ptid)) +
        theme_bw() + xlab(paste(comb.lion$gene_name[1], "Scaled Expression")) + ylab(paste0(comb.lion$gene_name[1], " log2(", best.gn$net_conc, ")"))
    
    ggsave(p1, file=paste0("figures/tmp_gene_net_perf_",comb.lion$gene_name[1], "_", best.gn$net_conc,".pdf"), width=5.5, height=3)
    
    paste0("figures/tmp_gene_net_perf_",comb.lion$gene_name[1], "_", best.gn$net_conc,".pdf")
}

plot.krt23.biopsy.mar <- function(preds, resp.cols){
    
    amtec.preds <- rbindlist(lapply(preds,function(x) x$preds))
    amtec.preds[,`:=`(cur_label=cur_labels, ptid=sapply(strsplit(ptid_biop, "\\."), "[[", 1))]
    
    cast.lion <- dcast(ptid+short_resp~biopsy, value.var="value", data=amtec.preds[gene_id == "ENSG00000108244.16" & variable == "MAR"])
    cast.lion[,cohort:="AMTEC"]
    
    p1 <- ggplot(data=cast.lion, mapping=aes(x=Bx1, y=Bx2, fill=short_resp)) + geom_point(shape=21, size=3) +
        scale_fill_manual(values=resp.cols, name="Best Response", guide="none") + 
        facet_wrap(~cohort) + 
        geom_abline(slope=1, intercept=0, linetype="dashed") +
        theme_bw() + xlab("KRT23 log2(MAR+1) Bx1") + ylab("KRT23 log2(MAR+1) Bx2")
    
    ggsave(p1, file="figures/net_conc_plot_KRT23_MAR.pdf", width=3, height=3)
    
    "figures/net_conc_plot_KRT23_MAR.pdf"
    
}

plot.krt23.delta.mar <- function(all.samps, valid.samps, preds, resp.cols){
    
    amtec.preds <- rbindlist(lapply(preds$amtec,function(x) x$preds))
    amtec.preds <- merge(unique(all.samps[,.(ptid, short_resp, cohort)]), amtec.preds, by="ptid")
    
    valid.preds <- rbindlist(lapply(preds$valid,function(x) x$preds))
    valid.preds <- merge(unique(valid.samps[,.(ptid, short_resp, cohort)]), valid.preds, by="ptid")
    
    comb.lion <- rbind(amtec.preds, valid.preds[,names(amtec.preds),with=F])
    
    best.fit <- preds$amtec[["ENSG00000108244.16.MAR.Mod1 (turquoise)"]]$fit
    
    y.break <- split_node(node_party(best.fit[[1]]))$breaks
    
    p1 <- ggplot(data=comb.lion[gene_id == "ENSG00000108244.16" & variable == "MAR"], mapping=aes(x=z, y=value, fill=short_resp)) + geom_point(shape=21, size=3) +
        scale_fill_manual(values=resp.cols, name="Best Response") + 
        facet_wrap(~cohort) + 
        geom_hline(yintercept=y.break, linetype="dashed") +
        theme_bw() + xlab("KRT23 Expression (Bx2-Bx1)") + ylab("KRT23 log2(MAR+1) (Bx2-Bx1)")
    
    ggsave(p1, file="figures/delta_net_conc_plot_KRT23_MAR.pdf", width=5.5, height=3)
    
    "figures/delta_net_conc_plot_KRT23_MAR.pdf"
}

