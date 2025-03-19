library(targets)
source("R/preprocess.R")
source("R/summarize.R")
source("R/figures.R")
source("R/validation.R")
tar_option_set(packages=c("data.table"))

list(
    
    #preprocessing
    
    ##colors
    
    tar_target(
        response_cols,
        setNames(RColorBrewer::brewer.pal(n=3, "Dark2"), c("SD", "PR", "PD"))
    ),
    
    tar_target(
        blia_blis_cols,
        setNames(c(scales::hue_pal()(2), "purple", scales::hue_pal()(2)), c("BLIA", "BLIS", "LAR", "Non-BLIS/LAR", "BLIS/LAR"))
    ),
    
    ##gene annotation
    
    tar_target(
        trans_info,
        parse.trans.info("data/annotation/gencode.v24.annotation.gtf.gz")
    ),
    
    tar_target(
        gene_info,
        gene.info.from.trans(trans_info)
    ),
    
    tar_target(
        blia_blis_cents,
        fread("data/burstein/TNBC_Ding_77_gene_signatures.txt")
    ),
    
    ##tcga
    
    tar_target(
        tpm_mat,
        read.matrix.from.file("data/tcga/TCGA_BRCA_tpm.tsv.gz")
    ),
    
    tar_target(
        counts_mat,
        read.matrix.from.file("data/tcga/TCGA_BRCA_counts.tsv.gz")
    ),
    
    tar_target(
        tcga_annots,
        process.subtypes("data/tcga/TCGA_ID_MAP.csv", "data/tcga/brca_tcga_pan_can_atlas_2018_clinical_data.tsv"),
    ),
    
    tar_target(
        filt_annots,
        filter.subtypes(tcga_annots)
    ),
    
    tar_target(
        tcga_dge,
        tcga.kallisto.abund(tpm_mat, counts_mat, tcga_annots),
        packages=c("tximport", "edgeR")
    ),
    
    ###expression filtering
    tar_target(
        filt_tcga_dge,
        filter.tcga.kallisto.tpm(tcga_dge, filt_annots),
        packages="edgeR"
    ),
    
    ##AMTEC main cohort
    
    tar_target(
        amtec_annot_file,
        "data/amtec/amtec_rna_mapping_deident.RData",
        format="file"
    ),
    
    tar_target(
        amtec_annots,
        local({get(load(amtec_annot_file))})
    ),
    
    tar_target(
        amtec_kallisto,
        local({mget(load("data/amtec/amtec_kallisto_counts.RData"))})
    ),
    
    tar_target(
        filt_amtec_dge,
        form.dge.kallisto(amtec_kallisto )
    ),
    
    ##Combining AMTEC and TCGA
    
    tar_target(
        amtec_tcga_annot,
        make.combined.annot(amtec_annots,filt_tcga_dge, filt_annots, "AMTEC")
    ),
    
    tar_target(
        comb_ltpm,
        combine.abundance(filt_tcga_dge, filt_amtec_dge, amtec_tcga_annot),
        packages="sva"
    ),
    
    ##AMTEC validation cohort
    
    tar_target(
        valid_samp_file,
        "data/amtec/validation_rna_mapping_deident.RData",
        format="file"
    ),
    
    tar_target(
        valid_exprs_file,
        "data/amtec/validation_kallisto_counts.RData",
        format="file"
    ),
    
    tar_target(
        valid_samps,
        local({get(load(valid_samp_file))})
    ),
    
    tar_target(
        valid_exprs,
        get.validation.exprs(valid_exprs_file)
    ),
    
    tar_target(
        valid_tcga_annot,
        make.combined.annot(valid_samps,filt_tcga_dge, filt_annots, "Validation")
    ),
    
    tar_target(
        valid_ltpm,
        combine.abundance(filt_tcga_dge, valid_exprs, valid_tcga_annot),
        packages="sva"
    ),
    
    tar_target(
        scaled_exprs,
        get.scaled.exprs(gene_info, comb_ltpm, amtec_tcga_annot, valid_ltpm, valid_tcga_annot)
    ),
    
    ##Burstein subtyping
    
    tar_target(
        comb_clin_bb_file,
        "data/amtec/amtec_manuscript_bb_calls_deident_v1.xlsx",
        format="file"
    ),
    
    tar_target(
        comb_clin_bb,
        data.table(openxlsx::read.xlsx(comb_clin_bb_file))
    ),
    
    tar_target(
        batch_blia_blis,
        batch.call.blia.blis(blia_blis_cents, comb_ltpm,amtec_tcga_annot , gene_info, "TCGA")
    ),
    
    ##Zhang et al 2021 genes
    
    tar_target(
        zhang_2021_sets,
        data.table(openxlsx::read.xlsx("data/zhang/bcell_wrkbk.xlsx"))[,.(signature, gene_name=genes)]
    ),
    
    tar_target(
        tamborero_2018_sets,
        compile.tamborero.sets("data/annotation/10780432ccr173509-sup-192911_3_supp_4675335_p6rxmz.xlsx")
    ),
    
    ##WGCNA
    
    ##these input files are derived from WGCNA.Rmd
    tar_target(
        mm,
        get.module.membership(gene_info, 
                              file.path("data", "wgcna", "run_twok_2p_v4.RData"), 
                              file.path("data", "wgcna", "used_expression.RData"), blia_blis_cents),
        packages="WGCNA"
    ),
    
    ##Table S1
    tar_target(
        written_mm,
        write.mm(mm),
        format="file"
    ),
    
    tar_target(
        mtnbc3s,
        mm[cur_label=="Mod3 (brown)" & kme > .9]
    ),
    
    tar_target(
        mtnbc3s_top3,
        mm[order(-kme)][1:3]
    ),
    
    tar_target(
        comb_mes,
        proc.mes(comb_ltpm, amtec_tcga_annot, mm, "AMTEC")
    ),
    
    #was valid_mods
    tar_target(
        valid_mes,
        proc.mes(valid_ltpm, valid_tcga_annot, mm, "Validation")
    ),
    
    ##LIONESS
    
    ###for processing on HPC--see lioness.Rmd
    tar_target(
        for_lioness,
        save(mm, comb_ltpm, amtec_tcga_annot, file="output/output_exprs_lioness.RData"),
        format="file"
    ),
    
    ###output from lioness.Rmd
    tar_target(
        lion_concepts,
        local({mget(load("data/wgcna/lioness_network_concepts.RData"))})
    ),
    
    ###output from lioness.Rmd
    tar_target(
        valid_lion_concepts,
        local({mget(load("data/wgcna/lioness_validation_network_concepts.RData"))})
    ),
    
    ###output from lioness.Rmd
    tar_target(
        lion_network,
        local({get(load("data/wgcna/lioness_array.RData"))})
    ),
    
    tar_target(
        tcga_bb_preds_output,
        fit.blia.blis.tcga(comb_mes, amtec_tcga_annot, batch_blia_blis$call),
        packages=c("tidymodels", "openxlsx"),
        format="file"
    ),
    
    
    #Summarization/model fitting
    
    ##Hallmark enrichment of modules
    tar_target(
        mod_enrich,
        enrich.modules(mm, "data/annotation/h.all.v7.5.1.symbols.gmt")
    ),
    
    ##Figure 2a
    
    tar_target(
        mod_enrich_plot,
        mm.hallmark.plot(mod_enrich),
        format="file",
        packages=c("ggplot2", "stringr")
    ),
    
    ##Figure 2b
    ###requires internal data/code--not provided
    
    ##Figure 2c
    tar_target(
        mods_by_bb_plot,
        tcga.bb.plot(comb_mes, batch_blia_blis$call, amtec_tcga_annot),
        packages="ggplot2",
        format = "file"
    ),
    
    ##Figure 2d
    tar_target(
        brown_kme_plot,
        brown.kme.figure(mm, amtec_tcga_annot, "data/annotation/Cell_marker_Human.xlsx"),
        format = "file",
        packages=c("ggplot2", "ggrepel")
    ),
    
    ##Figure 2e
    tar_target(
        brown_bcell_plot,
        brown.bcell.cors(tamborero_2018_sets, filt_amtec_dge, comb_mes, gene_info, amtec_tcga_annot, response_cols),
        format="file",
        packages=c("GSVA", "ggplot2")
    ),
    
    ##Figure 2f
    tar_target(
        brown_mihc_plot,
        brown.mihc.cors(comb_mes, amtec_tcga_annot, "data/amtec/Supplementary_Data_S6.xlsx"),
        format = "file",
        packages=c("ggplot2", "data.table")
    ),
    
    ##ctree fitting of MEs
    
    tar_target(
        me_preds,
        ctree.mes(comb_mes, amtec_tcga_annot),
        packages=c("caret", "partykit", "yardstick")
    ),
    
    ##Figure S1
    
    tar_target(
        all_me_resp_fig,
        plot.all.modules.by.response(amtec_tcga_annot, valid_tcga_annot, comb_mes,  valid_mes, response_cols, me_preds),
        format = "file",
        packages="patchwork"
    ),
    
    ##B cell repetoire analysis
    
    tar_target(
        imm_rep_sum,
        fread("data/amtec/amtec_trust4_summary.txt", colClasses=c(ptid="character"))
    ),
    
    tar_target(
        imm_rep_subs,
        fread("data/amtec/amtec_trust4_summary-subsampled.txt", colClasses=c(ptid="character"))
    ),
    
    ##Figure S2a
    tar_target(
        me_by_imm_rep_plot,
        immune.rep.by.me(comb_mes, imm_rep_sum, response_cols, amtec_tcga_annot),
        format = "file",
        packages=c("ggplot2")
    ),
    
    ##Figure S2b
    tar_target(
        ighg_diversity_plot,
        ighg.diversity.plot(imm_rep_sum, amtec_tcga_annot, response_cols),
        format = "file",
        packages=c("ggplot2")
    ),
    
    ##Figure S2c
    tar_target(
        raw_even_plot,
        ighg.eveness.plot(imm_rep_sum, amtec_tcga_annot, response_cols, "", "Raw"),
        format = "file",
        packages=c("ggplot2")
    ),
    
    ##Figure S2d
    tar_target(
        ssubs_even_plot,
        ighg.eveness.plot(imm_rep_subs, amtec_tcga_annot, response_cols, "_subs", "Downsampled"),
        format = "file",
        packages=c("ggplot2")
    ),
    
    ##Figure S3a
    tar_target(
        purity_resp_plot,
        tumor.purity.by.response(comb_clin_bb, response_cols),
        packages="ggplot2",
        format = "file"
    ),
    
    ##Figure S3b
    tar_target(
        purity_mod_plot,
        tumor.purity.tissue.by.mod.plot(comb_mes, comb_clin_bb, amtec_annots, response_cols),
        packages=c("ggplot2"),
        format = "file"
    ),
    
    ##Figure S3c
    tar_target(
        tissue_mod_plot,
        tissue.by.mod.plot(comb_mes, comb_clin_bb, amtec_annots, response_cols),
        packages=c("ggplot2"),
        format = "file"
    ),
    
    ##Figure S4
    tar_target(
        me_bb_plot,
        me.bb.heatmap(comb_mes, comb_clin_bb, amtec_annots, blia_blis_cols, response_cols),
        packages=c("ggplot2", "patchwork")
    ),
    
    ##Figure 3a
    tar_target(
        me_resp_fig_bx1,
        plot.modules.by.response(amtec_tcga_annot,comb_mes, response_cols, me_preds, "Bx1"),
        format = "file",
        packages="patchwork"
    ),
    
    ##Figure 3b
    tar_target(
        mzb1_eg_rel_plot,
        m3s.exprs.plot(comb_mes, comb_ltpm, amtec_tcga_annot, scaled_exprs, mm, mtnbc3s, mtnbc3s_top3),
        packages=c("ComplexHeatmap", "viridis"),
        format="file"
    ),
    
    ##Evaluation of MEs and other signatures
    
    tar_target(
        valid_me_preds,
        classify.valid.mes(valid_mes,  me_preds, valid_samps),
        packages=c("caret", "partykit", "yardstick")
    ),
    
    tar_target(
        mzb1_preds,
        ctree.mzb1(mm, amtec_kallisto, amtec_tcga_annot),
        packages=c("caret", "partykit", "yardstick")
    ),
    
    tar_target(
        valid_mzb1_preds,
        classify.valid.mzb1(mm, valid_samps, "data/amtec/validation_kallisto_counts.RData", mzb1_preds),
        packages=c("caret", "partykit", "yardstick")
    ),
    
    tar_target(
        ssgsea_preds,
        ctree.ssgsea(mm, amtec_kallisto, amtec_tcga_annot, mtnbc3s),
        packages=c("caret", "partykit", "GSVA", "yardstick")
    ),
    
    tar_target(
        valid_ssgsea_preds,
        classify.valid.ssgsea(mm, valid_samps, "data/amtec/validation_kallisto_counts.RData", ssgsea_preds, mtnbc3s),
        packages=c("caret", "partykit", "GSVA", "yardstick")
    ),
    
    tar_target(
        ssgsea_preds_top3,
        ctree.ssgsea(mm, amtec_kallisto, amtec_tcga_annot, mtnbc3s_top3),
        packages=c("caret", "partykit", "GSVA", "yardstick")
    ),
    
    tar_target(
        valid_ssgsea_top3_preds,
        classify.valid.ssgsea(mm, valid_samps, "data/amtec/validation_kallisto_counts.RData", ssgsea_preds_top3, mtnbc3s_top3),
        packages=c("caret", "partykit", "GSVA", "yardstick")
    ),
    
    tar_target(
        zhang_ctree_preds,
        classify.zhang.ctree(zhang_2021_sets, gene_info, amtec_tcga_annot, amtec_kallisto, valid_samps, "data/amtec/validation_kallisto_counts.RData"),
        packages=c("caret", "partykit", "yardstick")
    ),
    
    tar_target(
        comb_clin_ssg,
        merge.clin.ssg.preds(ssgsea_preds, valid_ssgsea_preds, comb_clin_bb, amtec_annots, valid_samps)
    ),
    
    tar_target(
        model_summary,
        summarize.models(comb_clin_ssg[(is.na(Consensus_Call)) | (Consensus_Call != "LAR")], me_preds, valid_me_preds, mzb1_preds, valid_mzb1_preds, 
                         ssgsea_preds, valid_ssgsea_preds, ssgsea_preds_top3, valid_ssgsea_top3_preds, zhang_ctree_preds),
        packages="stringr"
    ),
    
    tar_target(
        summary_matrix_list,
        generate.summary.matrix(comb_clin_ssg[(is.na(Consensus_Call)) | (Consensus_Call != "LAR")], me_preds, mzb1_preds, 
                                ssgsea_preds, ssgsea_preds_top3, zhang_ctree_preds, model_summary)
    ),
    
    ##Figure 3c
    tar_target(
        summary_mat_plot,
        plot.summary.matrix(summary_matrix_list, response_cols),
        packages=c("ggplot2"),
        format = "file"
        
    ),
    
    ##Validation
    ###Note this only shows METABRIC validation as CALGB is protected (No Figure 4a)
    tar_target(
        metabric_exprs,
        get.metabric.exprs("data/metabric/data_mrna_agilent_microarray.txt")
    ),
    
    tar_target(
        metabric_clin,
        get.metabric.clin("data/metabric/data_clinical_patient.txt", "data/metabric/data_clinical_sample.txt", metabric_exprs)
    ),
    
    tar_target(
        metabric_basal_exprs,
        subset.exprs(metabric_exprs, metabric_clin)
    ),
    
    tar_target(
        metabric_m3s,
        score.m3s(mm[,.(gene_id=gene_name)], mtnbc3s_top3[,.(gene_id=gene_name)], metabric_basal_exprs, "mTNBC3s_Top3"),
        packages=c("GSVA")
    ),
    
    tar_target(
        metabric_mzb1,
        compute.zscore(data.table(gene_id="MZB1", signature="MZB1"), metabric_basal_exprs)
    ),
    
    tar_target(
        combined_metabric_preds,
        combine.metabric.preds(metabric_clin, metabric_m3s, metabric_mzb1, 
                               ssgsea_preds_top3, mzb1_preds),
        packages="partykit"
    ),
    
    tar_target(
        metabric_full_models,
        fit.metabric.full.models(combined_metabric_preds),
        packages=c("survival", "data.table")
    ),
    
    tar_target(
        univ_meta_results,
        summarize.metabric.univ(combined_metabric_preds),
        packages=c("survival", "broom")
    ),
    
    ##Figure 4b
    tar_target(
        metabric_km_plot,
        plot.metabric.km(combined_metabric_preds),
        packages=c("survival", "ggsurvfit", "viridis"),
        format="file"
    ),
    
    ##Table S2b/c
    tar_target(
        written_full_models,
        write.full.models(metabric_full_models, univ_meta_results),
        packages="openxlsx",
        format="file"
    ),
    
    ##Differential expression vs network concepts
    
    tar_target(
        amtec_de,
        amtec.exprs.de.by.biopsy(amtec_kallisto, amtec_annots, gene_info),
        packages=c("edgeR", "limma")
    ),
    
    ##Table S3a
    tar_target(
        written_de_by_biopsy,
        write.de(amtec_de, "de_by_biopsy.xlsx"),
        package=c("openxlsx", "data.table"),
        format="file"
    ),
    
    tar_target(
        amtec_dnetc,
        amtec.dnet.by.biopsy(lion_concepts, amtec_tcga_annot, mm),
        packages=c("limma")
    ),
    
    ##Table S3b
    tar_target(
        written_dnetc_by_biopsy,
        write.de(amtec_dnetc, "dnetc_by_biopsy.xlsx"),
        package=c("openxlsx", "data.table"),
        format="file"
    ),
    
    tar_target(
        amtec_delta_de,
        amtec.exprs.de.delta(amtec_kallisto, amtec_annots, gene_info),
        packages=c("edgeR", "limma")
    ),
    
    ##Table S4a
    tar_target(
        written_delta_biopsy,
        write.de(amtec_delta_de, "delta_biopsy.xlsx"),
        package=c("openxlsx", "data.table"),
        format="file"
    ),
    
    tar_target(
        amtec_delta_dnetc,
        amtec.dnet.delta(amtec_tcga_annot, lion_concepts, mm),
        packages=c("limma")
    ),
    
    ##Table S4b
    tar_target(
        written_delta_dnetc,
        write.de(amtec_delta_dnetc, "delta_dnetc.xlsx"),
        package=c("openxlsx", "data.table"),
        format="file"
    ),
    
    ##Figure S5a
    tar_target(
        de_net_conc_plot_biopsy,
        de.by.net.conc.biopsy(amtec_de, amtec_dnetc),
        packages=c("ggplot2", "data.table"),
        format = "file"
    ),
    
    ##Figure S5b
    
    tar_target(
        de_net_conc_plot_delta,
        de.by.net.conc.delta(amtec_delta_de, amtec_delta_dnetc),
        packages=c("ggplot2", "data.table", "ggrepel")
    ),
    
    ##Figure 5b
    tar_target(
        ifi27_de_net_conc_plot_biopsy,
        ifi27.de.by.net.conc.biopsy(amtec_de, amtec_dnetc),
        packages=c("ggplot2", "data.table", "ggrepel"),
        format = "file"
    ),
    
    ##Figure 5c
    tar_target(
        lioness_heatmaps,
        make.lioness.mats(mm, lion_concepts, lion_network, amtec_tcga_annot),
        packages=c("ComplexHeatmap"),
        format="file"
    ),
    
    tar_target(
        gene_net_exprs_preds,
        ctree.gene.net.exprs(amtec_tcga_annot, lion_concepts, scaled_exprs),
        packages=c("caret", "partykit")
    ),
    
    tar_target(
        best_gene_net,
        data.table(desc="Bx1.Mod2 (blue).ENSG00000165949.12.Connectivity", cur_label="Mod2 (blue)", biopsy="Bx1", net_conc="Connectivity", gene_id="ENSG00000165949.12")
    ),
    
    ##Figure 5d
    tar_target(
        gene_net_perf_plot,
        plot.gene.net.results(mm, best_gene_net,amtec_tcga_annot, valid_tcga_annot, lion_concepts, valid_lion_concepts, gene_net_exprs_preds, scaled_exprs, response_cols),
        format="file",
        packages=c("ggplot2", "ggrepel")
    ),
    
    ##Figure 6a
    tar_target(
        krt23_de_net_conc_plot_biopsy,
        krt23.de.by.net.conc.delta(amtec_delta_de, amtec_delta_dnetc),
        packages=c("ggplot2", "data.table", "ggrepel")
    ),
    
    ##Figure 6b
    tar_target(
        krt23_bx_plot,
        plot.krt23.biopsy.mar(gene_net_exprs_preds, response_cols),
        format="file",
        packages=c("partykit","ggplot2")
    ),
    
    tar_target(
        delta_gene_net_exprs_preds,
        classify.gene.net.delta.ctree(mm, amtec_tcga_annot, lion_concepts, scaled_exprs, valid_tcga_annot, valid_lion_concepts),
        packages=c("caret", "partykit", "yardstick")
    ),
    
    ##Figure 6c
    tar_target(
        krt23_delta_plot,
        plot.krt23.delta.mar(amtec_tcga_annot, valid_tcga_annot, delta_gene_net_exprs_preds, response_cols),
        format="file",
        packages=c("partykit","ggplot2")
    )
    
    
    
    
)