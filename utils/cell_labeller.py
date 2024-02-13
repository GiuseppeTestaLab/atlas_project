
import scanpy as sc
import numpy as np

def assign_scores(adata):
    
    sc.tl.score_genes(adata, 
    ['WFDC2','CD24','CLDN3','KRT7','KRT8','KRT17','KRT18','KRT19','EPCAM','WT1','CLDN4','MSLN','FOLR1','MUC1'], 
    score_name = "CancerMSK", use_raw=True)

    sc.tl.score_genes(adata, 
    ['CLEC14A','PECAM1','VWF','CAV1','EMCN','CDH5','MCAM','IL3RA','IGFBP7','COL4A1','COL4A2','COL15A1','SPARCL1','A2M','HSPG2','PLVAP','AQP1','ENG','RAMP2','GNG11','EGFL7','CLDN5','INSR'], 
    score_name = "EndothelialMSK", use_raw=True)

    sc.tl.score_genes(adata, 
    ['COL1A1','COL3A1','WT1','ACTA2','CAV1','COL1A2','DCN','SPARC','COL6A1','CCDC80','LUM','COL6A2','COL6A3','CALD1','RARRES2','MGP','CTHRC1','AEBP1','POSTN','COL5A2','FBLN1','TAGLN','C1S','C1R','NNMT','MMP2','IGFBP5','TIMP1','FN1','IGFBP7','C3','COL5A1','LGALS1'], 
    score_name = "FibroblastsMSK", use_raw=True)

    sc.tl.score_genes(adata, 
    ['PTPRC'], 
    score_name = "HematopoieticMSK", use_raw=True)

    return adata

# Cell labelling strategy

def actual_labeller(adata):
    def positive_ratio(x):
        x[x<0]=0
        if (sum(x)-max(x))==0:
            return max(x)
        else:
            return max(x)/(sum(x)-max(x))

    adata.obs['cell_labels_ratio']= adata.obs[['EndothelialMSK', 'FibroblastsMSK', 'CancerMSK', 'HematopoieticMSK']].apply(lambda x:positive_ratio(x), axis=1)

    np.log(adata.obs.loc[adata.obs.cell_labels_ratio > 0.1, 'cell_labels_ratio']).hist(bins=30)

    adata.obs['max'] = adata.obs[['EndothelialMSK', 'FibroblastsMSK', 'CancerMSK','HematopoieticMSK']].idxmax(axis=1)

    percentage = 0.2
    counter = 0
    for dataset in adata.obs.dataset.unique():
        print("Processing the dataset: {}".format(dataset))
        adata_slice = sc.AnnData(X=adata[adata.obs.dataset == dataset].X, 
                                var=adata[adata.obs.dataset == dataset].var, 
                                obs = adata[adata.obs.dataset == dataset].obs) 
        sc.tl.pca(adata_slice, svd_solver='arpack')
        sc.pp.neighbors(adata_slice, n_neighbors=15, n_pcs=50)
        sc.tl.leiden(adata_slice,resolution = 1.8, key_added="leiden-{}".format('1.8'))
        sc.tl.umap(adata_slice)
        for leiden in adata_slice.obs["leiden-1.8"].unique():
            n_cells_xpopul = adata_slice.obs.loc[adata_slice.obs['leiden-1.8'] == leiden, 'max'].value_counts()[0]
            n_cells_total = adata_slice.obs.loc[adata_slice.obs['leiden-1.8'] == leiden, 'max'].shape[0]
            assignment = adata_slice.obs.loc[adata_slice.obs['leiden-1.8'] == leiden, 'max'].value_counts().index[0]
            if (n_cells_xpopul / n_cells_total) > percentage:
                adata_slice.obs["assignment"] = assignment
        #putting the assignment in the initial adata
        if counter == 0:
            adata.obs['assignment'] = adata_slice.obs['assignment']
            adata.obs['leiden-1.8'] = adata_slice.obs['leiden-1.8']
        else :
            adata.obs.loc[adata_slice.obs.index, 'assignment'] = adata_slice.obs['assignment']
            adata.obs['leiden-1.8'] = adata.obs['leiden-1.8'].astype('str')
            adata.obs.loc[adata_slice.obs.index, 'leiden-1.8'] = adata_slice.obs['leiden-1.8'].astype('str')
        counter = counter + 1

    adata = adata.raw.to_adata()

    return adata

def create_cancer_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'CancerMSK'].X, 
                            var=adata[adata.obs['max'] == "CancerMSK"].var, 
                            obs = adata[adata.obs['max'] == 'CancerMSK'].obs)
    
    return adata

def create_immune_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'HematopoieticMSK'].X, 
                            var=adata[adata.obs['max'] == "HematopoieticMSK"].var, 
                            obs = adata[adata.obs['max'] == 'HematopoieticMSK'].obs)
    
    return adata

def create_fibroblast_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)
    
    return adata

def create_endothelial_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'EndothelialMSK'].X, 
                            var=adata[adata.obs['max'] == "EndothelialMSK"].var, 
                            obs = adata[adata.obs['max'] == 'EndothelialMSK'].obs)
    
    return adata

# Immune cell labeller according to the Vasquez-Garcia et al. signatures genes

def immune_labeller(adata):
    adata_immune = adata[(adata.obs['max'] == 'HematopoieticMSK')]

    # Plasma cells markers
    sc.tl.score_genes(adata_immune, ['IGKC','IGHG1','CD79A','IGHG2','IGLC2','IGLC3','IGHG3','IGHG4','JCHAIN','MZB1','XBP1'], 
    score_name = "Plasma_cells", use_raw=False)

    # Mast cells markers
    sc.tl.score_genes(adata_immune, ['KIT','CPA3','CTSG','MS4A2','TPSAB1','TPSB2','HPGD','HPGDS','GATA2'], 
    score_name = "Mast_cells", use_raw=False)

    # B cells markers
    sc.tl.score_genes(adata_immune, ['MS4A1', 'CD79A', 'CD19', 'BANK1', 'IGKC', 'IGHM'], score_name = "B_cells", use_raw=False)

    # Macrophages M1 cells markers
    sc.tl.score_genes(adata_immune, ['S100A9','S100A8','VCAN','FCN1','S100A12','THBS1','LYZ','CD55','RETN','CD52','EREG','AC020656.1','CD300E','APOBEC3A','AREG','CFP','SMIM25',
                            'MCEMP1','S100A4','RIPOR2','CD36','S100A6','CYP1B1','SERPINA1','STXBP2','LST1','LILRA5','SLC25A37','FPR1','MNDA','LILRB2','TIMP1','FGR',
                            'NAMPT','COTL1','IRAK3','LYST','SH3BGRL3','NCF2','CCR2','FLNA','MXD1','LTA4H','RGS2','TUBA1A','H3F3A','CEBPB','CSF3R','CSTA','LCP1'], 
    score_name = "M1_macrophages", use_raw=False)

    # Myeloid cells markers
    sc.tl.score_genes(adata_immune, 
    ['CD14','FCER1G','FCGR3A','LYZ','CTSS','CD33','CD68','CD163','ITGAX','ITGAM','CD4','MRC1',
    'VSIG4','SPP1','APOE','C1QA','C1QB','C1QC','APOC1','FTL','S100A9','TYROBP','AIF1','CD74','PSAP','CTSB'], 
    score_name = "Myeloid_cells", use_raw=False)

    # Dendritic cells markers
    sc.tl.score_genes(adata_immune, 
    ['IL3RA','IRF7','IRF8','GZMB','CD4','CLEC4C','JCHAIN',
    'PTGDS','PLAC8','PLD4','TCF4','BCL11A','GPR183','CCDC50','LILRA4','TSPAN13','CLIC3','MPEG1'], 
    score_name = "Dendritic_cells", use_raw=False)

    # T CD4 naive cells markers
    sc.tl.score_genes(adata_immune, ['IL7R','CCR7','KLF2','EEF1B2','TPT1','EEF1A1','TCF7','MAL','CD40LG','GPR183','LDHB','SELL','SNHG8','NOSIP','PABPC1','NOP53','LEF1',
                            'LTB','EIF3E','RACK1','JUNB','NACA','SOCS3','TOMM7','UBA52','TMEM123','SERINC5','EEF2','FXYD5','TRABD2A','TSHZ2','SARAF','AQP3',
                            'ANK3','RIPOR2','AP3M2','TOB1','ZFAS1','LINC02273','EIF4B','ANXA1','NSA2','TNFRSF25','CTSL','SESN3','EEF1D','FAU','LDLRAP1','FLT3LG','TIMP1'], 
    score_name = "T_CD4_naive", use_raw=False)

    # T CD4 CXCL13 markers
    sc.tl.score_genes(adata_immune, ['CXCL13','NMB','NR3C1','FKBP5','IL6ST','MAF','ITM2A','CTLA4','TSHZ2','LIMS1','CD40LG','PDCD1','TNFRSF4','CD4','RNF19A','RBPJ','CORO1B','CPM',
                            'ZBED2','AHI1','ICA1','TOX2','DUSP4','AC004585.1','ARID5B','CCDC50','CD84','IGFL2','SRGN','BATF','CH25H','TNFRSF18','SPOCK2','CHN1','CD200',
                            'RGS1','RILPL2','ZNRF1','TNFRSF25','METTL8','SLA','SMCO4','BTLA','SESN3','NAP1L4','BHLHE40-AS1','MIR155HG','BIRC3','PTPN13','CYSLTR1'], 
    score_name = "T_CD4_CXCL13", use_raw=False)

    # T CD4 reg markers
    sc.tl.score_genes(adata_immune, ['TNFRSF4','IL2RA','FOXP3','CTLA4','LTB','RTKN2','BATF','TNFRSF18','SAT1','TBC1D4','TIGIT','GADD45A','TNFRSF1B','PMAIP1','UGP2','IKZF2','TNFRSF9',
                            'ICOS','SOX4','LINC01943','IL32','ARID5B','LAYN','CD27','BIRC3','CORO1B','TYMP','CD4','DUSP4','ENTPD1','CTSC','MIR4435-2HG','LINC02099','MAGEH1',
                            'SPOCK2','CARD16','PHACTR2','S100A4','STAM','SPATS2L','GLRX','AC005224.3','MAF','BTG3','PBXIP1','F5','SLAMF1','IL1R1','DNPH1','TRAC'], 
    score_name = "T_CD4_reg", use_raw=False)

    # T CD8 cytotoxic markers
    sc.tl.score_genes(adata_immune, ['GZMK','CD8A','CD8B','ITM2C','GZMH','CCL5','TRGC2','GZMA','KLRG1','CCL4','CRTAM','CST7','GZMM','DTHD1','HLA-DPB1','PPP1R14B','CD3G','THEMIS','EOMES',
                            'TC2N','DUSP2','LYAR','CD3D','PPP2R5C','SLF1','KIAA1551','CXCR6','YBX3','HLA-DPA1','CCL4L2','F2R','CXCR4','FAM102A','HLA-DRB1','SLAMF7','APOBEC3G',
                            'SH2D1A','CD84','STK17A','CCR5','TUBA4A','ARAP2','GPR174','PECAM1'], 
    score_name = "T_CD8_cytotoxic", use_raw=False)

    # T CD8 CXCL13 markers
    sc.tl.score_genes(adata_immune, ['CXCL13','GZMB','CCL4L2','MIR155HG','TNFRSF9','HAVCR2','RBPJ','LAG3','IFNG','PTMS','CCL3','CD8A','CRTAM','FABP5','PHLDA1','JAML','TIGIT','KRT86',
                            'CCL5','CXCR6','LINC01871','PDCD1','HLA-DRB1','TNIP3','GAPDH','CD63','FAM3C','GZMH','CTLA4','CCND2','SPRY1','CD8B','VCAM1','HLA-DRA','ID2','ITGAE',
                            'DUSP4','LYST','ENTPD1','SRGAP3','TNFSF4','NDFIP2','GOLIM4','AKAP5','CD27','HLA-DPA1','SNAP47','RGS1','ITM2A','APOBEC3C'], 
    score_name = "T_CD8_CXCL13", use_raw=False)

    # T CD8 ISG markers
    sc.tl.score_genes(adata_immune, ['IFIT3','ISG15','IFIT1','MX1','IFIT2','RSAD2','IFI6','MX2','ISG20','IFI44L','HERC5','OAS1','SAMD9L','TNFSF10','STAT1','EIF2AK2','GBP1','MT2A','OAS3',
                            'EPSTI1','SAMD9','PLSCR1','IFI35','XAF1','OASL','IFI44','USP18','LY6E','CMPK2','NT5C3A','DDX58','HELZ2','IRF7','TRIM22','PARP14','DDX60','LAG3',
                            'DDX60L','IFIH1','PPM1K','OAS2','RNF213','PARP9','PNPT1','SP110','SAT1','C19orf66','STAT2','BST2','LAP3'], 
    score_name = "T_CD8_ISG", use_raw=False)

    # Innate lymphoid cells markers
    sc.tl.score_genes(adata_immune, ['KLRB1','IL7R','IL4I1','CEBPD','LST1','TNFSF13B','LTB','FOS','SLC4A10''CCR6','NFKBIA','RORA','NCR3','TNFAIP3','AQP3','MYBL1','DUSP1','LINC01871','CCL20',
                            'JAML','CTSH','TMIGD2','TNFRSF25','ERN1','DPP4','KLRG1','SPOCK2','ZBTB16','IFNGR1','FKBP11','TPT1','MGAT4A','PDCD4','S100A4','SATB1','S100A6','CD40LG',
                            'B3GALT2','ABCB1','RUNX2','TLE1','EEF1A1','CERK','RORC','PERP','LTK','PLCB1','LTC4S','EEF1B2','KIF5C'], 
    score_name = "ILC", use_raw=False)

    # NK CD56 cells markers
    sc.tl.score_genes(adata_immune, ['GNLY','TYROBP','AREG','KLRC1','FCER1G','TRDC','KRT81','XCL1','KLRD1','IGFBP2','XCL2','CLIC3','KRT86','IL2RB','CEBPD','CTSW','TXK','MATK','KLRB1',
                            'CD7','CD63','NKG7','CCL3','TMIGD2','HOPX','TNFRSF18','CMC1','GSTP1','SRGAP3','KLRC2','LAT2','GZMB','LINC00996','NCAM1','PRF1','CXXC5','IFITM3',
                            'ZNF683','KLRF1','MCTP2','SH2D1B','ITGA1','IFITM2','CCL5','CD38','SLC16A3','ITGAX','CAPN12','CD247','SAMD3'], 
    score_name = "NK_CD56", use_raw=False)

    # NK cytotoxic cells markers
    sc.tl.score_genes(adata_immune, ['FGFBP2','FCGR3A','SPON2','PRF1','KLRF1','GNLY','KLRD1','NKG7','CX3CR1','GZMB','PLAC8','CLIC3','PLEK','TYROBP','PTGDS','EFHD2','FCER1G','CST7','GZMH',
                            'ADGRG1','CCL3','HOPX','ZEB2','IGFBP7','PRSS23','CD247','AKR1C3','C1orf21','MYBL1','AREG','S1PR5','CTSW','KLF2','ABHD17A','TTC38','PTPN12','KLRB1','CCL4',
                            'PTGDR','TRDC','ITGB2','XBP1','CEP78','CMC1','LITAF','BIN2','CHST2','CD300A','ARL4C','TXK'], 
    score_name = "NK_cytotoxic", use_raw=False)

    adata_immune.obs['cell_types'] = adata_immune.obs[['Plasma_cells','Mast_cells', 'B_cells', 'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells',
                                                        'T_CD4_naive', 'T_CD4_CXCL13', 'T_CD4_reg', 'T_CD8_cytotoxic', 'T_CD8_CXCL13', 'T_CD8_ISG',
                                                        'ILC', 'NK_CD56', 'NK_cytotoxic']].idxmax(axis=1)
    
    return adata_immune

# Fibroblast cell labeller according to literature derived gene signatures

def fibroblast_labeller(adata):
    adata_fibroblasts = adata[(adata.obs['max'] == 'FibroblastsMSK')]
    
    # inflammatory CAFs (iCAF)
    sc.tl.score_genes(adata_fibroblasts, ['IL6', 'CXCL12', 'VTN', 'CADM3', 'PLIN2', 'SERPINB2', 'KRT19', 'DES', 'CALB2', 'WT1', 'KRT7'], 
    score_name = "iCAF", use_raw=False)

    # matrix CAFs (mCAF)
    sc.tl.score_genes(adata_fibroblasts, ['ACTA2', 'COL11A1', 'MFAP5', 'SFRP2', 'ISLR', 'COL10A1', 'VIM', 'COL3A', 
                        'MMP11', 'PTHLH', 'FGF1', 'WNT7B', 'WNT2', 'TGFB3', 'THRC1', 'POSTN', 'VCAN', 'ZEB1'], 
    score_name = "mCAF", use_raw=False)

    # vascular cancer-associated fibroblasts CAFs (vCAF)
    sc.tl.score_genes(adata_fibroblasts, ['COX4I2', 'HIGD1B', 'PTP4A3', 'MCAM', 'PPP1R14A'], 
    score_name = "vCAF", use_raw=False)

    # STAR gene expressing CAFs (starCAF)
    sc.tl.score_genes(adata_fibroblasts, ['STAR', 'IGFBP5', 'TSPAN8', 'C7', 'ALDH1A1', 'LGR5'], 
    score_name = "starCAF", use_raw=False)

    adata_fibroblasts.obs['cell_types'] = adata_fibroblasts.obs[['iCAF', 'mCAF','vCAF', 'starCAF']].idxmax(axis=1)

    return adata_fibroblasts
