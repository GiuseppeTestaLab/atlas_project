# Running Liana on cancer and immune cells divided by treatment type

#%%
import scanpy as sc
# Only needed for visualization:
from scanpy.pl import umap
import liana as li
import pandas as pd
import plotnine as p9

## Inizialize directories
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'

#%%
adata = sc.read('atlas_embeddings_cell_labelled.h5ad')
adata
adata.obs
adata = adata[(adata.obs['max'] == 'CancerMSK') | (adata.obs['max'] == 'HematopoieticMSK')]

adata.obs = adata.obs.drop(columns=['cell_type', 'cell_subtype'])

#%%
adata_immune = adata[(adata.obs['max'] == 'HematopoieticMSK')]

#%%
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

#%%
cells = pd.DataFrame(adata_immune.obs['cell_types'], index=adata_immune.obs.index)
cells

adata.obs['cell_type'] = ''

for i in adata.obs.index:
    if i in cells.index:
        adata.obs['cell_type'][i] = cells['cell_types'][i]
    else:
        adata.obs['cell_type'][i] = 'CancerMSK'
        
#%%
values = []

for index, row in adata.obs.iterrows():
    if row['cell_type'] == 'CancerMSK' and row['tissue'] == 'Primary':
        values.append('Cancer_primary')
    elif row['cell_type'] == 'CancerMSK' and row['tissue'] == 'Ascites':
        values.append('Cancer_ascites')
    elif row['cell_type'] == 'CancerMSK' and row['tissue'] == 'Metastasis':
        values.append('Cancer_metastasis')
    elif row['cell_type'] == 'T_CD4_naive':
        values.append('T_CD4_cells')
    elif row['cell_type'] == 'T_CD4_CXCL13':
        values.append('T_CD4_cells')
    elif row['cell_type'] == 'T_CD4_reg':
        values.append('T_CD4_cells')
    elif row['cell_type'] == 'T_CD8_cytotoxic':
        values.append('T_CD8_cells')
    elif row['cell_type'] == 'T_CD8_CXCL13':
        values.append('T_CD8_cells')
    elif row['cell_type'] == 'NK_CD56':
        values.append('NK_cells')
    elif row['cell_type'] == 'NK_cytotoxic':
        values.append('NK_cells')
    else:
        values.append(row['cell_type'])

#%%
adata.obs['major_celltypes'] = values

#%%
adata_cht = adata[(adata.obs['treatment'] == 'CHT')]
adata_cht.obs['cell_type'] = adata_cht.obs['cell_type'].astype('category')
adata_naive = adata[(adata.obs['treatment'] == 'Naive')]
adata_naive.obs['cell_type'] = adata_naive.obs['cell_type'].astype('category')
adata_nact = adata[(adata.obs['treatment'] == 'NACT')]
adata_nact.obs['cell_type'] = adata_nact.obs['cell_type'].astype('category')

adata_cht = adata[(adata.obs['treatment'] == 'CHT')]
adata_cht.obs['major_celltypes'] = adata_cht.obs['major_celltypes'].astype('category')
adata_naive = adata[(adata.obs['treatment'] == 'Naive')]
adata_naive.obs['major_celltypes'] = adata_naive.obs['major_celltypes'].astype('category')
adata_nact = adata[(adata.obs['treatment'] == 'NACT')]
adata_nact.obs['major_celltypes'] = adata_nact.obs['major_celltypes'].astype('category')

#%%
## Methods
li.mt.show_methods()

#%%
## import liana's rank_aggregate
from liana.mt import rank_aggregate
## resources
li.resource.show_resources()

#%%
## import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

#%%
# Run rank_aggregate
li.mt.rank_aggregate(adata_cht, groupby='major_celltypes', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_naive, groupby='major_celltypes', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_nact, groupby='major_celltypes', expr_prop=0.1, verbose=True, use_raw=False)

#%%
adata_cht.uns['liana_res'].head()
adata_naive.uns['liana_res'].head()
adata_nact.uns['liana_res'].head()

rank_aggregate.describe()

#%%
my_plot = li.pl.dotplot(adata = adata_cht,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
                             'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
              target_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
                             'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(27, 10)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11)
 )
)

plot
plot.save('Figures/dotplot_rankaggr_cancer_subtypes_immune_subtypes_cht.png', limitsize=False)

#%%
my_plot = li.pl.dotplot(adata = adata_naive,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
                             'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
              target_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
                             'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(27, 10)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11)
 )
)

plot
plot.save('Figures/dotplot_rankaggr_cancer_subtypes_immune_subtypes_naive.png', limitsize=False)

#%%
my_plot = li.pl.dotplot(adata = adata_nact,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
                             'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
              target_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
                             'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(27, 10)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11)
 )
)

plot
plot.save('Figures/dotplot_rankaggr_cancer_subtypes_immune_subtypes_nact.png', limitsize=False)

#%%
adata_cht.write_h5ad('/group/testa/Project/OvarianAtlas/liana_results/liana_aggregaterank_cancer_subtypes_immune_subtypes_cht.h5ad')
adata_naive.write_h5ad('/group/testa/Project/OvarianAtlas/liana_results/liana_aggregaterank_cancer_subtypes_immune_subtypes_naive.h5ad')
adata_nact.write_h5ad('/group/testa/Project/OvarianAtlas/liana_results/liana_aggregaterank_cancer_subtypes_immune_subtypes_nact.h5ad')
