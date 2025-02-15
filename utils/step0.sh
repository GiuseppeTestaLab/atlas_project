#!/bin/sh
INIT_DIR=/group/testa/Project/OvarianAtlas/atlas_project/raw_data

source bash_ini_parser/read_ini.sh
read_ini config.ini

DEST_RAW=$INI__DEFAULT__rawPath

ORIGINAL_COUNTS=$INIT_DIR/original_counts_backup/original_counts
DEST_COUNTS=$DEST_RAW/original_counts

echo "Builing directory tree for raw data at:" $DEST_RAW

mkdir -p $DEST_COUNTS/Geistlinger2020/Adata
cp -R $ORIGINAL_COUNTS/Geistlinger2020/Samples $DEST_COUNTS/Geistlinger2020/

mkdir -p $DEST_COUNTS/Loret2022/Adata
cp -R $ORIGINAL_COUNTS/Loret2022/Counts $DEST_COUNTS/Loret2022/

mkdir -p $DEST_COUNTS/Olbrecht2021/Adata
cp -R $ORIGINAL_COUNTS/Olbrecht2021/10xCounts $DEST_COUNTS/Olbrecht2021/

mkdir -p $DEST_COUNTS/Qian2020/Adata
cp -R $ORIGINAL_COUNTS/Qian2020/Counts $DEST_COUNTS/Qian2020/

mkdir -p $DEST_COUNTS/Regner2021/Adata
cp -R $ORIGINAL_COUNTS/Regner2021/HGSOC_tum $DEST_COUNTS/Regner2021/

mkdir -p $DEST_COUNTS/Ren2022/Adata
cp -R $ORIGINAL_COUNTS/Ren2022/Counts $DEST_COUNTS/Ren2022/

mkdir -p $DEST_COUNTS/Vasquez2022/Adata
cp -R $ORIGINAL_COUNTS/Vasquez2022/GSE180661_matrix.h5ad $DEST_COUNTS/Vasquez2022/

mkdir -p $DEST_COUNTS/Xu2022/Adata
cp -R $ORIGINAL_COUNTS/Xu2022/Cancer $DEST_COUNTS/Xu2022/

mkdir -p $DEST_COUNTS/Zhang2022/Adata
cp -R $ORIGINAL_COUNTS/Zhang2022/GSE165897* $DEST_COUNTS/Zhang2022/

mkdir -p $DEST_COUNTS/Zheng2023/Adata
cp -R /group/testa/Project/OvarianAtlas/Zheng2023/Counts $DEST_COUNTS/Zheng2023/


mkdir -p $DEST_RAW/atlas_annotated
mkdir -p $DEST_RAW/metacells/endothelial/
mkdir -p $DEST_RAW/metacells/immune/
mkdir -p $DEST_RAW/metacells/fibroblasts
mkdir -p $DEST_RAW/metacells/cancer
mkdir -p $DEST_RAW/original_anndata/Zhang2022
mkdir -p $DEST_RAW/original_anndata/Ren2022
mkdir -p $DEST_RAW/original_anndata/Xu2022
mkdir -p $DEST_RAW/original_anndata/Qian2020
mkdir -p $DEST_RAW/original_anndata/Regner2021
mkdir -p $DEST_RAW/original_anndata/Geistlinger2020
mkdir -p $DEST_RAW/original_anndata/Olbrecht2021
mkdir -p $DEST_RAW/original_anndata/Loret2022
mkdir -p $DEST_RAW/original_anndata/Vasquez2022
mkdir -p $DEST_RAW/integration/metacells/cancer
mkdir -p $DEST_RAW/integration/metacells/endothelial
mkdir -p $DEST_RAW/integration/metacells/fibroblasts
mkdir -p $DEST_RAW/integration/metacells/immune
mkdir -p $DEST_RAW/integration/cells/cancer
mkdir -p $DEST_RAW/integration/cells/endhotelial
mkdir -p $DEST_RAW/integration/cells/fibroblasts
mkdir -p $DEST_RAW/integration/cells/immune
mkdir -p $DEST_RAW/metacells/cancer
mkdir -p $DEST_RAW/metacells/endothelial
mkdir -p $DEST_RAW/metacells/fibroblasts
mkdir -p $DEST_RAW/metacells/immune


DEST_PLOT=$INI__DEFAULT__figPath
echo "Builing directory tree for plots at:" $DEST_PLOT
mkdir -p $DEST_PLOT/atlas_annotated
mkdir -p $DEST_PLOT/integration/metacells/cancer
mkdir -p $DEST_PLOT/integration/metacells/endothelial
mkdir -p $DEST_PLOT/integration/metacells/fibroblasts
mkdir -p $DEST_PLOT/integration/metacells/immune
mkdir -p $DEST_PLOT/integration/cells/cancer
mkdir -p $DEST_PLOT/integration/cells/endhotelial
mkdir -p $DEST_PLOT/integration/cells/fibroblasts
mkdir -p $DEST_PLOT/integration/cells/immune
