#!/bin/sh
URL=https://repo.bioserver.ieo.it/GT/original_counts/
source bash_ini_parser/read_ini.sh
read_ini config.ini

DEST_RAW=$INI__DEFAULT__rawPath

DEST_COUNTS=$DEST_RAW/original_counts

echo "Builing directory tree for raw data at:" $DEST_RAW
cd $DEST_RAW
wget -mpEk --no-parent --cut-dirs=1 -nH --reject="index.html*" $URL --directory-prefix $DEST_RAW

mkdir -p $DEST_COUNTS/Geistlinger2020/Adata
mkdir -p $DEST_COUNTS/Loret2022/Adata
mkdir -p $DEST_COUNTS/Olbrecht2021/Adata
mkdir -p $DEST_COUNTS/Qian2020/Adata
mkdir -p $DEST_COUNTS/Regner2021/Adata
mkdir -p $DEST_COUNTS/Ren2022/Adata
mkdir -p $DEST_COUNTS/Vasquez2022/Adata
mkdir -p $DEST_COUNTS/Xu2022/Adata
mkdir -p $DEST_COUNTS/Zhang2022/Adata
mkdir -p $DEST_COUNTS/Zheng2023/Adata

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
mkdir -p $DEST_RAW/integration/cells/endothelial
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
mkdir -p $DEST_PLOT/integration/cells/endothelial
mkdir -p $DEST_PLOT/integration/cells/fibroblasts
mkdir -p $DEST_PLOT/integration/cells/immune
