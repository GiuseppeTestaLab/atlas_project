#!/bin/sh
INIT_DIR=/group/testa/Project/OvarianAtlas/atlas_project/raw_data
DEST_DIR=$1

ORIGINAL_COUNTS=$INIT_DIR/original_counts_backup/original_counts
DEST_COUNTS=$DEST_DIR/original_counts

mkdir -p $DEST_DIR
mkdir -p $DEST_COUNTS/Geistlinger2020/
cp -R $ORIGINAL_COUNTS/Geistlinger2020/Samples $DEST_COUNTS/Geistlinger2020/

mkdir -p $DEST_COUNTS/Loret2022/
cp -R $ORIGINAL_COUNTS/Loret2022/Counts $DEST_COUNTS/Loret2022/

mkdir -p $DEST_COUNTS/Olbrecht2021/
cp -R $ORIGINAL_COUNTS/Olbrecht2021/10xCounts $DEST_COUNTS/Olbrecht2021/

mkdir -p $DEST_COUNTS/Qian2020/
cp -R $ORIGINAL_COUNTS/Qian2020/Counts $DEST_COUNTS/Qian2020/

mkdir -p $DEST_COUNTS/Regner2021/
cp -R $ORIGINAL_COUNTS/Regner2021/HGSOC_tum $DEST_COUNTS/Regner2021/

mkdir -p $DEST_COUNTS/Ren2022/
cp -R $ORIGINAL_COUNTS/Ren2022/Counts $DEST_COUNTS/Ren2022/

mkdir -p $DEST_COUNTS/Vasquez2022/
cp -R $ORIGINAL_COUNTS/Vasquez2022/GSE180661_matrix.h5ad $DEST_COUNTS/Vasquez2022/

mkdir -p $DEST_COUNTS/Xu2022/
cp -R $ORIGINAL_COUNTS/Xu2022/Cancer $DEST_COUNTS/Xu2022/

mkdir -p $DEST_COUNTS/Zhang2022/
cp -R $ORIGINAL_COUNTS/Zhang2022/GSE165897* $DEST_COUNTS/Zhang2022/

mkdir -p $DEST_COUNTS/Zheng2023/
cp -R /group/testa/Project/OvarianAtlas/Zheng2023/Counts $DEST_COUNTS/Zheng2023/