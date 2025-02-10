#!/bin/bash
sbatch cells/cancer/cancer_cells_scGen_HDG.sh
sbatch cells/cancer/cancer_cells_scGen.sh
sbatch cells/cancer/cancer_cells_scVI_scANVI_2500.sh
sbatch cells/cancer/cancer_cells_scVI_scANVI_HDG.sh
