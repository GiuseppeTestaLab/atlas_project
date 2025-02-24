#!/bin/bash
find /group/testa/Project/OvarianAtlasTestStep0/ -iname *HDG*.h5ad -type f -exec sbatch to_pandas.sh \{\} \;
find /group/testa/Project/OvarianAtlas/ -iname *HDG*.h5ad -type f -exec sbatch to_pandas.sh \{\} \;