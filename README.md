# Pandora Paper
This repository contains the code required to reproduce all datasets we used in the publication of our novel [Pandora tool](https://github.com/tschuelia/Pandora.git) for estimating the uncertainty of dimensionality reduction on population genetics data.

To execute the code, you need the Pandora installed on your system. See the [Pandora documentation](https://pandorageno.readthedocs.io/) for installation instructions.
We highly recommend the installation using conda.
You further need to install jupyter notebooks. If you installed Pandora using conda you can simply run 
```
conda install notebook -c conda-forge
```

## Reproducing datasets
1. Start the jupyter server by running `jupyter notebook` in your terminal. This should automatically start a new tab in your browser.
2. Open the jupyter file for the dataset you want to recreate. Note that for the Mathieson and Çayönü datasets, you first need to run the `NearEastPublic.ipynb` dataset as this provides the required modern samples for Mathieson and Çayönü.
3. Execute the jupyter notebook.
4. Done 🙂

You can now use the dataset to run the benchmark pipeline as described in the next section.

## Running the benchmark pipeline
In this repo we further provide the Snakemake pipeline we used to generate all Pandora results in the paper.

1. To run the pipeline, you need to install snakemake and the EIGENSOFT software package:
    ```
    conda install snakemake eigensoft -c bioconda
    ```
   
    Note: if you are working on a MacBook with an M1 or M2 chip, installing EIGENSOFT from conda won't work. 
    In this case, only install snakemake using conda and follow the instructions in the [Pandora wiki](https://pandorageno.readthedocs.io/en/latest/install.html#installing-eigensoft-on-macbooks-with-m1-m2-chips) to install EIGENSOFT.

2. Next, you need to configure the snakemake pipeline according to the dataset you want to analyze. You can do so by modifying the `config.yaml` file. See the provided `config.yaml` for hints on the allowed settings.
3. Run the pipeline in dry-run mode to see what rules will be run and what files created: `snakemake -n`
4. Finally, run the pipeline using `snakemake -c1`. Note that using `c1` ensures that all rules are executed sequentially. This is required to make sure that all benchmarks are correct and don't influence each other.
5. You can load the resulting `results.parquet` file using pandas:
   ```python
   import pandas as pd
   df = pd.read_parquet("results.parquet")
   ```



## Results and Pandora Reference
TODO: add link to preprint once available
TODO: add link to results of the benchmark runs once available
