# Pandora Paper

This repository contains the code required to reproduce all datasets we used in the publication of our
novel [Pandora tool](https://github.com/tschuelia/Pandora.git) for estimating the uncertainty of dimensionality
reduction on population genetics data.

## Setting up the environment

To run the code, you need a few packages installed. You can use the provided `environment.yaml` file and conda:

```
conda env create --file environment.yaml
```

Note: if you are working on a MacBook with an M1 or M2 chip, installing EIGENSOFT from conda won't work.
In this case, remove the `- eigensoft` line from the environment file before executing the `conda` command above.
Install EIGENSOFT by following the instructions in
the [Pandora documentation](https://pandorageno.readthedocs.io/en/latest/install.html#installing-eigensoft-on-macbooks-with-m1-m2-chips).

## Reproducing datasets

We implemented Snakemake pipelines for the empirical dataset analyses, as well as for the data simulation and subsequent
analysis of the simulated datasets used in our paper.

### Empirical datasets

1. Start the jupyter server by running `jupyter notebook` in your terminal. This should automatically start a new tab in
   your browser.
2. Open the jupyter file for the dataset you want to recreate. Note that for the Çayönü datasets, you
   first need to run the `NearEastPublic.ipynb` dataset as this provides the required modern samples the ancient
   Çayönü data.
3. Execute the jupyter notebook.
4. Done 🙂

You can now use the dataset to run the benchmark pipeline:

1. First, you need to configure the snakemake pipeline according to the dataset you want to analyze. You can do so by
   modifying the `Snakefile` file in the `empirical` directory. The config section is marked by a `# CONFIG` comment.
   Note that this pipeline first copies the empirical dataset to a local temporary directory. This was required for our
   lab server setup due to unreliable benchmark results caused by the NFS.
2. Once you setup the config according to your data, run the pipeline in dry-run mode to see what rules will be run and
   what files created: `snakemake -n`
3. Finally, run the pipeline using `snakemake -c1`. Note that using `c1` ensures that all rules are executed
   sequentially. This is required to make sure that benchmarks are correct and don't influence each other. Fun fact: you
   can't make use of more than one core in this analysis anyway as there are no rules that can be run in parallel :-)
4. The final output of the pipeline will be a single `.parquet` file containing a summary of the Pandora run and
   results.

### Simulated Datasets

You can repeat our simulation analysis using the provided Snakemake pipeline in the `popgen_simulations` directory.
This will generate all simulated datasets as described in our paper and will also run Pandora under the configured
settings per dataset.

1. Modify the `# CONFIG` section in the `Snakefile` according to your setup. For our paper, we executed the pipeline
   using the following combination of settings:
    1. `SEQ_LENGTHS = ["1e5", "1e6", "1e7", "1e8"]`, `MISSING = []`, `NOISE = []`
    2. `SEQ_LENGTHS = ["1e8"]`, `MISSING = [0.01, 0.05, 0.1, 0.2, 0.5]`, `NOISE = [0.1, 0.2, 0.5]`
2. Similar to the above instructions for the empirical pipeline, first make a dry-run to see what data will be generated
   and analyzed: `snakemake -n`
3. Run the pipeline using `snakemake -c1`. This will generate all datasets and run Pandora on them. Again, we recommend
   using `-c1` to ensure that all rules are executed sequentially for correct benchmark results.
4. The final output of the pipeline will be a single `.parquet` file containing a summary of the Pandora runs and
   results.


## Download the results

Instead of rerunning the analyses yourself, you can also download the simulated datasets and results
via [our lab server](https://cme.h-its.org/exelixis/material/Pandora_supplementary_data.tar.gz).
The directory follows this structure:

```text
├── empirical_datasets 
│   ├── HO-Cayonu
    │   ├── HO-Cayonu.pandora.parquet
    │   ├── bootstrap_{3,11}.{eval,evec}
    │   ├── HO-Cayonu.{eval,evec}
    │   └── modern.poplist.txt  # names of the HO-West modern populations
    └── DogGenomes
        ├── HighlyAdmixedDogs.pandora.parquet
        └── PureBredDogs.pandora.parquet
└── simulated_datasets 
    ├── {model} / {seq_len}
    ├── dataset.{ind,geno,snp}
    ├── {missing}  # only for seq_len=1e8, missing in [0.01, 0.05, 0.1, 0.2, 0.5]
    │   ├── missing_0.01.{ind,geno,snp}
    │   ├── missing_0.05.{ind,geno,snp}
    │   ├── missing_0.1.{ind,geno,snp}
    │   ├── missing_0.2.{ind,geno,snp}
    │   └── missing_0.5
    ├── {noise}  # only for seq_len=1e8, noise in [0.1, 0.2, 0.5]
    │   ├── noise_0.1.{ind,geno,snp}
    │   ├── noise_0.2.{ind,geno,snp}
    │   └── noise_0.5.{ind,geno,snp}
    └── results
        ├──{PCA/MDS}
           ├── no_convergence.parquet
           ├── convergence_5.parquet
           ├── convergence_1.parquet
           ├── PSVs_no_convergence.parquet
           ├── PSVs_convergence_5.parquet
           └── PSVs_convergence_1.parquet
```



## Publication

The paper explaining the details of Pandora is available as preprint on bioRxiv:

Haag, J., Jordan A. I. & Stamatakis, A. (2024). **Pandora: A Tool to Estimate Dimensionality Reduction Stability of
Genotype Data.** *bioRxiv*. [https://doi.org/10.1101/2024.03.14.584962](https://doi.org/10.1101/2024.03.14.584962)
