# Pandora Paper
This repository contains the code required to reproduce all datasets we used in the publication of our novel [Pandora tool](https://github.com/tschuelia/Pandora.git) for estimating the uncertainty of dimensionality reduction on population genetics data.


## Setting up the environment
To run the code, you need a few packages installed. You can use the provided `environment.yaml` file and conda:
```
conda env create --file environment.yaml
```

Note: if you are working on a MacBook with an M1 or M2 chip, installing EIGENSOFT from conda won't work.
In this case, remove the `- eigensoft` line from the environment file before executing the `conda` command above. 
Install EIGENSOFT by following the instructions in the [Pandora documentation](https://pandorageno.readthedocs.io/en/latest/install.html#installing-eigensoft-on-macbooks-with-m1-m2-chips).


## Reproducing datasets
Since the empirical datasets are based on genotype files in EIGENSTRAT format, we implemented a Snakemake pipeline for benchmarking
making use of the Pandora command line interface.

The simulated datasets we used in our paper are `RData` exports which we have to load as Pandora `NumpyDatasets`, thus using Pandora as Python library.

### Empirical datasets
1. Start the jupyter server by running `jupyter notebook` in your terminal. This should automatically start a new tab in your browser.
2. Open the jupyter file for the dataset you want to recreate. Note that for the Mathieson and Çayönü datasets, you first need to run the `NearEastPublic.ipynb` dataset as this provides the required modern samples for Mathieson and Çayönü.
3. Execute the jupyter notebook.
4. Done 🙂

You can now use the dataset to run the benchmark pipeline:

1. First, you need to configure the snakemake pipeline according to the dataset you want to analyze. You can do so by modifying the `config.yaml` file. See the provided `config.yaml` for hints on the allowed settings.
2. Run the pipeline in dry-run mode to see what rules will be run and what files created: `snakemake -n`
3. Finally, run the pipeline using `snakemake -c1`. Note that using `c1` ensures that all rules are executed sequentially. This is required to make sure that all benchmarks are correct and don't influence each other.
4. You can load the resulting `results.parquet` file using pandas:
   ```python
   import pandas as pd
   df = pd.read_parquet("results.parquet")
   ```

### Simulated Datasets
To reproduce the simulated dataset benchmarks, you first need to download the data as instructed in the first few lines of the `benchmark_simulated_data.py` script.
You can then adapt the `# CONFIG` section and set the analysis to `PCA` or `MDS` and run the benchmarks:
```bash
python bebenchmark_simulated_data.py
```    
    
Finally, you can load the resulting `results.parquet` file using pandas:
```python
import pandas as pd
df = pd.read_parquet("results.parquet")
```


## Download the results
Instead of rerunning the analyses yourself, you can also download all (log)files from our experiments via [our lab server](https://cme.h-its.org/exelixis/material/Pandora_supplementary_data.tar.gz).
The directory follows this structure:
```text
├── empirical_datasets [1]
│   ├── pca  # PCA results
│   │   ├── {dataset}  # one such directory for each empirical dataset 
│   │   │   ├── no_convergence  # Logs and (intermediate) results with convergence check disabled
│   │   │   ├── convergence_5p  # Logs and (intermediate) results with 5% convergence tolerance
│   │   │   ├── convergence_1p  # Logs and (intermediate) results with 1% convergence tolerance
│   │   │   ├── convergence_5p_10threads  # Logs and (intermediate) results with 5% convergence tolerance and only 10 threads
│   ├── mds  # MDS results
│   │   ├── {dataset}  # one such directory for each empirical dataset 
│   │   │   ├── no_convergence  # Logs and (intermediate) results with convergence check disabled
│   │   │   ├── convergence_5p  # Logs and (intermediate) results with 5% convergence tolerance
│   │   │   ├── convergence_1p  # Logs and (intermediate) results with 1% convergence tolerance
├── simulated_datasets [2]
│   ├── pca  # PCA results
│   │   ├── {dataset}  # one such directory for each simulated dataset 
│   ├── mds  # MDS results (FST population distance)
│   │   ├── {dataset}  # one such directory for each simulated dataset 
│   ├── mds_sample_dist  # MDS (missing corrected Hamming distance)
│   │   ├── {dataset}  # one such directory for each simulated dataset 
├── sliding_window_hgdp  # Sliding-window results
│   ├── 12_windows
│   ├── 50_windows
```

[1] Empirical datasets includes the following:
- HO-WE (*smartpca* settings: 5 outlier iterations, shrinking disabled)
- HO-WE-no_outlier (*smartpca* settings: 0 outlier iterations, shrinking disabled)
- HO-WE-shrink (*smartpca* settings: 5 outlier iterations, shrinking enabled)
- HO-WE-230
- HO-WE-Cayonu (*smartpca* settings: 5 outlier iterations, shrinking disabled)
- HO-WE-Cayonu-shrink (*smartpca* settings: 5 outlier iterations, shrinking enabled)
- HO-Glob
- Goats
- Sheep
- Panel{i} for i in 1, ..., 13

[2] Simulated datasets includes the following:
- Cline
- Cline_1p (1% missing data)
- Cline_10p (10% missing data)
- Cline_20p (20% missing data)
- Island
- Island_1p (1% missing data)
- Island_10p (10% missing data)
- Island_20p (20% missing data)
- P3
- P3_1p (1% missing data)
- P3_10p (10% missing data)
- P3_20p (20% missing data)
- P3-mig50
- P3-mig50_1p (1% missing data)
- P3-mig50_10p (10% missing data)
- P3-mig50_20p (20% missing data)

Within the individual Pandora result directories, the directory structure follows the specification as stated in the [Pandora documentation](https://pandorageno.readthedocs.io/en/latest/cli_config.html#output-files).


## Publication
The paper explaining the details of Pandora is available as preprint on bioRxiv:

Haag, J., Jordan A. I. & Stamatakis, A. (2024). **Pandora: A Tool to Estimate Dimensionality Reduction Stability of Genotype Data.** *bioRxiv*. [https://doi.org/10.1101/2024.03.14.584962](https://doi.org/10.1101/2024.03.14.584962)
