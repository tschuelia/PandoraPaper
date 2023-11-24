"""
# Simulated Data Benchmark

Data obtained from:
Yi, X., & Latch, E. K. (2022). Nonrandom missing data can bias Principal Component Analysis inference of population genetic structure.
Molecular Ecology Resources, 22, 602– 611. https://doi.org/10.1111/1755-0998.13498.

To download the required data, run the following set of commands in your terminal:

cd ../datasets && mkdir simulated_data  && cd simulated_data
wget https://github.com/xuelingyi/missing_data_PCA/raw/main/simulation.zip && unzip simulation.zip && rm simulation.zip
wget https://github.com/xuelingyi/missing_data_PCA/raw/main/MISSdata.zip && unzip MISSdata.zip && rm MISSdata.zip
"""
import pathlib
import pickle
import time
from typing import Any, Callable, Dict, List

import pandas as pd
import pyreadr
from pandora.bootstrap import bootstrap_and_embed_multiple_numpy
from pandora.custom_types import EmbeddingAlgorithm
from pandora.dataset import NumpyDataset
from pandora.distance_metrics import fst_population_distance
from pandora.embedding_comparison import BatchEmbeddingComparison
from scipy.stats import pearsonr

# CONFIG
N_REPLICATES = 100
N_THREADS = 20
EMBEDDING = EmbeddingAlgorithm.PCA
REDO = False
IDENTIFIER = "pca"
PCA_N_COMPONENTS = 2
# CONFIG END

SIMULATED_DATA = pathlib.Path("../datasets/simulated_data/simulation")
SIMULATED_MISSING_DATA = pathlib.Path("../datasets/simulated_data/MISSdata")

OUTDIR = pathlib.Path("results_benchmark") / "simulated_data" / IDENTIFIER
OUTDIR.mkdir(exist_ok=True, parents=True)


def _pop_and_id(joint_name: str) -> (str, str):
    # joint_name is of form "pop1.indiv1"
    # since indivX is not unique, so we actually keep the joint name as sample ID
    population, _ = joint_name.split(".")
    return population.strip(), joint_name.strip()


def load_simulated_data(r_path: pathlib.Path) -> NumpyDataset:
    data = pyreadr.read_r(str(r_path))
    raw_data = data["SNP"]
    populations, sample_ids = zip(*[_pop_and_id(name) for name in raw_data.index])
    geno_data = raw_data.to_numpy()
    dataset = NumpyDataset(geno_data, pd.Series(sample_ids), pd.Series(populations))
    return dataset


def load_missing_data(r_path: pathlib.Path) -> NumpyDataset:
    data = pyreadr.read_r(str(r_path))
    raw_data = data["rand"]
    populations, sample_ids = zip(*[_pop_and_id(name) for name in raw_data.index])
    geno_data = raw_data.to_numpy()
    dataset = NumpyDataset(geno_data, pd.Series(sample_ids), pd.Series(populations))
    return dataset


def _get_bootstraps_pca(
    dataset: NumpyDataset,
    threads: int,
    bootstrap_convergence_check: bool,
    bootstrap_convergence_confidence_level: float,
) -> List[NumpyDataset]:
    return bootstrap_and_embed_multiple_numpy(
        dataset=dataset,
        n_bootstraps=N_REPLICATES,
        embedding=EmbeddingAlgorithm.PCA,
        n_components=PCA_N_COMPONENTS,
        seed=0,
        threads=threads,
        bootstrap_convergence_check=bootstrap_convergence_check,
        bootstrap_convergence_confidence_level=bootstrap_convergence_confidence_level,
    )


def _get_bootstraps_mds(
    dataset: NumpyDataset,
    threads: int,
    bootstrap_convergence_check: bool,
    bootstrap_convergence_confidence_level: float,
) -> List[NumpyDataset]:
    return bootstrap_and_embed_multiple_numpy(
        dataset=dataset,
        n_bootstraps=N_REPLICATES,
        embedding=EmbeddingAlgorithm.MDS,
        n_components=2,
        seed=0,
        threads=threads,
        distance_metric=fst_population_distance,
        bootstrap_convergence_check=bootstrap_convergence_check,
        bootstrap_convergence_confidence_level=bootstrap_convergence_confidence_level,
    )


def _compare_support_values(psv_no_convergence, psv_with_convergence):
    shared_sample_ids = list(
        set(psv_no_convergence.index) & set(psv_with_convergence.index)
    )

    psv_no_convergence = psv_no_convergence.loc[shared_sample_ids]
    psv_with_convergence = psv_with_convergence.loc[shared_sample_ids]

    return pearsonr(psv_no_convergence.values, psv_with_convergence.values)


def _get_bootstraps(dataset: NumpyDataset, threads: int, bootstrap_convergence_check: bool,
    bootstrap_convergence_confidence_level: float):
    if EMBEDDING == EmbeddingAlgorithm.PCA:
        bootstraps = _get_bootstraps_pca(
            dataset,
            threads,
            bootstrap_convergence_check,
            bootstrap_convergence_confidence_level,
        )
    elif EMBEDDING == EmbeddingAlgorithm.MDS:
        bootstraps = _get_bootstraps_mds(
            dataset,
            threads,
            bootstrap_convergence_check,
            bootstrap_convergence_confidence_level,
        )
    else:
        raise ValueError("Unsupported Embedding: ", EMBEDDING)

    return bootstraps


def run_pandora(
    dataset: NumpyDataset,
    threads: int,
    bootstrap_convergence_check: bool,
    bootstrap_convergence_confidence_level: float,
    key_suffix: str,
    outdir: pathlib.Path
) -> (Dict[str, Any], pd.Series):
    data = dict()

    bootstrap_pickle = outdir / f"bootstrap_{key_suffix}.pickle"
    support_values_csv = outdir / f"support_values_{key_suffix}.csv"
    runtime_log = outdir / f"runtime_{key_suffix}.txt"

    do_bootstrap = REDO or not bootstrap_pickle.exists() or not runtime_log.exists() or not support_values_csv.exists()
    bs_start = time.perf_counter()

    if do_bootstrap:
        bootstraps = _get_bootstraps(dataset, threads, bootstrap_convergence_check, bootstrap_convergence_confidence_level)
    else:
        bootstraps = pickle.load(bootstrap_pickle.open("rb"))

    if EMBEDDING == EmbeddingAlgorithm.PCA:
        comparison = BatchEmbeddingComparison([b.pca for b in bootstraps])
    elif EMBEDDING == EmbeddingAlgorithm.MDS:
        comparison = BatchEmbeddingComparison([b.mds for b in bootstraps])
    else:
        raise ValueError("Unsupported Embedding: ", EMBEDDING)

    stability = comparison.compare(threads=threads)
    cluster_stability = comparison.compare_clustering(threads=threads, kmeans_k=3)
    support_values = comparison.get_sample_support_values(threads=threads)
    bs_end = time.perf_counter()

    if do_bootstrap:
        # only use the measured runtime if we actually had to compute the bootstraps
        bs_runtime = bs_end - bs_start
        runtime_log.write_text(str(bs_runtime))
    else:
        # otherwise load it from data
        bs_runtime = float(runtime_log.open().read().strip())

    data[f"stability_{key_suffix}"] = stability
    data[f"cluster_stability_{key_suffix}"] = cluster_stability
    data[f"n_bootstraps_{key_suffix}"] = len(bootstraps)
    data[f"threads_{key_suffix}"] = threads
    data[f"level_{key_suffix}"] = bootstrap_convergence_confidence_level
    data[f"runtime_{key_suffix}"] = bs_runtime

    # store the support values, as we will need them for subsequent analyses
    support_values.to_csv(support_values_csv)
    return data, support_values


def run_dataset(
    dataset_r_path: pathlib.Path,
    load_function: Callable[[pathlib.Path], NumpyDataset],
    outdir: pathlib.Path,
    redo: bool = False
) -> pd.DataFrame:
    """
    Perform and benchmark the following runs:
    1. No convergence check, compute all N_REPLICATES bootstraps
    2. Convergence check with low confidence limit (0.05) and N_THREADS
    3. Convergence check with high confidence limit (0.01) and N_THREADS
    4. Convergence check with low confidence limit (0.05) and N_THREADS // 2
    """
    outfile = outdir / "results.parquet"

    if outfile.exists() and not redo:
        return pd.read_parquet(outfile)

    dataset = load_function(dataset_r_path)
    results = {
        "dataset": dataset_r_path.stem,
        "n_snps": dataset.input_data.shape[1],
        "n_samples": dataset.sample_ids.shape[0],
        "n_populations": dataset.populations.unique().shape[0],
    }

    # baseline results
    results_no_convergence, psv_no_convergence = run_pandora(
        dataset=dataset,
        threads=N_THREADS,
        bootstrap_convergence_check=False,
        bootstrap_convergence_confidence_level=0.05,
        key_suffix="no_conv",
        outdir=outdir
    )
    results.update(results_no_convergence)

    # results lower confidence level 0.05
    results_with_convergence_low, psv_with_convergence_low = run_pandora(
        dataset=dataset,
        threads=N_THREADS,
        bootstrap_convergence_check=True,
        bootstrap_convergence_confidence_level=0.05,
        key_suffix="conv_low",
        outdir=outdir
    )
    results.update(results_with_convergence_low)
    results["speedup_conv_low"] = (
        results["runtime_no_conv"] / results["runtime_conv_low"]
    )
    results["ps_deviation_low"] = (
        abs(results["stability_no_conv"] - results["stability_conv_low"])
        / results["stability_no_conv"]
    )

    psv_corr_low = _compare_support_values(psv_no_convergence, psv_with_convergence_low)
    results["pearson_correlation_statistic_conv_low"] = psv_corr_low.statistic
    results["pearson_correlation_pvalue_conv_low"] = psv_corr_low.pvalue

    # results higher confidence level 0.01
    results_with_convergence_high, psv_with_convergence_high = run_pandora(
        dataset=dataset,
        threads=N_THREADS,
        bootstrap_convergence_check=True,
        bootstrap_convergence_confidence_level=0.01,
        key_suffix="conv_high",
        outdir=outdir
    )
    results.update(results_with_convergence_high)
    results["speedup_conv_high"] = (
        results["runtime_no_conv"] / results["runtime_conv_high"]
    )
    results["ps_deviation_high"] = (
        abs(results["stability_no_conv"] - results["stability_conv_high"])
        / results["stability_no_conv"]
    )

    psv_corr_high = _compare_support_values(
        psv_no_convergence, psv_with_convergence_high
    )
    results["pearson_correlation_statistic_conv_high"] = psv_corr_high.statistic
    results["pearson_correlation_pvalue_conv_high"] = psv_corr_high.pvalue

    # results lower confidence level 0.05 and fewer cores
    results_with_convergence_low_fewer_cores, _ = run_pandora(
        dataset=dataset,
        threads=N_THREADS // 2,
        bootstrap_convergence_check=True,
        bootstrap_convergence_confidence_level=0.05,
        key_suffix="conv_low_fewer_cores",
        outdir=outdir
    )
    results.update(results_with_convergence_low_fewer_cores)
    results["ps_deviation_low_fewer_cores"] = (
        abs(results["stability_conv_low"] - results["stability_conv_low_fewer_cores"])
        / results["stability_conv_low"]
    )

    # write data as dataframe parquet
    results = pd.DataFrame(data=results, index=[0])
    results.to_parquet(outfile)
    return results


if __name__ == "__main__":
    dfs = []
    for directory, load_method, exclusion_criterion in [
        (SIMULATED_DATA, load_simulated_data, lambda pth: "replicate" in pth.name),
        (SIMULATED_MISSING_DATA, load_missing_data, lambda pth: "rand" not in pth.name)
    ]:
        for pth in directory.iterdir():
            if exclusion_criterion(pth):
                continue
            print("RUNNING: ", pth.stem)
            outdir = OUTDIR / pth.stem
            outdir.mkdir(exist_ok=True, parents=True)
            try:
                results = run_dataset(pth, load_method, outdir, REDO)
                dfs.append(results)
            except Exception as e:
                print("error running ", pth.stem, e)

    df = pd.concat(dfs)
    df.to_parquet(OUTDIR / f"results_{IDENTIFIER}.parquet")