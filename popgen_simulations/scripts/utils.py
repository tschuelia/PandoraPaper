import allel
import datetime
import math
import pandas as pd
import pathlib
import shutil
import textwrap
import time
import uuid
import yaml

from snakemake import shell

from pandora.bootstrap import bootstrap_and_embed_multiple_numpy
from pandora.custom_types import EmbeddingAlgorithm
from pandora.converter import get_filenames, FileFormat
from pandora.dataset import numpy_dataset_from_eigenfiles
from pandora.embedding_comparison import BatchEmbeddingComparison


def write_ind_file(ts, ind_file: pathlib.Path):
    pop_mapping = {pop.id: pop.metadata["description"].replace(" ", "_") for pop in ts.populations()}

    with ind_file.open("w") as f:
        for ind, pop in zip(ts.individuals(), ts.individuals_population):
            f.write(f"indiv_{ind.id} U {pop_mapping[pop]}\n")


def write_geno_file(ts, geno_file: pathlib.Path):
    haps = allel.HaplotypeArray(ts.genotype_matrix())
    gns = haps.to_genotypes(ploidy=2)
    all_geno = gns.to_n_ref(fill=9)
    with geno_file.open("w") as f:
        for row in all_geno:
            f.write("".join(row.astype(str)) + "\n")


def write_snp_file(ts, snp_file: pathlib.Path):
    with snp_file.open("w") as f:
        for var in ts.variants():
            ref, alt, *_ = var.alleles
            f.write(f"snp_{var.site.id} 1 0.0 {int(var.position)} {ref} {alt}\n")


def write_pandora_config(
        eigen_prefix: pathlib.Path,
        result_dir: pathlib.Path,
        config_file: pathlib.Path,
        n_bootstraps: int,
        n_threads: int,
        seed: int,
        smartpca: pathlib.Path,
        n_populations: int,
        algo: str,
        n_components: int,
        bootstrap_convergence_check: bool,
        convergence_tolerance: float
):
    config = {
        "dataset_prefix": str(eigen_prefix.absolute()),
        "result_dir": str(result_dir.absolute()),
        "n_replicates": n_bootstraps,
        "threads": n_threads,
        "seed": seed,
        "kmeans_k": n_populations,
        "smartpca": str(smartpca.absolute()),
        "smartpca_optional_settings": {
            "numoutlieriter": 0
        },
        "embedding_algorithm": algo,
        "n_components": n_components,
        "bootstrap_convergence_check": bootstrap_convergence_check,
        "bootstrap_convergence_tolerance": convergence_tolerance
    }

    yaml.dump(config, config_file.open("w"))


def get_psv_summary(support_values: pathlib.Path):
    psvs = pd.read_csv(support_values, index_col=0)
    psvs["sample_id"] = psvs.index
    psvs.reset_index(drop=True, inplace=True)
    return {
        "psv_mean": psvs["PSV"].mean(),
        "psv_median": psvs["PSV"].median(),
        "psv_std": psvs["PSV"].std(),
        "psv_min": psvs["PSV"].min(),
        "psv_max": psvs["PSV"].max(),
        "psv_q10": psvs["PSV"].quantile(0.10),
        "psv_q5": psvs["PSV"].quantile(0.05),
        "psv_q1": psvs["PSV"].quantile(0.01),
        "psvs": [psvs["PSV"].values],
    }


def get_pandora_results(pandora_log: pathlib.Path, support_values: pathlib.Path):
    data = {}

    for line in pandora_log.open().readlines():
        line = line.strip()
        if line.startswith("Pandora Stability"):
            data["PS"] = float(line.split(":")[1].strip())
        elif line.startswith("Total runtime"):
            # Total runtime: 0:53:06 (3186 seconds)
            _, runtime = line.rsplit("(")
            runtime = runtime.strip("seconds)")
            data["runtime"] = int(runtime)
        elif line.startswith("> Number of replicates computed"):
            data["n_replicates"] = int(line.split(":")[1].strip())
        elif line.startswith("> Number of Kmeans clusters"):
            data["n_clusters"] = int(line.split(":")[1].strip())

    data.update(get_psv_summary(support_values))
    data = pd.DataFrame(data, index=[0])
    return data


def get_n_samples_n_snps(ind_file: pathlib.Path, snp_file: pathlib.Path):
    n_samples = sum(1 for line in ind_file.open() if line.strip())
    n_snps = sum(1 for line in snp_file.open() if line.strip())
    return n_samples, n_snps


def collect_data_for_experiment(
        pandora_log: pathlib.Path,
        support_values: pathlib.Path,
        ind_file: pathlib.Path,
        snp_file: pathlib.Path,
        geno_file: pathlib.Path,
        missing_fraction: float,
        noise_fraction: float,
        ld_pruned: bool,
        convergence: bool,
        convergence_tolerance: float,
):
    results = get_pandora_results(pandora_log, support_values)
    results["missing"] = missing_fraction
    results["noise"] = noise_fraction
    results["n_samples"], results["n_snps"] = get_n_samples_n_snps(ind_file, snp_file)
    results["ld_pruned"] = ld_pruned
    results["convergence"] = convergence
    results["convergence_tolerance"] = convergence_tolerance

    return results


def execute_pandora_config_pca(config):
    dataset_prefix = pathlib.Path(config["dataset_prefix"])
    geno, snp, ind = get_filenames(dataset_prefix, FileFormat.EIGENSTRAT)

    result_dir = pathlib.Path(config["result_dir"])

    # 1: Copy the input files to the local directory (required for correct benchmarking)
    local_tmp_dir = pathlib.Path(f"{dataset_prefix.name}_{uuid.uuid4()}")
    local_tmp_dir.mkdir()

    shutil.copy(geno, local_tmp_dir)
    shutil.copy(snp, local_tmp_dir)
    shutil.copy(ind, local_tmp_dir)
    local_dataset_prefix = local_tmp_dir / (geno.stem)

    config["dataset_prefix"] = str(local_dataset_prefix.absolute())

    local_result_dir = local_tmp_dir / result_dir.name
    outfile = local_tmp_dir / "pandora.out"
    errorfile = local_tmp_dir / "pandora.err"

    config["result_dir"] = str(local_result_dir.absolute())

    local_config_file = local_tmp_dir / "pandora_config.yaml"
    yaml.dump(config, local_config_file.open("w"))

    # 2. Run Pandora with the corrected config and write the results to a local directory
    shell(f"pandora -c {local_config_file} > {outfile} 2> {errorfile}")

    # 3. Move the files back to the original target destination
    shell(f"mv {local_result_dir}/* {result_dir}/")
    shell(f"rm -rf {local_tmp_dir}")


def execute_pandora_config_mds(config):
    dataset_prefix = pathlib.Path(config["dataset_prefix"])
    results_dir = pathlib.Path(config["result_dir"])
    logfile = results_dir / "pandora.log"
    sample_support_values_csv = results_dir / "pandora.supportValues.csv"
    result_file = results_dir / "pandora.txt"

    with logfile.open("w") as f:
        f.write(f"Starting analysis: {dataset_prefix} (MDS)\n")
        _start = time.perf_counter()
        dataset = numpy_dataset_from_eigenfiles(dataset_prefix)

        f.write("Bootstrap and embed multiple\n")
        f.write(f"Number of bootstraps: {config['n_replicates']}\n")
        f.write(f"Convergence: {config['bootstrap_convergence_check']} with tolerance {config['bootstrap_convergence_tolerance']}\n")
        bootstraps = bootstrap_and_embed_multiple_numpy(
            dataset=dataset,
            n_bootstraps=config["n_replicates"],
            embedding=EmbeddingAlgorithm.MDS,
            n_components=config["n_components"],
            seed=config["seed"],
            threads=config["threads"],
            bootstrap_convergence_check=config["bootstrap_convergence_check"],
            bootstrap_convergence_tolerance=config["bootstrap_convergence_tolerance"]
        )

        f.write(f"Number of bootstraps: {len(bootstraps)}\n")
        f.write("Compare embeddings\n")

        comp = BatchEmbeddingComparison([b.mds for b in bootstraps])
        ps = comp.compare()
        psvs = comp.get_sample_support_values(threads=config["threads"])

        _end = time.perf_counter()
        runtime = math.ceil(_end - _start)

        f.write("Done\n")

        psvs.to_csv(sample_support_values_csv)

        f.write(f"Pandora Stability: {ps}\n")
        f.write(f"Total runtime: {datetime.timedelta(seconds=runtime)} ({runtime} seconds)\n")
        f.write(f"> Number of replicates computed: {len(bootstraps)}\n")
        f.write(f"> Number of Kmeans clusters: {config['kmeans_k']}\n")

    results_string = textwrap.dedent(
        f"""
            > Performed Analysis: MDS
            > Number of replicates computed: {len(bootstraps)}
            > Number of Kmeans clusters: {config['kmeans_k']}

            ------------------
            Results
            ------------------
            Pandora Stability: {round(ps, 2)}"""
    )

    result_file.write_text(results_string)

def execute_pandora_config(config_file: pathlib.Path):
    config = yaml.safe_load(config_file.open())

    if config["embedding_algorithm"] == "PCA":
        execute_pandora_config_pca(config)
    else:
        execute_pandora_config_mds(config)