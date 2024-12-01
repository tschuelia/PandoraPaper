import allel
import numpy as np
import pandas as pd
import pathlib
import yaml


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
        seed: int
):
    config = {
        'dataset_prefix': str(eigen_prefix.absolute()),
        'result_dir': str(result_dir.absolute()),
        'n_replicates': n_bootstraps,
        'threads': n_threads,
        'seed': seed,
        'plot_results': True,
        'kmeans_k': 2,
        'keep_replicates': False,
        'smartpca_optional_settings': {
            'numoutlieriter': 0
        }
    }

    yaml.dump(config, config_file.open("w"))


def get_pandora_results(pandora_log: pathlib.Path, support_values: pathlib.Path):
    ps = None
    with pandora_log.open() as f:
        for line in f:
            if line.startswith("Pandora Stability"):
                ps = float(line.rsplit(maxsplit=1)[1])

    psvs = pd.read_csv(support_values, index_col=0)
    psvs["sample_id"] = psvs.index
    psvs.reset_index(drop=True, inplace=True)
    psvs["PS"] = ps
    return psvs


def get_n_samples_n_snps(ind_file: pathlib.Path, snp_file: pathlib.Path):
    n_samples = sum(1 for line in ind_file.open() if line.strip())
    n_snps = sum(1 for line in snp_file.open() if line.strip())
    return n_samples, n_snps


def get_missing_per_sample(geno_file: pathlib.Path, ind_file: pathlib.Path):
    sample_ids = []
    populations = []

    for line in ind_file.open():
        sid, _, pop = line.strip().split()
        sample_ids.append(sid.strip())
        populations.append(pop.strip())

    geno_data = []

    for line in geno_file.open():
        geno_data.append([int(c) for c in line.strip()])

    geno_data = np.asarray(geno_data).T

    assert geno_data.shape[0] == len(sample_ids) == len(populations), "Incorrect dimensions of geno data"

    missing = np.sum(geno_data == 9, axis=1) / geno_data.shape[1]
    return pd.DataFrame({"sample_id": sample_ids, "population": populations, "missing_per_sample": missing})


def collect_data_for_experiment(
        pandora_result: pathlib.Path,
        support_values: pathlib.Path,
        ind_file: pathlib.Path,
        snp_file: pathlib.Path,
        geno_file: pathlib.Path,
        missing_fraction: float,
        ld_pruned: bool
):
    results = get_pandora_results(pandora_result, support_values)
    results["missing"] = missing_fraction
    results["n_samples"], results["n_snps"] = get_n_samples_n_snps(ind_file, snp_file)
    results["ld_pruned"] = ld_pruned

    missing_per_sample = get_missing_per_sample(geno_file, ind_file)
    results = pd.merge(results, missing_per_sample, on="sample_id")
    return results