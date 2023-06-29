import pathlib
import subprocess
import tempfile
import textwrap
import yaml

import pandas as pd

from pandora.custom_types import *

from .config import *


def merge_datasets(
    prefix_ds1: FilePath, prefix_ds2: FilePath, prefix_out: FilePath, redo: bool = False
):
    outfiles = [
        pathlib.Path(f"{prefix_out}.geno"),
        pathlib.Path(f"{prefix_out}.snp"),
        pathlib.Path(f"{prefix_out}.ind"),
    ]
    if all(o.exists() for o in outfiles) and not redo:
        return

    with tempfile.NamedTemporaryFile("w") as tmpfile:
        content = textwrap.dedent(
            f"""
            geno1: {prefix_ds1}.geno
            snp1: {prefix_ds1}.snp
            ind1: {prefix_ds1}.ind
            geno2: {prefix_ds2}.geno
            snp2: {prefix_ds2}.snp
            ind2: {prefix_ds2}.ind
            genooutfilename: {prefix_out}.geno
            snpoutfilename: {prefix_out}.snp
            indoutfilename: {prefix_out}.ind
        """
        )
        tmpfile.write(content)
        tmpfile.flush()
        subprocess.check_output([MERGEIT, "-p", tmpfile.name])


def filter_dataset(
    prefix_in: FilePath, prefix_out: FilePath, poplistname: FilePath, redo: bool = False
):
    outfiles = [
        pathlib.Path(f"{prefix_out}.geno"),
        pathlib.Path(f"{prefix_out}.snp"),
        pathlib.Path(f"{prefix_out}.ind"),
    ]
    if all(o.exists() for o in outfiles) and not redo:
        return

    with tempfile.NamedTemporaryFile("w") as tmpfile:
        content = textwrap.dedent(
            f"""
            genotypename: {prefix_in}.geno
            snpname: {prefix_in}.snp
            indivname: {prefix_in}.ind
            poplistname: {poplistname}
            genotypeoutname: {prefix_out}.geno
            snpoutname: {prefix_out}.snp
            indivoutname: {prefix_out}.ind
        """
        )
        tmpfile.write(content)
        tmpfile.flush()
        subprocess.check_output([CONVERTF, "-p", tmpfile.name])


def indfile_to_dataframe(indfile: FilePath):
    cols = ["sampleID", "sex", "population"]
    df = pd.read_table(indfile, delimiter="\s+", skipinitialspace=True, names=cols)
    for c in cols:
        df[c] = df[c].str.strip()
    return df


def get_pop_set_from_string(s: str):
    return [v.strip() for v in s.split(",")]


def save_pop_set(pop_set: List[str], outfile: FilePath):
    with outfile.open("w") as f:
        f.write("\n".join(pop_set))


def save_pandora_config(
    dataset_prefix: FilePath,
    result_dir: FilePath,
    smartpca_optional_settings: Dict[str, str],
    configfile: FilePath,
    **extra
):
    config = dict(
        dataset_prefix=str(dataset_prefix),
        file_format="ANCESTRYMAP",
        result_dir=str(result_dir),
        bootstrap=True,
        n_bootstraps=100,
        keep_bootstraps=False,
        seed=0,
        n_pcs=20,
        smartpca=SMARTPCA,
        convertf=CONVERTF,
        smartpca_optional_settings=smartpca_optional_settings,
        **extra
    )

    config_yaml = yaml.safe_dump(config)
    configfile.open(mode="w").write(config_yaml)
