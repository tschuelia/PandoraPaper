import pathlib
import subprocess
import tempfile
import textwrap

import pandas as pd

MERGEIT = "/Users/julia/micromamba/envs/pca_micromamba/bin/mergeit"
CONVERTF = "/Users/julia/micromamba/envs/pca_micromamba/bin/convertf"


def merge_datasets(prefix_ds1, prefix_ds2, prefix_out, redo=False):
    outfiles = [
        pathlib.Path(f"{prefix_out}.geno"),
        pathlib.Path(f"{prefix_out}.snp"),
        pathlib.Path(f"{prefix_out}.ind"),
    ]
    if all(o.exists() for o in outfiles) and not redo:
        return

    with tempfile.NamedTemporaryFile("w") as tmpfile:
        content = textwrap.dedent(f"""
            geno1: {prefix_ds1}.geno
            snp1: {prefix_ds1}.snp
            ind1: {prefix_ds1}.ind
            geno2: {prefix_ds2}.geno
            snp2: {prefix_ds2}.snp
            ind2: {prefix_ds2}.ind
            genooutfilename: {prefix_out}.geno
            snpoutfilename: {prefix_out}.snp
            indoutfilename: {prefix_out}.ind
        """)
        tmpfile.write(content)
        tmpfile.flush()
        subprocess.check_output([MERGEIT, "-p", tmpfile.name])


def filter_dataset(prefix_in, prefix_out, poplistname, redo=False):
    outfiles = [
        pathlib.Path(f"{prefix_out}.geno"),
        pathlib.Path(f"{prefix_out}.snp"),
        pathlib.Path(f"{prefix_out}.ind"),
    ]
    if all(o.exists() for o in outfiles) and not redo:
        return

    with tempfile.NamedTemporaryFile("w") as tmpfile:
        content = textwrap.dedent(f"""
            genotypename: {prefix_in}.geno
            snpname: {prefix_in}.snp
            indivname: {prefix_in}.ind
            poplistname: {poplistname}
            genotypeoutname: {prefix_out}.geno
            snpoutname: {prefix_out}.snp
            indivoutname: {prefix_out}.ind
        """)
        tmpfile.write(content)
        tmpfile.flush()
        subprocess.check_output([CONVERTF, "-p", tmpfile.name])


def indfile_to_dataframe(indfile):
    cols = ["sampleID", "sex", "population"]
    df = pd.read_table(indfile, delimiter="\s+", skipinitialspace=True, names=cols)
    for c in cols:
        df[c] = df[c].str.strip()
    return df


def get_pop_set_from_string(s):
    return [v.strip() for v in s.split(",")]


def save_pop_set(pop_set, outfile):
    with outfile.open("w") as f:
        f.write("\n".join(pop_set))


