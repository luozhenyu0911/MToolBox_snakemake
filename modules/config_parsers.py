#!/usr/bin/env python
import os
from typing import List

import pandas as pd
from snakemake.io import expand

def parse_config_tabs(analysis_tab_file=None, reference_tab_file=None, datasets_tab_file=None):
    analysis_tab = pd.read_table("data/analysis.tab", sep = "\t", comment='#')
    reference_tab = (pd.read_table("data/reference_genomes.tab", sep = "\t", comment='#')
                     .set_index("ref_genome_mt", drop=False))
    datasets_tab = pd.read_table("data/datasets.tab", sep = "\t", comment='#')
    return analysis_tab, reference_tab, datasets_tab
    
def get_analysis_species(ref_genome_mt, ref_genome_mt_species_dict=None, config_species=None):
    """ Return a string of species used for analysis.

    Args:
        ref_genome_mt: ref_genome_mt as parsed from analysis_tab 
        config_species: 
    Returns:
        string
    """
    if config_species:
        species = config_species
    else:
        species = ref_genome_mt_species_dict[ref_genome_mt]
    return species

# TODO: infolder is not used anywhere
def get_datasets_for_symlinks(df, sample=None, library=None, d=None,
                              infolder="data/reads", outfolder="data/reads"):
    dataset_file = None
    for row in df.itertuples():
        if (getattr(row, "sample") == sample and
                getattr(row, "library") == int(library)):
            dataset_file = os.path.join(outfolder, getattr(row, d))

    return dataset_file


# TODO: infolder is not used anywhere
def get_symlinks(df, analysis_tab=None,
                 infolder="data/reads", outfolder="data/reads"):
    outpaths = []
    # TODO: convert .iterrows() to .itertuples() for efficiency
    for i, l in df.iterrows():
        if l["sample"] in list(analysis_tab["sample"]):
            outpaths.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}.R1.fastq.gz".format(
                        sample=l["sample"], library=l["library"])
                )
            )
            outpaths.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}.R2.fastq.gz".format(
                        sample=l["sample"], library=l["library"])
                )
            )
    return outpaths


def get_genome_single_vcf_files(df, res_dir="results", ref_genome_mt=None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
            outpaths.append(("{results}/{sample}/{sample}_{ref_genome_mt}_"
                             "{ref_genome_n}.vcf.gz").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=getattr(row, "ref_genome_mt"),
                ref_genome_n=getattr(row, "ref_genome_n")))

    return list(set(outpaths))


# TODO: library is not used anywhere
def get_sample_bamfiles(df, res_dir="results", sample=None, library=None,
                        ref_genome_mt=None, ref_genome_n=None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "sample") == sample:
            bam_file = ("{sample}_{library}_{ref_genome_mt}_"
                        "{ref_genome_n}_OUT-sorted.final.bam").format(
                sample=sample,
                library=getattr(row, "library"),
                ref_genome_mt=ref_genome_mt,
                ref_genome_n=ref_genome_n)
            out_folder = "OUT_{base}".format(
                base=bam_file.replace("_OUT-sorted.final.bam", ""))
            outpaths.append("{results}/{sample}/map/{out_folder}/{bam_file}".format(
                results=res_dir,
                bam_file=bam_file,
                sample=sample,
                out_folder=out_folder))

    return outpaths


def get_genome_single_vcf_index_files(df, res_dir="results", ref_genome_mt=None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
            outpaths.append(
                ("{results}/{sample}/{sample}_{ref_genome_mt}_"
                 "{ref_genome_n}.vcf.gz.csi").format(
                    results=res_dir,
                    sample=getattr(row, "sample"),
                    ref_genome_mt=getattr(row, "ref_genome_mt"),
                    ref_genome_n=getattr(row, "ref_genome_n")))

    return list(set(outpaths))


def get_genome_vcf_files(df: pd.DataFrame,
                         annotation = False,
                         res_dir: str = "results/vcf") -> List[str]:
    """ Return a list of output filenames where VCF files will be stored.

    Args:
        df: input pandas DataFrame
        annotation: are these annotated vcfs?
        res_dir: output directory name

    Returns:
        list of paths
    """
    outpaths = set()
    if annotation:
        annotated = ".annotated"
    else:
        annotated = ""
    # TODO: this is inefficient as it goes through every row and
    #   then removes duplicates, there is a better way for this
    for row in df.itertuples():
        outpaths.add(
            "{results}/{ref_genome_mt}_{ref_genome_n}{annotated}.vcf".format(
                results=res_dir,
                ref_genome_mt=row.ref_genome_mt,
                ref_genome_n=row.ref_genome_n,
                annotated=annotated
            )
        )
    return list(outpaths)


def get_bed_files(df: pd.DataFrame,
                  res_dir: str = "results") -> List[str]:
    """ Return a list of output filenames where BED files will be stored.

    Args:
        df: input pandas DataFrame
        res_dir: output directory name

    Returns:
        list of paths
    """
    outpaths = []
    for row in df.itertuples():
        outpaths.append(
            ("{results}/{sample}/{sample}_"
             "{ref_genome_mt}_{ref_genome_n}.bed").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=row.ref_genome_mt,
                ref_genome_n=row.ref_genome_n
            )
        )
    return outpaths


def get_fasta_files(df: pd.DataFrame,
                    res_dir: str = "results") -> List[str]:
    """ Return a list of output filenames where fasta files will be stored.

    Args:
        df: input pandas DataFrame
        res_dir: output directory name

    Returns:
        list of paths
    """
    outpaths = []
    for row in df.itertuples():
        outpaths.append(
            ("{results}/{sample}/{sample}_"
             "{ref_genome_mt}_{ref_genome_n}.fasta").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=row.ref_genome_mt,
                ref_genome_n=row.ref_genome_n
            )
        )
    return outpaths


def get_haplo_prediction_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append(
            ("{results}/{sample}/{sample}_"
             "{ref_genome_mt}_{ref_genome_n}.csv").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=getattr(row, "ref_genome_mt"),
                ref_genome_n=getattr(row, "ref_genome_n")))
    return outpaths


def get_genome_files(df: pd.DataFrame,
                     ref_genome_mt: str,
                     field: str) -> List[str]:
    """ Return a list of output filenames where fna files will be stored.

    Args:
        df: input pandas DataFrame
        ref_genome_mt: wildcard from snakefile
        field: column name

    Returns:
        list of paths
    """
    return expand(df.loc[ref_genome_mt, field])


def get_mt_genomes(df: pd.DataFrame) -> List[str]:
    """ Return a list of unique mt genome identifiers from the
    given dataframe.

    Args:
        df: input pandas DataFrame

    Returns:
        list of str
    """
    return df["ref_genome_mt"].unique().tolist()


def get_mt_fasta(df, ref_genome_mt, field):
    return df.loc[df['ref_genome_mt'] == ref_genome_mt, field][0]


def fastqc_outputs(datasets_tab: pd.DataFrame,
                   analysis_tab: pd.DataFrame,
                   out: str = "raw") -> List[str]:
    """ Return a list of output filenames where FastQC results will be stored.

    Args:
        datasets_tab: input pandas DataFrame with fastq filenames
        analysis_tab: input pandas DataFrame with analysis details
        out: either 'raw' or 'filtered', determines the output
            directory where FastQC results will be stored

    Returns:
        list of paths
    """
    if out == "raw":
        outfolder = "results/fastqc_raw"
        suffixes = {"R1" : ".R1_fastqc.html", "R2" : ".R2_fastqc.html", "U" : ".U_fastqc.html"}
    elif out == "filtered":
        outfolder = "results/fastqc_filtered"
        suffixes = {"R1" : "_qc_R1_fastqc.html", "R2" : "_qc_R2_fastqc.html", "U" : "_qc_U_fastqc.html"}
    else:
        raise ValueError(f"{out} is not a valid argument")

    fastqc_out = []
    samples = analysis_tab["sample"].tolist()
    for row in datasets_tab.itertuples():
        # TODO: using getattr for the sample column since sample
        #   is already a method name in pandas dataframes, possibly
        #   need to change that column name
        if getattr(row, "sample") in samples:
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}{suffix}".format(
                        sample=getattr(row, "sample"),
                        library=row.library,
                        suffix=suffixes["R1"]
                    )
                )
            )
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}{suffix}".format(
                        sample=getattr(row, "sample"),
                        library=row.library,
                        suffix=suffixes["R2"]
                    )
                )
            )
            if out == "filtered":
                fastqc_out.append(
                    os.path.join(
                        outfolder,
                        "{sample}_{library}{suffix}".format(
                            sample=getattr(row, "sample"),
                            library=row.library,
                            suffix=suffixes["U"]
                        )
                    )
                )
    return fastqc_out
