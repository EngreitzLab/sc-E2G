import os
import shutil

import click
import pandas as pd


def create_intermediate_dir(results_dir):
    intermediate_dir = os.path.join(results_dir, "intermediate_files")
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    return intermediate_dir

def delete_intermediate_dir(intermediate_dir):
    shutil.rmtree(intermediate_dir)


def add_midpoint(df):
    df["midpoint"] = ((df["start"] + df["end"]) / 2).astype("int")


@click.command()
@click.option("--enhancer_list")
@click.option("--abc_predictions")
@click.option("--ref_gene_tss")
@click.option("--chr_sizes")
@click.option("--results_dir")
def main(enhancer_list, abc_predictions, ref_gene_tss, chr_sizes, results_dir):
    intermediate_dir = create_intermediate_dir(results_dir)

    ##
    pred_df = pd.read_csv(abc_predictions, sep="\t")
    pred_df = pred_df[pred_df["class"] != "promoter"]
    if len(pred_df) == 0:
        raise Exception("Did not find any enhancers in the Predictions file")
    add_midpoint(pred_df)

    determine_num_candidate_enh_gene(pred_df, results_dir)
    determine_num_tss_enh_gene(pred_df, ref_gene_tss, results_dir, intermediate_dir)
    generate_num_sum_enhancers(abc_predictions, enhancer_list, chr_sizes, results_dir, intermediate_dir)

    # Intermediate directory takes up a lot of space. Delete after usage
    delete_intermediate_dir(intermediate_dir)

def _populate_enhancer_count_from_tss(df, enhancers, is_upstream):
    enh_indexes = enhancers.index
    if is_upstream:
        # start counting from the enhancer closest to TSS
        enh_indexes = reversed(enh_indexes)
    
    count_from_tss = 0
    for enh_idx in enh_indexes:
        count_from_tss += 1
        df.loc[enh_idx, "NumCandidateEnhGene"] = count_from_tss


def determine_num_candidate_enh_gene(pred_df, results_dir):
    # Need df to be sorted by midpoint for each chromosome
    df = pred_df.sort_values(by=['chr', 'midpoint'], ascending=True).reset_index(drop=True)

    gene_groups = df.groupby(["TargetGene", "TargetGeneTSS"])
    for (gene, tss), indexes in gene_groups.groups.items():
        enhancers = df.loc[indexes]
        upstream_enh = enhancers[enhancers["midpoint"] < tss]
        downstream_enh = enhancers[enhancers["midpoint"] > tss]
        _populate_enhancer_count_from_tss(df, upstream_enh, is_upstream=True)
        _populate_enhancer_count_from_tss(df, downstream_enh, is_upstream=False)

    df["NumCandidateEnhGene"] = df["NumCandidateEnhGene"].astype("int")
    df[["name", "TargetGene", "NumCandidateEnhGene"]].to_csv(
        os.path.join(results_dir, "NumCandidateEnhGene.tsv"),
        sep="\t",
        index=False,
    )
    print("Saved num candidate enhancers")


def determine_num_tss_enh_gene(pred_df, ref_gene_tss, results_dir, intermediate_dir):
    ##### Make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
    pred_df["new_end"] = (pred_df["midpoint"] + pred_df["distance"]).astype("int")

    ## If gene is located upstream of enhancer, modify the start to be the beginning of the TargetGeneTSS and the end be the midpoint of the enhancer
    downstream_enh = pred_df[pred_df["TargetGeneTSS"] < pred_df["midpoint"]]
    pred_df.loc[downstream_enh.index, "new_end"] = pred_df.loc[
        downstream_enh.index, "end"
    ]
    pred_df.loc[downstream_enh.index, "start"] = pred_df.loc[
        downstream_enh.index, "TargetGeneTSS"
    ]

    # This file will be used to intersect gene annotations to count How many protein-coding TSSs away is the enhancer from the promoter?  (i.e., how many protein-coding gene TSSs are located between the enhancer and promoter?  0 = closest TSS)
    # Intersect this file with gene TSS
    extended_enhregions = os.path.join(intermediate_dir, "EnhancerRegions_extended.txt")
    pred_df[["chr", "start", "new_end", "name", "TargetGene"]].to_csv(
        extended_enhregions, sep="\t", index=False
    )
    ## Intersect midpoint of enhancer to target gene regions with GeneTSS
    ## This will include overlaps with the TargetGene
    filename = os.path.join(
        intermediate_dir, "EnhancerRegions_extended_RefSeq.intersected.tsv.gz"
    )
    os.system(
        "sed '1d' {} | bedtools intersect -a stdin -b {} -wa -wb | cut -f4,5 | pigz > {}".format(
            extended_enhregions, ref_gene_tss, filename
        )
    )
    print("Reading in {}".format(filename))
    ## Read in intersected midpoint of enhancer to target gene regions with GeneTSS file
    predictions = pd.read_csv(filename, sep="\t", names=["class", "gene"])
    # Calculate the number of TSS regions that fall within the enhancer to target gene regions.
    num_tss_between_enh_and_gene = (
        predictions.groupby(["class", "gene"]).size().reset_index()
    )
    num_tss_between_enh_and_gene.to_csv(
        os.path.join(results_dir, "NumTSSEnhGene.tsv"),
        sep="\t",
        index=False,
    )
    print("Saved num TSS between enh and gene")


def generate_num_sum_enhancers(
    pred_file, enhancer_list, chr_sizes, results_dir, intermediate_dir
):
    enh_list_df = pd.read_csv(
        enhancer_list, sep="\t", usecols=["chr", "start", "end", "name"]
    )
    add_midpoint(enh_list_df)

    ############ Generate Num/Sum Enhancers within 5kb/10kb ############
    ##### Make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
    enh_list_midpoint_file = os.path.join(intermediate_dir, "EnhancerList_midpoint.bed")
    enh_list_df[["chr", "midpoint", "midpoint", "name"]].to_csv(
        enh_list_midpoint_file,
        sep="\t",
        index=False,
        header=False,
    )
    slop_5kb_file = os.path.join(intermediate_dir, "EnhancerRegions.slop_5kb.bed")
    os.system(
        "bedtools slop -b 5000 -i {} -g {} > {}".format(
            enh_list_midpoint_file, chr_sizes, slop_5kb_file
        )
    )
    slop_10kb_file = os.path.join(intermediate_dir, "EnhancerRegions.slop_10kb.bed")
    os.system(
        "bedtools slop -b 10000 -i {} -g {} > {}".format(
            enh_list_midpoint_file, chr_sizes, slop_10kb_file
        )
    )
    prediction_slim_file = os.path.join(
        intermediate_dir, "EnhancerPredictionsAllPutative.slim.bed"
    )
    os.system(
        "zcat {} | csvtk cut -t -f chr,start,end,name,activity_base | sed '1d' > {}".format(
            pred_file, prediction_slim_file
        )
    )
    os.system(
        "bedtools intersect -a {} -b {} -wa -wb > {} | sort -u ".format(
            slop_5kb_file,
            prediction_slim_file,
            os.path.join(intermediate_dir, "NumEnhancers5kb.txt"),
        )
    )
    os.system(
        "bedtools intersect -a {} -b {} -wa -wb | sort -u > {}".format(
            slop_10kb_file,
            prediction_slim_file,
            os.path.join(intermediate_dir, "NumEnhancers10kb.txt"),
        )
    )
    data = pd.read_csv(
        os.path.join(intermediate_dir, "NumEnhancers5kb.txt"),
        sep="\t",
        header=None,
    )
    data1 = data[data[3] != data[7]]
    data2 = data1.groupby([3]).size().reset_index(name="count")
    data3 = data1.groupby([3])[8].sum().reset_index(name="sum")
    data2.to_csv(
        os.path.join(results_dir, "NumEnhancersEG5kb.txt"),
        sep="\t",
        header=False,
        index=False,
    )
    data3.to_csv(
        os.path.join(results_dir, "SumEnhancersEG5kb.txt"),
        sep="\t",
        header=False,
        index=False,
    )
    data = pd.read_csv(
        os.path.join(intermediate_dir, "NumEnhancers10kb.txt"),
        sep="\t",
        header=None,
    )
    data1 = data[data[3] != data[7]]
    data2 = data1.groupby([3]).size().reset_index(name="count")
    data3 = data1.groupby([3])[8].sum().reset_index(name="sum")
    data2.to_csv(
        os.path.join(results_dir, "NumEnhancersEG10kb.txt"),
        sep="\t",
        header=False,
        index=False,
    )
    data3.to_csv(
        os.path.join(results_dir, "SumEnhancersEG10kb.txt"),
        sep="\t",
        header=False,
        index=False,
    )


if __name__ == "__main__":
    main()
