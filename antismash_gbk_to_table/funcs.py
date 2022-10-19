#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from collections import OrderedDict
import csv
from pathlib import Path
from typing import Dict, List, Union

OUT_COLUMNS = [
    "locus",
    "tool",
    "type",
    "start",
    "end",
    "strand",
    "category",
    "contig_edge",
    "product",
    "protocluster_number",
    "candidate_cluster_number",
    "candidate_cluster_numbers",
    "region_number",
    "kind",
]
ANTISMASH_INFO = [
    "category",
    "contig_edge",
    "product",
    "protocluster_number",
    "candidate_cluster_number",
    "candidate_cluster_numbers",
    "region_number",
    "tool",
    "kind",
]  #'detection_rule',
ANTISMASH_TYPES = ("region", "cand_cluster", "protocluster", "proto_core")


def out_dict(x: Dict):
    # Creates a uniform, OrderedDict for the output and maps the input diction
    temp = OrderedDict({k: None for k in OUT_COLUMNS[1:]})
    for k, v in x.items():
        temp[k] = v
    return temp


def unlist(x: Union[List, str]):
    # Some antismash fields are lists, for simplicity they are collapsed into a single string
    if isinstance(x, list):
        if len(x) > 1:
            to_return = "; ".join(x)
        else:
            to_return = x[0]
    if isinstance(to_return, list):
        raise ValueError
    return to_return


def get_seqio_start(seq_feature: SeqFeature, offset: int):
    # otherwise starts are zero-indexed which is fine if you know, but not the same as in the genbank file
    return seq_feature.location.start.real + 1 + offset


def get_seqio_end(seq_feature: SeqFeature, offset: int):
    # otherwise starts are zero-indexed which is fine if you know, but not the same as in the genbank file
    return seq_feature.location.end.real + offset


def parse_and_write(gbk_path_list, outpath, append=False, header=False):
    if append:
        write_type = "a"
    else:
        write_type = "w"
    with open(outpath, write_type) as handle:
        tsv_writer = csv.writer(handle, delimiter="\t")
        if header:
            tsv_writer.writerow(OUT_COLUMNS)
        counter = 1
        input_len = len(gbk_path_list)
        for gbk_path in gbk_path_list:
            print(f"{counter}/{input_len}", end="\r")
            counter += 1
            parse_gen = SeqIO.parse(gbk_path, "genbank")
            for seq_record in parse_gen:
                offset = int(
                    seq_record.annotations["structured_comment"]["antiSMASH-Data"][
                        "Orig. start"
                    ]
                )

                extracted_features = (
                    i for i in seq_record.features if i.type in ANTISMASH_TYPES
                )
                extracted_antismash_feature_info = (
                    {
                        "type": i.type,
                        "start": get_seqio_start(i, offset),
                        "end": get_seqio_end(i, offset),
                        "strand": i.location.strand,
                    }
                    | {
                        k: unlist(v)
                        for k, v in i.qualifiers.items()
                        if k in ANTISMASH_INFO
                    }
                    for i in extracted_features
                )
                results = (out_dict(i) for i in extracted_antismash_feature_info)
                for i in results:
                    tsv_writer.writerow([seq_record.id] + [v for v in i.values()])
