import multiprocessing
import os
import re
from typing import List

import GEOparse
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from tqdm import tqdm

from utils import genbank_to_gene_symbol


class GeneExpressionDataset:

    def __init__(
        self, data: np.array, labels: np.array, gene_symbols: List, label_names: List
    ):
        self.data = data.astype(np.float32)
        self.labels = labels.astype(np.float32)
        self.gene_symbols = gene_symbols
        self.label_names = label_names


def load_GSE49541(destdir="./data", return_deseq=True):

    # Gather raw dataset
    gse = GEOparse.get_GEO(geo="GSE49541", destdir=destdir, silent=True)
    data = []
    metadata = []
    for i, (gsm_name, gsm) in enumerate(gse.gsms.items()):
        metadata.append(gsm.metadata)
        table = gsm.table
        table.columns = ["ID_REF", gsm_name]
        data.append(table if i == 0 else table[gsm_name])
    data = pd.concat(data, axis=1)
    metadata = pd.DataFrame(metadata)

    # Process gene expression data
    probes = list(gse.gpls.values())[0].table
    keep = ~(probes["Gene Symbol"].isna().values)  # Remove rows with missing gene names
    samples = data.columns[1:]
    data = data.loc[keep].T[1:]
    data.insert(0, "sample", samples)
    data = data.set_index("sample")
    symbols = probes["Gene Symbol"].values[keep]
    symbols = [s.split(" /// ") for s in symbols]

    # Process labels
    stages = [
        re.search(r"Stage:\s*(\w+)", m[0]).group(1)
        for m in metadata["characteristics_ch1"]
    ]
    mapping = {"mild": "NAFL", "advanced": "NASH"}
    stages = [mapping[stage] for stage in stages]
    labels = pd.DataFrame({"sample": samples, "condition": stages})
    labels = labels.set_index("sample")

    if return_deseq:
        data = (
            data.astype(float).multiply(1000).astype(int)
        )  # Convert to integer (~counts)
        dataset = DeseqDataSet(
            counts=data, metadata=labels, design_factors="condition", quiet=True
        )
        return dataset, symbols
    else:
        labels = pd.Categorical(stages).codes
        label_names = np.unique(stages).tolist()
        data = data.values
        dataset = GeneExpressionDataset(data, labels, symbols, label_names)
        return dataset


def load_GSE48452(destdir="./data", return_deseq=True):

    # Gather raw dataset
    gse = GEOparse.get_GEO(geo="GSE48452", destdir=destdir, silent=True)
    data = []
    metadata = []
    for i, (gsm_name, gsm) in enumerate(gse.gsms.items()):
        metadata.append(gsm.metadata)
        table = gsm.table
        table.columns = ["ID_REF", gsm_name]
        data.append(table if i == 0 else table[gsm_name])
    data = pd.concat(data, axis=1)
    metadata = pd.DataFrame(metadata)

    # Process gene expression data
    probes = list(gse.gpls.values())[0].table
    keep = ~(
        probes["GB_LIST"].isna().values | (probes["gene_assignment"] == "---").values
    )  # Remove rows with missing gene names
    samples = data.columns[1:]
    data = data.loc[keep].T[1:]
    data.insert(0, "sample", samples)
    data = data.set_index("sample")

    # Gene symbols
    symbols = []
    for gbs, assignment in probes[["GB_LIST", "gene_assignment"]].values[keep]:
        gbs = gbs.split(",")
        parsed = []
        for gb in gbs:
            try:
                parsed.append(re.search(rf"{gb} // (\w+) //", assignment).group(1))
            except:
                pass
        parsed = list(set(parsed))
        symbols.append(parsed)

    # Process labels
    stages = [
        next((item.split(": ")[1] for item in m if item.startswith("group:")), None)
        for m in metadata["characteristics_ch1"]
    ]
    mapping = {
        "Control": "control",
        "Healthy obese": "healthy obese",
        "Steatosis": "NAFL",
        "Nash": "NASH",
    }
    stages = [mapping[stage] for stage in stages]
    labels = pd.DataFrame({"sample": samples, "condition": stages})
    labels = labels.set_index("sample")

    if return_deseq:
        data = (
            data.astype(float).multiply(1000).astype(int)
        )  # Convert to integer (~counts)
        dataset = DeseqDataSet(
            counts=data, metadata=labels, design_factors="condition", quiet=True
        )
        return dataset, symbols
    else:
        labels = pd.Categorical(stages).codes
        label_names = np.unique(stages).tolist()
        data = data.values
        dataset = GeneExpressionDataset(data, labels, symbols, label_names)
        return dataset


def load_GSE89632(destdir="./data", return_deseq=True):

    # Gather raw dataset
    gse = GEOparse.get_GEO(geo="GSE89632", destdir=destdir, silent=True)
    data = []
    metadata = []
    for i, (gsm_name, gsm) in enumerate(gse.gsms.items()):
        metadata.append(gsm.metadata)
        table = gsm.table[["ID_REF", "VALUE"]]
        table.columns = ["ID_REF", gsm_name]
        data.append(table if i == 0 else table[gsm_name])
    data = pd.concat(data, axis=1)
    metadata = pd.DataFrame(metadata)

    # Process gene expression data
    probes = list(gse.gpls.values())[0].table
    keep = ~(probes["Symbol"].isna().values)  # Remove rows with missing gene names
    samples = data.columns[1:]
    data = data.loc[keep].T[1:]
    data.insert(0, "sample", samples)
    data = data.set_index("sample")
    symbols = probes["Symbol"].values[keep].tolist()

    # Process labels
    stages = [
        next((item.split(": ")[1] for item in m if item.startswith("diagnosis:")), None)
        for m in metadata["characteristics_ch1"]
    ]
    mapping = {"HC": "control", "SS": "NAFL", "NASH": "NASH"}
    stages = [mapping[stage] for stage in stages]
    labels = pd.DataFrame({"sample": samples, "condition": stages})
    labels = labels.set_index("sample")

    if return_deseq:
        data = (
            data.astype(float).multiply(1000).astype(int)
        )  # Convert to integer (~counts)
        dataset = DeseqDataSet(
            counts=data, metadata=labels, design_factors="condition", quiet=True
        )
        return dataset, symbols
    else:
        labels = pd.Categorical(stages).codes
        label_names = np.unique(stages).tolist()
        data = data.values
        dataset = GeneExpressionDataset(data, labels, symbols, label_names)
        return dataset


def load_GSE83452(destdir="./data", return_deseq=True):

    # Gather raw dataset
    gse = GEOparse.get_GEO(geo="GSE83452", destdir=destdir, silent=True)
    data = []
    metadata = []
    for i, (gsm_name, gsm) in enumerate(gse.gsms.items()):
        time = next(
            (
                item.split(": ")[1]
                for item in gsm.metadata["characteristics_ch1"]
                if item.startswith("time:")
            ),
            None,
        )
        if time != "baseline":
            continue
        metadata.append(gsm.metadata)
        table = gsm.table[["ID_REF", "VALUE"]]
        table.columns = ["ID_REF", gsm_name]
        data.append(table if i == 0 else table[gsm_name])
    data = pd.concat(data, axis=1)
    metadata = pd.DataFrame(metadata)

    # Gather gene symbols
    probes = list(gse.gpls.values())[0].table
    f = f"{destdir}/GSE83452_symbols.npy"
    if os.path.exists(f):
        symbols = np.load(f)
    else:

        def convert(task):
            index = task[0]
            value = task[1]
            if not isinstance(value, str) and np.isnan(value):
                return index, ""
            gene_symbol = genbank_to_gene_symbol(value)
            if gene_symbol is None:
                return index, ""
            return index, gene_symbol

        tasks = probes["GB_ACC"].values
        tasks = list(zip(range(len(tasks)), tasks))
        pool = multiprocessing.Pool(processes=8)
        symbols = np.empty(len(tasks), dtype=object)
        for index, gene_symbol in tqdm(
            pool.imap_unordered(convert, tasks), total=len(tasks)
        ):
            symbols[index] = gene_symbol
        symbols = symbols.astype(str)
        np.save(f, symbols)

    # Process gene expression data
    symbols = symbols[probes["ID"].isin(data["ID_REF"].values.astype(str))]
    keep = symbols != ""  # Remove rows with missing gene names
    samples = data.columns[1:]
    data = data.loc[keep].T[1:]
    data.insert(0, "sample", samples)
    data = data.set_index("sample")
    symbols = symbols[keep].tolist()

    # Process labels
    stages = [
        next(
            (item.split(": ")[1] for item in m if item.startswith("liver status:")),
            None,
        )
        for m in metadata["characteristics_ch1"]
    ]
    mapping = {"undefined": "other", "no NASH": "control", "NASH": "NASH"}
    stages = [mapping[stage] for stage in stages]
    labels = pd.DataFrame({"sample": samples, "condition": stages})
    labels = labels.set_index("sample")

    if return_deseq:
        data = (
            data.astype(float).multiply(1000).astype(int)
        )  # Convert to integer (~counts)
        dataset = DeseqDataSet(
            counts=data, metadata=labels, design_factors="condition", quiet=True
        )
        return dataset, symbols
    else:
        labels = pd.Categorical(stages).codes
        label_names = np.unique(stages).tolist()
        data = data.values
        dataset = GeneExpressionDataset(data, labels, symbols, label_names)
        return dataset
