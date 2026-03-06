import argparse
from pathlib import Path

import pandas as pd


def map_gistic(row):
    if row["nMajor"] == 0 and row["nMinor"] == 0:
        return -2
    if row["nMajor"] == 0 or row["nMinor"] == 0:
        return -1
    if "AMPLIFICATION" in row["alteration"]:
        return 2
    return pd.NA


def resolve_calls(values):
    vals = set(pd.Series(values).dropna().astype(int).tolist())
    if 2 in vals:
        return 2
    if -2 in vals:
        return -2
    if -1 in vals:
        return -1
    return pd.NA


def build_matrix(input_path: str, output_path: str) -> pd.DataFrame:
    cna = pd.read_csv(input_path, sep="\t")

    required_cols = {"hugoSymbol", "sample_id", "alteration", "nMajor", "nMinor"}
    missing_cols = required_cols - set(cna.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {sorted(missing_cols)}")

    cna["nMajor"] = pd.to_numeric(cna["nMajor"], errors="coerce")
    cna["nMinor"] = pd.to_numeric(cna["nMinor"], errors="coerce")
    cna["alteration"] = cna["alteration"].astype(str).str.upper()

    cna["gistic"] = cna.apply(map_gistic, axis=1)
    cna = cna.dropna(subset=["gistic"])

    matrix = (
        cna.groupby(["hugoSymbol", "sample_id"], as_index=False)["gistic"]
        .agg(resolve_calls)
        .pivot(index="hugoSymbol", columns="sample_id", values="gistic")
        .reset_index()
        .rename(columns={"hugoSymbol": "Gene Symbol"})
    )

    sample_cols = [c for c in matrix.columns if c != "Gene Symbol"]
    for col in sample_cols:
        matrix[col] = matrix[col].astype("Int64")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_path, sep="\t", index=False)
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert CNA table to OncoKB annotator matrix format."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input TSV with columns: hugoSymbol, sample_id, alteration, nMajor, nMinor",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output TSV path for OncoKB annotator CNA matrix",
    )
    args = parser.parse_args()

    matrix = build_matrix(args.input, args.output)
    print("Output:", args.output)
    print("Shape:", matrix.shape)


if __name__ == "__main__":
    main()