import argparse
from pathlib import Path

import pandas as pd


def convert_snvs(mapping_path: str, local_path: str, output_path: str) -> pd.DataFrame:
    mapping = pd.read_csv(mapping_path)
    local = pd.read_csv(local_path, sep="\t")

    # Ensure mapping has expected columns
    if not {"local", "oncokb"}.issubset(mapping.columns):
        raise ValueError("Mapping CSV must have 'local' and 'oncokb' columns")

    # Initialize oncokb with all mapped OncoKB columns
    oncokb_columns = mapping["oncokb"].dropna().unique().tolist()
    oncokb = pd.DataFrame(index=local.index, columns=oncokb_columns)

    # Fill oncokb columns from local according to mapping
    for _, row in mapping.iterrows():
        local_col = row["local"]
        oncokb_col = row["oncokb"]
        if local_col in local.columns:
            oncokb[oncokb_col] = local[local_col]
        else:
            # Column missing in local input; keep as NA and later convert to empty string.
            oncokb[oncokb_col] = pd.NA

    # Normalize Chromosome by removing optional chr prefix.
    if "Chromosome" in oncokb.columns:
        chrom = oncokb["Chromosome"].astype(str).str.strip()
        chrom = chrom.str.replace(r"(?i)^chr", "", regex=True)
        oncokb["Chromosome"] = chrom.str.upper()

    # Optional transformation: derive HGVSp_Short from AAChangeMANE if available.
    if "AAChangeMANE" in local.columns:
        oncokb["HGVSp_Short"] = local["AAChangeMANE"].astype(str).str.split(".").str[-1]

    if "HGVSp_Short" in oncokb.columns:
        oncokb["HGVSp_Short"] = oncokb["HGVSp_Short"].fillna("")

    # Derive End_Position from local: position + (len(reference_allele) - len(sample_allele)).
    if all(col in local.columns for col in ["position", "reference_allele", "sample_allele"]):
        ref_len = local["reference_allele"].fillna("").astype(str).str.len()
        alt_len = local["sample_allele"].fillna("").astype(str).str.len()
        pos = pd.to_numeric(local["position"], errors="coerce")
        oncokb["End_Position"] = pos + (ref_len - alt_len)

    oncokb = oncokb.fillna("")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    oncokb.to_csv(output_path, sep="\t", index=False)
    return oncokb


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert local SNV table to OncoKB annotator SNV format."
    )
    parser.add_argument(
        "--mapping",
        default="src/local_to_oncokb_mapping.csv",
        help="Mapping CSV with 'local' and 'oncokb' columns",
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input local SNV TSV path",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output TSV path for OncoKB annotator",
    )
    args = parser.parse_args()

    oncokb = convert_snvs(args.mapping, args.input, args.output)
    print("Output:", args.output)
    print("oncokb columns:", sorted(oncokb.columns.tolist())[:12], "...")


if __name__ == "__main__":
    main()
