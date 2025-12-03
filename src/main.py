import argparse
import subprocess

'''
    Usage: main.py [OPTIONS]
    
    Options:
      --annotator <local|external|both>  Choose which annotator to run (required).
      --output <str>                     Path to output file (required).
      --somatic_variants <str>           Path to somatic variants file.
      --copy_number_alterations <str>    Path to copy number alterations file.
      --ascatestimates <str>             Path to ASCAT estimates file (required).
      --cn_annotations <str>             Path to filtered and annotated CNAs.
      --tumortype <str>                  Tumor type identifier (default: HGSOC).
      --refgen <str>                     Reference genome version (default: GRCh38).
      --cores <int>                      Number of cores to use for processing (default: 0).
      --post_annotate                    Perform post-annotation.
      --ploidy_threshold <float>         Threshold for filtering by ploidy (default: 2.5).
      --pid <str>                        Patient ID.
      --homogeneity_threshold <float>    Homogeneity threshold for variant filtering (default: 0.05).
      --rf_score_threshold <float>       Random Forest score threshold for variant filtering (default: 0.95).
      --ada_score_threshold <float>      AdaBoost score threshold for variant filtering (default: 0.95).
    
    Examples:
      python main.py --annotator local --output path/to/output --somatic_variants path/to/somatic_mutations.tsv --ascatestimates path/to/ascat.tsv
      python main.py --annotator local --output path/to/output --copy_number_alterations path/to/cnas.tsv --ascatestimates path/to/ascat.tsv
      
      External annotator examples. Appends annotation columns to existing file:
      python main.py --annotator external --copy_number_alterations path/to/cnas.tsv --ascatestimates path/to/ascat.tsv
      python main.py --annotator external --oncokbcna --copy_number_alterations path/to/cnas.tsv
      python main.py --annotator external --oncokbsomatic_mutation --somatic_variants path/to/somatic_mutations.tsv
      python main.py --annotator external --cgiquery --somatic_variants path/to/somatic_mutations.tsv
      python main.py --annotator external --cgiquery --copy_number_alterations path/to/cnas.tsv
      '''

def run_local_annotator(args):
    command = ["python", "local_annotator.py"]
    for key, value in vars(args).items():
        if value is not None:
            if isinstance(value, bool) and value:
                command.append(f"--{key}")
            else:
                command.extend([f"--{key}", str(value)])
    subprocess.run(command)


def run_external_annotator(args):
    command = ["python", "external_annotator.py"]
    for key, value in vars(args).items():
        if value is not None:
            if isinstance(value, bool) and value:
                command.append(f"--{key}")
            else:
                command.extend([f"--{key}", str(value)])
    subprocess.run(command)


def main():
    parser = argparse.ArgumentParser(description="Run local and external annotators")
    parser.add_argument("--annotator", choices=["local", "external", "both"], required=True,
                        help="Choose which annotator to run")

    # Common arguments
    parser.add_argument("--output", type=str, help="Path to output file", required=True)
    parser.add_argument("--somatic_variants", type=str, help="Path to somatic variants file")
    parser.add_argument("--copy_number_alterations", type=str, help="Path to copy number alterations file")
    parser.add_argument("--ascatestimates", type=str, help="Path to ASCAT estimates file", required=True)
    parser.add_argument("--cn_annotations", type=str, help="Path to filtered and annotated CNAs")
    parser.add_argument("--tumortype", type=str, help="Tumor type identifier", default="HGSOC")
    parser.add_argument("--refgen", type=str, help="Reference genome version", default="GRCh38")
    parser.add_argument("--cores", type=int, help="Number of cores to use for processing", default=0)
    parser.add_argument("--post_annotate", action="store_true", help="Perform post-annotation")
    parser.add_argument("--ploidy_threshold", type=float, help="Threshold for filtering by ploidy", default=2.5)
    parser.add_argument("--pid", type=str, help="Patient ID")
    parser.add_argument("--homogeneity_threshold", type=float, help="Homogeneity threshold for variant filtering",
                        default=0.05)
    parser.add_argument("--rf_score_threshold", type=float, help="Random Forest score threshold for variant filtering",
                        default=0.95)
    parser.add_argument("--ada_score_threshold", type=float, help="AdaBoost score threshold for variant filtering",
                        default=0.95)
    # External annotator specific arguments
    parser.add_argument("--oncokbcna", action="store_true", help="Query OncoKB for copy number alterations")
    parser.add_argument("--oncokbsomatic_mutation", action="store_true", help="Query OncoKB for somatic mutations")
    parser.add_argument("--cgiquery", action="store_true", help="Query Cancer Genome Interpreter")
    parser.add_argument("--cgijobid", type=str, help="Download results from CGI by jobid")

    args = parser.parse_args()

    if args.annotator in ["local", "both"]:
        run_local_annotator(args)

    if args.annotator in ["external", "both"]:
        run_external_annotator(args)


if __name__ == "__main__":
    main()