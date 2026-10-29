"""Command-line interface: interactive questions and non-interactive argument parsing."""
import argparse
import sys
from pathlib import Path

from .config import Config
from .io_report import load_reports


def _ask(prompt, default=None, cast=str):
    sfx = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{prompt}{sfx}: ").strip()
        if not raw and default is not None:
            return default
        if not raw:
            print("  (answer required)")
            continue
        try:
            return cast(raw)
        except Exception:
            print("  (invalid value, try again)")


def _ask_path(prompt, default=None, must_exist=True):
    while True:
        p = _ask(prompt, default)
        path = Path(p).expanduser()
        if not must_exist or path.exists():
            return path
        print(f"  (path not found: {path})")


def interactive_config():
    print("=" * 64)
    print(" DArTseq_InDel_Finder — interactive configuration")
    print("=" * 64)
    reports = _ask_path("Folder with the DArTseq reports (SNP_mapping + SilicoDArT)",
                        "reports")
    # genomes
    print("\nReference genome(s) for your species (FASTA .fa/.fa.gz).")
    print("You may give several separated by commas (e.g. two haplotypes). Enter = none.")
    graw = input("Genome file(s): ").strip()
    genome = [Path(g.strip()).expanduser() for g in graw.split(",") if g.strip()] or None
    if genome:
        for g in genome:
            if not g.exists():
                print(f"  WARNING: could not find {g}")

    meta = input("\nCSV of sample names (columns code,name) — Enter to skip: ").strip()
    metadata = Path(meta).expanduser() if meta else None

    mode = _ask("\nMode: 'tree' (tree for one group) or 'panel' (panel for all)", "tree")
    while mode not in ("tree", "panel"):
        mode = _ask("Type 'tree' or 'panel'", "tree")

    # available samples
    try:
        data = load_reports(reports)
        codes = data["cult"]
        print(f"\nSamples in report ({len(codes)}): {', '.join(codes)}")
    except Exception as e:
        print(f"  (could not list samples now: {e})")
        codes = []

    samples = None
    if mode == "tree":
        sraw = input("Sample codes for the group (comma-separated) or 'all': ").strip()
        if sraw and sraw.lower() != "all":
            samples = [s.strip() for s in sraw.split(",") if s.strip()]

    print("\nConditions of your HRM kit/master mix (Enter = default):")
    mg = _ask("  [Mg2+] mM", 3.0, float)
    mono = _ask("  [Mono+] mM", 20.0, float)

    indel_max = _ask("\nMaximum indel size to search (bp)", 45, int)
    flank = _ask("Flank extracted from the genome for primer design (bp each side)", 180, int)
    out = _ask("\nOutput folder", "indel_hrm_output")

    cfg = Config(reports_dir=reports, genome_fna=genome, out_dir=Path(out),
                 metadata_csv=metadata, indel_max=indel_max, flank=flank,
                 melt_mg=mg, melt_mono=mono)
    return cfg, samples, mode


def parse_args(argv):
    ap = argparse.ArgumentParser(description="DArTseq_InDel_Finder")
    ap.add_argument("--reports", help="folder with the DArTseq reports")
    ap.add_argument("--genome", nargs="*", help="reference genome FASTA(s) (.fa/.fa.gz)")
    ap.add_argument("--metadata", help="CSV with columns code,name to label the samples")
    ap.add_argument("--samples", help="sample codes (comma-separated) for 'tree' mode; omit = all")
    ap.add_argument("--mode", choices=["tree", "panel"], default="tree")
    ap.add_argument("--out", default="indel_hrm_output", help="output folder")
    ap.add_argument("--melt-mg", type=float, default=3.0, help="[Mg2+] mM of the HRM master mix")
    ap.add_argument("--melt-mono", type=float, default=20.0, help="[Mono+] mM of the HRM master mix")
    ap.add_argument("--indel-max", type=int, default=45, help="maximum indel size (bp)")
    ap.add_argument("--flank", type=int, default=180, help="genomic flank for primer (bp)")
    ap.add_argument("--n-trees", type=int, default=3, help="number of alternative trees")
    ap.add_argument("--min-exact-anchor", type=int, default=0,
                    help="minimum exact run flanking the gap, each side (0 = off)")
    ap.add_argument("--insilico-filter", action="store_true",
                    help="drop markers whose primer pair is not in-silico specific")
    ap.add_argument("-i", "--interactive", action="store_true", help="force interactive mode")
    return ap.parse_args(argv)


def build_cfg(a):
    genome = [Path(g) for g in a.genome] if a.genome else None
    cfg = Config(reports_dir=Path(a.reports), genome_fna=genome, out_dir=Path(a.out),
                 metadata_csv=Path(a.metadata) if a.metadata else None,
                 indel_max=a.indel_max, flank=a.flank,
                 melt_mg=a.melt_mg, melt_mono=a.melt_mono,
                 min_exact_anchor=a.min_exact_anchor, insilico_filter=a.insilico_filter)
    samples = None
    if a.samples and a.samples.lower() != "all":
        samples = [s.strip() for s in a.samples.split(",") if s.strip()]
    return cfg, samples, a.mode, a.n_trees


def main_cli(argv=None):
    from .pipeline import run as run_pipeline
    from .report_hrm import write_report, load_names

    argv = sys.argv[1:] if argv is None else argv
    a = parse_args(argv)
    if a.interactive or not a.reports:
        cfg, samples, mode = interactive_config()
        n_trees = a.n_trees
    else:
        cfg, samples, mode, n_trees = build_cfg(a)

    result = run_pipeline(cfg, samples=samples, mode=mode, n_trees=n_trees)
    names = load_names(cfg.metadata_csv)
    outputs = write_report(result, cfg, cfg.out_dir, names=names)

    print("\n" + "=" * 64)
    print("DONE. Outputs in:", cfg.out_dir.resolve())
    print("  workbook :", outputs["xlsx"].name)
    print("  tables   : candidates.tsv, primers.tsv, per_pair.tsv")
    print("  summary  :", outputs["summary_json"].name, "(parameters + chosen tree, machine-readable)")
    print("  amplicons:", outputs["fasta"].name, "(for checking in uMelt)")
    print("  curves   :", outputs["images"].name + "/ (per-marker images)")
    trees = result.get("trees", [])
    if trees:
        t = trees[0]
        print(f"\nRecommended tree: root {t['node']['chrom'].split('_chr')[-1]}:{t['node']['pos']} "
              f"| {'STRONG' if t['forte'] else 'uses N'} | weakest link {t['wmin']}% "
              f"| up to {t['depth']} analyses")
    elif mode == "tree":
        print("\nNo tree built — check whether the reference genome was provided.")
    return result, outputs
