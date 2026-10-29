#!/usr/bin/env python3
"""indel-hrm-finder — ponto de entrada.

Uso interativo (faz as perguntas):
    python run.py

Uso por linha de comando (não interativo):
    python run.py --reports reports --genome genomaA.fa.gz genomaB.fa.gz \
                  --samples 7011,7041,7115,7167 --mode tree --out saida \
                  --metadata amostras.csv --melt-mg 3 --melt-mono 20
"""
import argparse
import sys
from pathlib import Path

from indel_hrm.config import Config
from indel_hrm.pipeline import run as run_pipeline
from indel_hrm.report_hrm import write_report, load_names


def parse_args(argv):
    ap = argparse.ArgumentParser(description="indel-hrm-finder")
    ap.add_argument("--reports", help="pasta com os relatórios DArTseq")
    ap.add_argument("--genome", nargs="*", help="FASTA(s) do genoma de referência (.fa/.fa.gz)")
    ap.add_argument("--metadata", help="CSV com colunas code,name para rotular as amostras")
    ap.add_argument("--samples", help="códigos das amostras (vírgula) para o modo 'tree'; omitir = todas")
    ap.add_argument("--mode", choices=["tree", "panel"], default="tree")
    ap.add_argument("--out", default="indel_hrm_output", help="pasta de saída")
    ap.add_argument("--melt-mg", type=float, default=3.0, help="[Mg2+] mM do master mix de HRM")
    ap.add_argument("--melt-mono", type=float, default=20.0, help="[Mono+] mM do master mix de HRM")
    ap.add_argument("--indel-max", type=int, default=45, help="tamanho máximo de indel (pb)")
    ap.add_argument("--flank", type=int, default=180, help="flanco genômico p/ primer (pb)")
    ap.add_argument("--n-trees", type=int, default=3, help="nº de árvores alternativas")
    ap.add_argument("-i", "--interactive", action="store_true", help="forçar modo interativo")
    return ap.parse_args(argv)


def build_cfg(a):
    genome = [Path(g) for g in a.genome] if a.genome else None
    cfg = Config(reports_dir=Path(a.reports), genome_fna=genome, out_dir=Path(a.out),
                 metadata_csv=Path(a.metadata) if a.metadata else None,
                 indel_max=a.indel_max, flank=a.flank,
                 melt_mg=a.melt_mg, melt_mono=a.melt_mono)
    samples = None
    if a.samples and a.samples.lower() != "all":
        samples = [s.strip() for s in a.samples.split(",") if s.strip()]
    return cfg, samples, a.mode, a.n_trees


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    a = parse_args(argv)
    if a.interactive or not a.reports:
        from indel_hrm.cli import interactive_config
        cfg, samples, mode = interactive_config()
        n_trees = a.n_trees
    else:
        cfg, samples, mode, n_trees = build_cfg(a)

    result = run_pipeline(cfg, samples=samples, mode=mode, n_trees=n_trees)
    names = load_names(cfg.metadata_csv)
    outputs = write_report(result, cfg, cfg.out_dir, names=names)

    print("\n" + "=" * 64)
    print("CONCLUÍDO. Saídas em:", cfg.out_dir.resolve())
    print("  planilha :", outputs["xlsx"].name)
    print("  amplicons:", outputs["fasta"].name, "(para conferência no uMelt)")
    print("  curvas   :", outputs["images"].name + "/ (imagens por marcador)")
    trees = result.get("trees", [])
    if trees:
        t = trees[0]
        print(f"\nÁrvore recomendada: raiz {t['node']['chrom'].split('_chr')[-1]}:{t['node']['pos']} "
              f"| {'FORTE' if t['forte'] else 'usa N'} | elo mais fraco {t['wmin']}% | até {t['depth']} análises")
    elif mode == "tree":
        print("\nNenhuma árvore montada — verifique se o genoma de referência foi fornecido.")


if __name__ == "__main__":
    main()
