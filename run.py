#!/usr/bin/env python3
"""DArTseq_InDel_Finder — entry point.

Interactive use (asks the questions):
    python run.py

Command-line use (non-interactive):
    python run.py --reports reports --genome genomeA.fa.gz genomeB.fa.gz \
                  --samples 7011,7041,7115,7167 --mode tree --out output \
                  --metadata samples.csv --melt-mg 3 --melt-mono 20

Installed as a package, the same CLI is available as `dartseq-indel-finder`
or `python -m indel_hrm`.
"""
from indel_hrm.cli import main_cli

if __name__ == "__main__":
    main_cli()
