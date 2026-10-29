"""Leitura dos relatórios DArTseq -> marcadores com coordenadas e genótipo por cultivar.

Suporta o formato 1-row 'SNP_mapping' (escores 0/1/2/-) e 'SilicoDArT' (presença 1/0/-).
As triplicatas (colunas com o mesmo rótulo 'genotype <code>') são colapsadas em
um consenso por cultivar, com flag de inconsistência entre réplicas.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path


def _collapse(labels, vals):
    groups = defaultdict(list)
    for lab, v in zip(labels, vals):
        code = lab.replace('"', '').replace('genotype ', '').strip()
        groups[code].append(v)
    call, incons = {}, {}
    for code, vs in groups.items():
        obs = [x for x in vs if x not in ('-', '', 'NN', '--')]
        if not obs:
            call[code] = '-'; incons[code] = False; continue
        c = Counter(obs); top, topn = c.most_common(1)[0]
        call[code] = top
        incons[code] = (len(c) > 1 and topn < len(obs))
    return call, incons


def _parse(fn, token, seqc, chromc, posc, crc, repc, picc, sstart, src):
    markers = []
    with open(fn) as fh:
        r = csv.reader(fh); header = None
        for row in r:
            if row and row[0] == token:
                header = row; break
        if header is None:
            raise ValueError(f"Cabeçalho '{token}' não encontrado em {fn}")
        labels = header[sstart:]
        for row in r:
            if not row or not row[0] or row[0] == '*':
                continue
            chrom = row[chromc].strip()
            if not chrom:
                continue  # só marcadores alinhados ao genoma
            try:
                pos = int(float(row[posc]))
            except ValueError:
                continue
            try:
                cr = float(row[crc]); rep = float(row[repc]); pic = float(row[picc])
            except (ValueError, IndexError):
                cr = rep = pic = 0.0
            seq = row[seqc].strip().upper()
            if len(seq) < 30:
                continue
            call, incons = _collapse(labels, row[sstart:])
            markers.append(dict(id=row[0], chrom=chrom, pos=pos, seq=seq,
                                cr=cr, rep=rep, pic=pic, call=call, incons=incons, src=src))
    return labels, markers


def load_reports(reports_dir: Path):
    reports_dir = Path(reports_dir)
    snp_files = sorted(reports_dir.glob("*SNP_mapping*.csv"))
    sil_files = sorted(reports_dir.glob("*SilicoDArT*.csv"))
    if not snp_files or not sil_files:
        raise FileNotFoundError(
            f"Esperado um *SNP_mapping*.csv e um *SilicoDArT*.csv em {reports_dir}")
    # índices de coluna conforme o cabeçalho padrão DArT
    _, snp = _parse(snp_files[0], 'AlleleID', 2, 4, 5, 10, 21, 18, 22, 'SNP')
    _, sil = _parse(sil_files[0], 'CloneID', 1, 3, 4, 7, 13, 9, 14, 'SIL')
    cult = sorted({c for m in snp for c in m['call']})
    return dict(snp=snp, sil=sil, cult=cult)
