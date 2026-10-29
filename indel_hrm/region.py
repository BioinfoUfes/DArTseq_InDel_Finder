"""Build the primer-design template around each InDel.

Ways of obtaining the genomic context:
  - remap : locate the tag in the genome BY SEQUENCE (an anchor of remap_seed bp next to
            the InDel), ignoring ChromPos. Works with ANY genome/assembly. It also counts
            how many places the anchor matches (specificity): > remap_max_hits => the
            marker is non-specific.
  - coord : use the report's Chrom_/ChromPos_ (only valid on the assembly the data were
            aligned to) and extract a +-flank bp window, anchoring the InDel inside it.
  - tag   : no genome available — use the tag itself (~69 bp) as the template (fallback).

Returns, per marker: the template (the "long" allele, carrying the InDel), the InDel
position and length in the template, the "short" allele (without the InDel), the template
source and the number of genome hits.
"""
from difflib import SequenceMatcher


def _locate_indel(seqA, seqB):
    """Align the two tags; return (long, short, indel_start_in_long, indel_len)."""
    sm = SequenceMatcher(None, seqA, seqB, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == 'delete' and (i2 - i1) >= 1:
            return seqA, seqB, i1, i2 - i1
        if tag == 'insert' and (j2 - j1) >= 1:
            return seqB, seqA, j1, j2 - j1
    return seqA, seqB, 0, 0


def _open_maybe_gz(path):
    import gzip
    path = str(path)
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _load_fasta_index(fna_paths, hap_tags=('SWO-M', 'SWO-P', 'M', 'P')):
    """Accept one path or several (e.g. the two haplotypes of a phased assembly) and
    .fa/.fa.gz. Each sequence gets a haplotype label (derived from the file name) so the
    specificity check can tell homeologs (1x in each haplotype) from paralogs (2x in the
    same haplotype)."""
    from Bio import SeqIO
    if isinstance(fna_paths, (str, bytes)) or hasattr(fna_paths, '__fspath__'):
        fna_paths = [fna_paths]
    idx = {}            # id -> (seq, haplotype)
    for fp in fna_paths:
        name = str(fp)
        hap = next((h for h in hap_tags if h in name), name.split('/')[-1])
        with _open_maybe_gz(fp) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                base = rec.id.split()[0]
                seq = str(rec.seq).upper()
                idx[base] = (seq, hap)
                idx.setdefault(base.split('.')[0], (seq, hap))
    return idx


def _count_hits(genome_index, seed, cap=6):
    """Count anchor occurrences (+ strand) PER HAPLOTYPE, up to 'cap' total.
    Returns (hits_per_hap: dict, first_hit: (chrom, pos, hap))."""
    from collections import defaultdict
    per_hap = defaultdict(int)
    first = None
    total = 0
    seen = set()
    for cid, (seq, hap) in genome_index.items():
        if id(seq) in seen:            # skip the version-less alias
            continue
        seen.add(id(seq))
        start = 0
        while True:
            p = seq.find(seed, start)
            if p < 0:
                break
            per_hap[hap] += 1; total += 1
            if first is None:
                first = (cid, p, hap)
            if total >= cap:
                return dict(per_hap), first
            start = p + 1
    return dict(per_hap), first


def _anchor(genome_index, seq_long, istart, ilen, cfg):
    """Anchor the InDel in the genome using seeds from BOTH flanks (SNP-tolerant: a single
    SNP-free seed is enough; covers InDels near either tag end).
    Returns ((chrom, astart), per_hap, paralog_info); astart = InDel position in the genome."""
    K = min(cfg.remap_seed, 20)
    paralog = None
    cands = []
    # left flank: seed ends at 'end' (<= istart); astart = gpos + len + (istart - end)
    for end in (istart, istart - 6, istart - 12, istart - 20, istart - 30, istart - 42):
        if end >= 14:
            cands.append(('L', max(0, end - K), end))
    # right flank: seed starts at 'start' (>= istart+ilen); astart = gpos - (start - (istart+ilen))
    rbase = istart + ilen
    for start in (rbase, rbase + 6, rbase + 12, rbase + 20, rbase + 30):
        if start + 14 <= len(seq_long):
            cands.append(('R', start, min(len(seq_long), start + K)))
    for side, s, e in cands:
        seed = seq_long[s:e]
        if len(seed) < 14:
            continue
        per_hap, hit = _count_hits(genome_index, seed)
        if not hit:
            continue
        if max(per_hap.values()) <= cfg.remap_max_hits:
            cid, gpos, _ = hit
            if side == 'L':
                astart = gpos + len(seed) + (istart - e)
            else:
                astart = gpos - (s - rbase)
            if 0 <= astart <= len(genome_index[cid][0]):
                return (cid, astart), per_hap, None
        else:
            paralog = per_hap
    return None, None, paralog


def build_regions(selected, cfg, genome_index=None):
    regions = []
    for o in selected:
        seq_long, seq_short, istart, ilen = _locate_indel(o['seqA'], o['seqB'])
        source, template, tstart, tlen, nhits = 'tag', seq_long, istart, ilen, None

        if genome_index is not None and ilen > 0:
            anchor = None
            if cfg.remap:
                anchor, per_hap, paralog = _anchor(genome_index, seq_long, istart, ilen, cfg)
                if anchor is None and paralog is not None:
                    source = f'tag(paralog: {paralog})'
                    nhits = sum(paralog.values())
            else:
                chrom_base = o['chrom'].split('_chr')[0]
                entry = genome_index.get(chrom_base) or genome_index.get(chrom_base.split('.')[0])
                if entry:
                    seed = seq_long[max(0, istart - cfg.remap_seed):istart]
                    lo = max(0, o['pos'] - cfg.flank); hi = min(len(entry[0]), o['pos'] + cfg.flank)
                    p = entry[0][lo:hi].find(seed) if seed else -1
                    if p >= 0:
                        anchor = (chrom_base, lo + p + len(seed)); per_hap = {entry[1]: 1}

            if anchor:
                cid, astart = anchor
                genome, hap = genome_index[cid]
                lo = max(0, astart - cfg.flank); hi = min(len(genome), astart + cfg.flank)
                left = genome[lo:astart]; right = genome[astart:hi]
                template = left + seq_long[istart:istart + ilen] + right
                tstart = len(left); tlen = ilen
                source = f'genome:{hap}'; nhits = per_hap

        alt = template[:tstart] + template[tstart + tlen:]
        regions.append(dict(**o, template=template, indel_start=tstart, indel_len=tlen,
                            alt_template=alt, region_source=source, genome_hits=nhits))
    return regions
