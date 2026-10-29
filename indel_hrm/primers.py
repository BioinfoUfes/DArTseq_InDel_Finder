"""HRM primer design with primer3-py.

Forces the pair to flank the InDel (SEQUENCE_TARGET) and restricts the amplicon to the
short range typical of HRM. Produces both amplicons (with and without the InDel) for the
melting prediction.
"""
import primer3


def _design(seq_args, global_args):
    # compatible with new and old primer3-py
    if hasattr(primer3, "design_primers"):
        return primer3.design_primers(seq_args, global_args)
    return primer3.bindings.designPrimers(seq_args, global_args)


def design_for_region(region, cfg):
    template = region['template']
    istart, ilen = region['indel_start'], region['indel_len']
    # target = the InDel (if ilen=0, target the midpoint so primers still flank it)
    if ilen > 0:
        target = [istart, ilen]
    else:
        target = [max(1, len(template) // 2), 1]

    seq_args = {
        'SEQUENCE_ID': f"{region['chrom']}_{region['pos']}",
        'SEQUENCE_TEMPLATE': template,
        'SEQUENCE_TARGET': target,
    }
    global_args = {
        'PRIMER_NUM_RETURN': cfg.n_primer_pairs,
        'PRIMER_OPT_SIZE': cfg.primer_opt_size,
        'PRIMER_MIN_SIZE': cfg.primer_min_size,
        'PRIMER_MAX_SIZE': cfg.primer_max_size,
        'PRIMER_OPT_TM': cfg.primer_opt_tm,
        'PRIMER_MIN_TM': cfg.primer_min_tm,
        'PRIMER_MAX_TM': cfg.primer_max_tm,
        'PRIMER_MIN_GC': cfg.primer_min_gc,
        'PRIMER_MAX_GC': cfg.primer_max_gc,
        'PRIMER_PAIR_MAX_DIFF_TM': cfg.primer_max_pair_tm_diff,
        'PRIMER_PRODUCT_SIZE_RANGE': [[cfg.hrm_product_min, cfg.hrm_product_max]],
        'PRIMER_PRODUCT_OPT_SIZE': cfg.hrm_product_opt,
        'PRIMER_SALT_MONOVALENT': cfg.mono_conc,
        'PRIMER_SALT_DIVALENT': cfg.div_conc,
        'PRIMER_DNTP_CONC': cfg.dntp_conc,
    }
    try:
        res = _design(seq_args, global_args)
    except Exception as e:
        return {'ok': False, 'reason': f'primer3 error: {e}', 'pairs': []}

    n = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
    if not n:
        # relaxed retry: wider Tm and product-size ranges
        relaxed = dict(global_args)
        relaxed.update({
            'PRIMER_MIN_TM': cfg.primer_min_tm - 3, 'PRIMER_MAX_TM': cfg.primer_max_tm + 3,
            'PRIMER_PAIR_MAX_DIFF_TM': cfg.primer_max_pair_tm_diff + 1.5,
            'PRIMER_MIN_GC': 20.0, 'PRIMER_MAX_GC': 80.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[max(40, cfg.hrm_product_min - 10), cfg.hrm_product_max + 40]],
            'PRIMER_MAX_SIZE': min(27, cfg.primer_max_size + 2),
        })
        try:
            res = _design(seq_args, relaxed)
            n = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
        except Exception:
            n = 0
    if not n:
        expl = res.get('PRIMER_PAIR_EXPLAIN', '')
        return {'ok': False, 'reason': f'no valid pair ({expl})', 'pairs': []}

    pairs = []
    for i in range(n):
        lstart, llen = res[f'PRIMER_LEFT_{i}']
        rstart, rlen = res[f'PRIMER_RIGHT_{i}']   # rstart = 3' base on the top strand
        prod = res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
        amp_long = template[lstart:rstart + 1]
        # amplicon of the allele without the InDel (drop the InDel bases, if inside the amplicon)
        if lstart <= istart and (istart + ilen) <= rstart + 1 and ilen > 0:
            amp_short = template[lstart:istart] + template[istart + ilen:rstart + 1]
        else:
            amp_short = amp_long  # InDel outside the amplicon (should not happen with TARGET)
        pairs.append(dict(
            idx=i,
            left_seq=res[f'PRIMER_LEFT_{i}_SEQUENCE'],
            right_seq=res[f'PRIMER_RIGHT_{i}_SEQUENCE'],
            left_tm=round(res[f'PRIMER_LEFT_{i}_TM'], 2),
            right_tm=round(res[f'PRIMER_RIGHT_{i}_TM'], 2),
            left_gc=round(res[f'PRIMER_LEFT_{i}_GC_PERCENT'], 1),
            right_gc=round(res[f'PRIMER_RIGHT_{i}_GC_PERCENT'], 1),
            product_size=prod,
            amplicon_long=amp_long,
            amplicon_short=amp_short))
    return {'ok': True, 'reason': '', 'pairs': pairs}
