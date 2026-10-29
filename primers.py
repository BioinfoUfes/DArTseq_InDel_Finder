"""Desenho de primers para HRM com primer3-py.

Força o par a flanquear o indel (SEQUENCE_TARGET) e restringe o amplicon à faixa
curta típica de HRM. Gera os dois amplicons (com e sem o indel) para a predição de melting.
"""
import primer3


def _design(seq_args, global_args):
    # compatível com primer3-py novo e antigo
    if hasattr(primer3, "design_primers"):
        return primer3.design_primers(seq_args, global_args)
    return primer3.bindings.designPrimers(seq_args, global_args)


def design_for_region(region, cfg):
    template = region['template']
    istart, ilen = region['indel_start'], region['indel_len']
    # alvo = o indel (se ilen=0, alveja o ponto central para ainda flanquear)
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
        return {'ok': False, 'reason': f'primer3 erro: {e}', 'pairs': []}

    n = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
    if not n:
        # retry relaxado: Tm e tamanho de produto mais amplos
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
        return {'ok': False, 'reason': f'sem par válido ({expl})', 'pairs': []}

    pairs = []
    for i in range(n):
        lstart, llen = res[f'PRIMER_LEFT_{i}']
        rstart, rlen = res[f'PRIMER_RIGHT_{i}']   # rstart = base 3' na fita top
        prod = res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
        amp_long = template[lstart:rstart + 1]
        # amplicon do alelo sem o indel (remove as bases do indel, se dentro do amplicon)
        if lstart <= istart and (istart + ilen) <= rstart + 1 and ilen > 0:
            amp_short = template[lstart:istart] + template[istart + ilen:rstart + 1]
        else:
            amp_short = amp_long  # indel fora do amplicon (não deveria ocorrer com TARGET)
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
