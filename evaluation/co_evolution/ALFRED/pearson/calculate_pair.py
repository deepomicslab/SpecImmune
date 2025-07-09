import os
import glob
import random
import re
import pandas as pd
from scipy.stats import pearsonr

def extract_chrom_from_path(filepath):
    parts = os.path.normpath(filepath).split(os.sep)
    for part in parts[::-1]:
        if part.startswith('SI') or part.endswith('.txt'):
            continue
        return part
    m = re.search(r'chr(\w+)', filepath)
    if m:
        return m.group(1)
    return "unknown"

def read_alfred_freq_file(filename):
    records = []
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('Population|'):
            header_idx = i
            break
    if header_idx is None:
        return pd.DataFrame()
    data_start = header_idx + 1
    while data_start < len(lines) and lines[data_start].strip().startswith('-'):
        data_start += 1
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith('-'):
            continue
        parts = line.split('|')
        if len(parts) >= 6:
            pop, sid, date, n2, ref_str, alt_str = [x.strip() for x in parts[:6]]
        elif len(parts) == 5:
            pop, sid, date, n2, freq_str = [x.strip() for x in parts]
            if freq_str == "1.0":
                ref_str, alt_str = "1.0", "0.0"
            elif freq_str == "0.0":
                ref_str, alt_str = "0.0", "1.0"
            else:
                continue
        else:
            continue
        try:
            ref_freq = float(ref_str)
            alt_freq = float(alt_str)
        except Exception:
            continue
        records.append({
            'Population': pop,
            'ref_freq': ref_freq,
            'alt_freq': alt_freq
        })
    return pd.DataFrame(records)

def calc_population_corr(df1, df2):
    pops1 = set(df1['Population'])
    pops2 = set(df2['Population'])
    overlap = pops1 & pops2
    if not overlap:
        return None, None, 0, []
    df1_sub = df1[df1['Population'].isin(overlap)].groupby('Population').mean().reset_index()
    df2_sub = df2[df2['Population'].isin(overlap)].groupby('Population').mean().reset_index()
    merged = pd.merge(df1_sub, df2_sub, on='Population', suffixes=('_1', '_2'))
    if merged.empty:
        return None, None, 0, []
    try:
        corr, pval = pearsonr(merged["ref_freq_1"], merged["ref_freq_2"])
    except Exception:
        corr, pval = None, None
    return corr, pval, len(merged), list(merged['Population'])

def main():
    data_root = "alfred_freqs2"
    n_random = 200

    files = glob.glob(os.path.join(data_root, "*", "*.txt"))
    if len(files) < 2:
        print("At least 2 SNP files required.")
        return

    random.seed(42)
    files_sample = random.sample(files, min(n_random, len(files)))
    file_chrom = {f: extract_chrom_from_path(f) for f in files_sample}
    dfs = {f: read_alfred_freq_file(f) for f in files_sample}

    results = []
    for i in range(len(files_sample)):
        for j in range(i+1, len(files_sample)):
            f1, f2 = files_sample[i], files_sample[j]
            chrom1, chrom2 = file_chrom[f1], file_chrom[f2]
            if chrom1 == chrom2:
                continue  # Skip same chromosome
            df1, df2 = dfs[f1], dfs[f2]
            corr, pval, n_overlap, overlap_pops = calc_population_corr(df1, df2)
            if n_overlap < 5:
                continue  # Skip if less than 5 overlapping populations
            results.append({
                "file1": os.path.basename(f1),
                "chrom1": chrom1,
                "file2": os.path.basename(f2),
                "chrom2": chrom2,
                "n_overlap_pops": n_overlap,
                "corr_ref": corr,
                "pval": pval,
                "overlap_pops": ",".join(overlap_pops)
            })
            print(f"{os.path.basename(f1)}({chrom1}) vs {os.path.basename(f2)}({chrom2}) | overlap_pops: {n_overlap} | corr_ref: {corr}")

    out_tsv = "alfred_200_random_interchrom_pairwise_population_corr.tsv"
    pd.DataFrame(results).to_csv(out_tsv, sep='\t', index=False)
    print(f"Results saved to {out_tsv}")

if __name__ == "__main__":
    main()