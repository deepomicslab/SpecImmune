#!/usr/bin/env python3
"""
Read HLA G-group and P-group files, optionally truncate allele names
to a given digit resolution (2,4,6,8), and build mappings:
  - allele -> list of G-groups
  - allele -> list of P-groups

Usage:
  python3 parse_gp_group.py [-d DIGIT] hla_nom_g.txt hla_nom_p.txt
"""

import sys
import argparse
import re
import json

def truncate_allele(allele, res_digit):
    """
    Truncate the allele string to `res_digit` resolution digits:
      2 -> 1 field, 4 -> 2 fields, 6 -> 3 fields, 8 -> 4 fields.
    Preserve any trailing suffix letters (N/L/Q).
    """
    if res_digit is None:
        return allele
    if res_digit not in (2, 4, 6, 8):
        return allele
    n_fields = res_digit // 2
    m = re.match(r'^([0-9:]+)([A-Za-z]*)$', allele)
    if not m:
        return allele
    fields_part, suffix = m.groups()
    parts = fields_part.split(':')
    truncated = ':'.join(parts[:n_fields])
    return truncated + suffix

def parse_group_file(path, group_suffix, res_digit):
    """
    Parse an HLA group file (G or P) and return a dict mapping each (possibly truncated)
    allele to a list of full group names (locus+group_suffix).
    """
    mapping = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(';', 2)
            if len(parts) != 3:
                continue
            locus, allele_list, group_name = parts
            group_name = group_name.strip()
            if not group_name.endswith(group_suffix):
                continue
            full_group = locus + group_name
            for allele in allele_list.split('/'):
                allele = allele.strip()
                if not allele:
                    continue
                trunc = truncate_allele(allele, res_digit)
                full_allele = locus + trunc
                mapping.setdefault(full_allele, [])
                if full_group not in mapping[full_allele]:
                    mapping[full_allele].append(full_group)
    return mapping

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--digit', type=int, choices=[2,4,6,8], default=None,
                        help="allele resolution digits: 2,4,6,8")
    parser.add_argument('g_file', help='path to hla_nom_g.txt')
    parser.add_argument('p_file', help='path to hla_nom_p.txt')
    args = parser.parse_args()

    g_map = parse_group_file(args.g_file, 'G', args.digit)
    p_map = parse_group_file(args.p_file, 'P', args.digit)

    print(f"Parsed {len(g_map)} alleles in G-groups (resolution={args.digit})")
    print(f"Parsed {len(p_map)} alleles in P-groups (resolution={args.digit})")

    with open('g_group_map.json', 'w', encoding='utf-8') as fg:
        json.dump(g_map, fg, indent=2, ensure_ascii=False)
    with open('p_group_map.json', 'w', encoding='utf-8') as fp:
        json.dump(p_map, fp, indent=2, ensure_ascii=False)
    print("Wrote g_group_map.json and p_group_map.json")

if __name__ == '__main__':
    main()