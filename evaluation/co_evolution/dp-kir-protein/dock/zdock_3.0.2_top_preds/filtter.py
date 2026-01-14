#!/usr/bin/env python3
import sys

STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}

def field(line, a, b):
    return (line[a:b] + " "*(b-a))[:(b-a)]

def main(inp, outp):
    kept = 0
    with open(inp, "r", errors="ignore") as f, open(outp, "w") as g:
        for line in f:
            rec = field(line, 0, 6)
            if rec in ("ATOM  ", "HETATM"):
                resn = field(line, 17, 20).strip().upper()
                if resn not in STANDARD_AA:
                    continue
                # 统一成 ATOM
                if rec != "ATOM  ":
                    line = "ATOM  " + line[6:]
                g.write(line.rstrip("\n") + "\n")
                kept += 1
            elif rec == "TER   ":
                continue
            elif rec.strip() == "END":
                continue
        g.write("END\n")
    if kept == 0:
        raise SystemExit("No protein atoms kept; check input PDB.")
    print(f"[OK] kept {kept} atom lines -> {outp}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: pdb_keep_standard_protein.py in.pdb out.pdb")
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])