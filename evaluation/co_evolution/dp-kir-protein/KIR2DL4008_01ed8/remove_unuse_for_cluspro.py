#!/usr/bin/env python3
import argparse

def fix_mse_to_met(line: str) -> str:
    """
    Convert MSE residue to MET in-place for a PDB ATOM/HETATM line.
    Also converts atom name 'SE' to 'SD' and element 'SE' to 'S' when present.
    """
    # PDB columns (1-based):
    # 18-20 resName, 13-16 atomName, 77-78 element
    resn = line[17:20]
    if resn != "MSE":
        return line

    # Convert record to ATOM (ClusPro tends to like standard protein ATOM records)
    if line.startswith("HETATM"):
        line = "ATOM  " + line[6:]

    # resName: MSE -> MET
    line = line[:17] + "MET" + line[20:]

    # atom name: "SE" -> "SD" (MET has SD)
    atom_name = line[12:16]
    if atom_name.strip() == "SE":
        # keep alignment width=4
        new_atom = atom_name.replace("SE", "SD")
        line = line[:12] + new_atom + line[16:]

    # element: "SE" -> "S"
    if len(line) >= 78:
        elem = line[76:78]
        if elem.strip().upper() == "SE":
            # element field is 2 chars, right-justified typically
            line = line[:76] + " S" + line[78:]

    return line

def main():
    ap = argparse.ArgumentParser(
        description="Clean PDB for ClusPro: remove unknown residues (e.g., ESP) and optionally drop HETATM."
    )
    ap.add_argument("input_pdb", help="Input PDB file")
    ap.add_argument("output_pdb", help="Output cleaned PDB file")
    ap.add_argument("--remove", default="ESP",
                    help="Comma-separated residue names to remove anywhere (default: ESP). Example: ESP,SEP,TPO")
    ap.add_argument("--drop-hetatm", action="store_true", default=True,
                    help="Drop all HETATM records (default: True)")
    ap.add_argument("--keep-hetatm", action="store_false", dest="drop_hetatm",
                    help="Keep HETATM records (disables --drop-hetatm)")
    ap.add_argument("--mse-to-met", action="store_true", default=True,
                    help="Convert MSE to MET (default: True)")
    ap.add_argument("--no-mse-to-met", action="store_false", dest="mse_to_met",
                    help="Do not convert MSE to MET")
    ap.add_argument("--altloc", default="A",
                    help="Keep only this altLoc or blank. Default keeps ' ' and 'A'. Set to e.g. B to keep ' ' and 'B'.")
    args = ap.parse_args()

    remove_set = {r.strip().upper() for r in args.remove.split(",") if r.strip()}

    kept_any = 0
    with open(args.input_pdb, "r", encoding="utf-8", errors="ignore") as fin, \
         open(args.output_pdb, "w", encoding="utf-8") as fout:

        for line in fin:
            rec = line[:6]

            if rec in ("ATOM  ", "HETATM"):
                resn = line[17:20].strip().upper()

                # altLoc filter (column 17)
                alt = line[16:17]
                if alt not in (" ", args.altloc):
                    continue

                # remove specified residue names
                if resn in remove_set:
                    continue

                # drop all HETATM (unless keeping)
                if rec == "HETATM" and args.drop_hetatm:
                    # optionally keep/convert MSE
                    if args.mse_to_met and resn == "MSE":
                        line = fix_mse_to_met(line)
                    else:
                        continue
                else:
                    if args.mse_to_met and resn == "MSE":
                        line = fix_mse_to_met(line)

                fout.write(line)
                kept_any += 1
                continue

            # drop ANISOU (often tied to removed atoms; simplest is to discard)
            if rec == "ANISOU":
                continue

            # drop TER; not necessary for ClusPro, and can become inconsistent after filtering
            if rec == "TER   ":
                continue

            # keep other records (HEADER, TITLE, REMARK, END, etc.)
            fout.write(line)

        if kept_any == 0:
            raise SystemExit("No atoms kept! Check your input and --remove/--drop-hetatm options.")

if __name__ == "__main__":
    main()