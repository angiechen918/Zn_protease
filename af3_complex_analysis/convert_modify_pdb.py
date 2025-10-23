import os
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

# Function to modify Zn1 to ZN1 and LIG_C to LIG
def modify_hetatm_line(input_pdb_file, output_pdb_file):
    with open(input_pdb_file, 'r') as infile, open(output_pdb_file, 'w') as outfile:
        for line in infile:
            if line.startswith("HETATM"):
                line = line[:12] + " ZN1 LIG " + line[23:]
            outfile.write(line)

def modify_pdb_and_rename(pdb_path):
    output_file = pdb_path.removesuffix(".pdb") + "_mod.pdb"
    modify_hetatm_line(pdb_path, output_file)
    os.remove(pdb_path)
    os.rename(output_file, pdb_path)

def process_pdb_file(row):
    pdb_path = row['pdb_path']
    cif_path = pdb_path.replace('.pdb', '.cif')

    if not os.path.exists(cif_path):
        print(f"WARNING: CIF file not found: {cif_path}")
        return
    if os.path.exists(pdb_path):
        print(f"WARNING: pdb file already exists: {pdb_path}, skipping")
        return

    try:
        parser = MMCIFParser(QUIET=True)
        io = PDBIO()
        structure = parser.get_structure('structure', cif_path)
        io.set_structure(structure)
        io.save(pdb_path)
        modify_pdb_and_rename(pdb_path)
    except Exception as e:
        print(f"Error processing {cif_path}: {e}")

def main(csv_file, start_idx, end_idx):
    df = pd.read_csv(csv_file)
    df = df.iloc[start_idx:end_idx]

    for idx, row in df.iterrows():
        process_pdb_file(row)

    print(f"rows {start_idx} to {end_idx} processed")


if __name__ == "__main__":
    import sys
    csv_file = sys.argv[1]
    start_idx = int(sys.argv[2])
    end_idx = int(sys.argv[3])
    main(csv_file, start_idx, end_idx)
