import json
import argparse

parser = argparse.ArgumentParser(description='Removes all the unwanted residues from a protein and keeps only the ones'
                                             'in tje provided json')
parser.add_argument(dest="pdb_input_file")
parser.add_argument(dest="pdb_output_file")
parser.add_argument(dest="keptResiduals")

def remove_residuals(pdb_input_file: str, pdb_output_file: str, keptResiduals: str):
    CHAIN_ID_SLICE = slice(21, 22)
    RESIDUE_ID_SLICE = slice(23, 26)
    with open(pdb_input_file, "r") as input_pdb:
        with open(pdb_output_file, "w") as output_pdb:
            with open(keptResiduals, "r") as residual_file:
                residuals = json.load(residual_file)
                for line in input_pdb.readlines():
                    if "ATOM" in line:
                        chain_id = line[CHAIN_ID_SLICE]
                        residue_id = line[RESIDUE_ID_SLICE]

                        if residue_id not in residuals[chain_id]:
                            continue

                    output_pdb.write(line)



if __name__ == '__main__':
    args = parser.parse_args()
    remove_residuals(args.pdb_input_file, args.pdb_output_file, args.KeptResiduals)
