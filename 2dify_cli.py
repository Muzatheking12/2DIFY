import os
import argparse
from twodimensional_plot import Two_Dimensional_Interactions

def main():
    parser = argparse.ArgumentParser(
        description="Generate 2D Protein-Ligand Interactions plot."
    )

    parser.add_argument("--resname", required=True, help="Ligand residue name (e.g., UNL)")
    parser.add_argument("--pdb", required=True, help="Path to protein-ligand PDB file")
    parser.add_argument("--report", required=True, help="Path to docking report file")
    parser.add_argument("--template", required=True, help="Ligand template file (.sdf, .mol, .mol2, or SMILES)")
    parser.add_argument("--out", required=True, help="Output directory for results")
    parser.add_argument("--pad", type=float, default=0.2, help="Padding for 2D drawing (default: 0.2)")

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.out, exist_ok=True)

    # Call your main function
    Two_Dimensional_Interactions(
        resname=args.resname,
        pdb_path=args.pdb,
        report_file=args.report,
        template=args.template,
        output_dir=args.out,
        pad_setter=args.pad
    )

if __name__ == "__main__":
    main()