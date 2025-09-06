import numpy as np
from rdkit import Chem
from io import StringIO
import matplotlib.pyplot as plt
from collections import defaultdict
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem

from PIL import Image
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor, GetSymmSSSR
import matplotlib.patches as mpatches
import pandas as pd

# Please Cite Us. Developed by Muzammil Kabier


def compute_ring_centroids(mol, coord_map):
    """
    Compute centroid of each aromatic ring using 2D coordinates.
    Returns dict: {frozenset(atom_indices): centroid_xy}
    """
    ring_centroids = {}
    ring_info = GetSymmSSSR(mol)
    for ring in ring_info:
        atom_indices = list(ring)
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_indices):
            coords = np.array([coord_map[i] for i in atom_indices])
            centroid = coords.mean(axis=0)
            ring_centroids[frozenset(atom_indices)] = centroid
    return ring_centroids

def build_residue_positions(mol, atom_dict, coord_map, interaction_dict):
    """
    Build residue positions for interactions.
    Handles normal interactions and ring centroid (pi-pi, pi-cation).
    """

    residue_positions = {}

    for (atom_symbol, transformed_idx), og_idx in atom_dict.items():
        if og_idx not in coord_map:
            continue

      
        interaction_info = interaction_dict.get(og_idx, {})
        ring_idxs = interaction_info.get("ring_idxs", None)
        interaction_type = interaction_info.get("type", None)

       
        atom_xy = coord_map[og_idx]

       
        if ring_idxs and interaction_type in ("pi-pi", "pi-cation"):
            coords = [coord_map[i] for i in ring_idxs if i in coord_map]
            if coords:  
                atom_xy = np.mean(coords, axis=0)

        
        residue_name = interaction_info.get("residue", f"res_{og_idx}")
        if residue_name not in residue_positions:
            residue_positions[residue_name] = []
        residue_positions[residue_name].append(atom_xy)

    return residue_positions

def Two_Dimensional_Interactions(resname, pdb_path, report_file, template, output_dir, pad_setter=0.2, show_plot=True):
  """
  Generate 2D Interactions:\n
  1st arg: residue name: eg UNL\n
  2nd arg: pdb-complex path\n
  3rd arg: report.txt generated from PLIP\n
  4th arg: template file or SMILES for extracting bond orders\n
  5th arg: output directory\n
  """
  try:
    print(f"PATH CHECK: {pdb_path}, {report_file}, {template}")
    print("IGNORE ERRORS, CONDITIONAL STATEMENT")
    try:
     molx = Chem.MolFromSmiles(template)
    except Exception as e:
        print("TEMPLATE IS NOT SMILES")
    if molx:
        print("SMILES TEMPLATE SET")
    else:
        print("TEMPLATE 3D FILE")

    serial_map = {}
    new_serial = 1

    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith(("HETATM", "ATOM  ")) and resname == line[17:20].strip():
                old_serial = int(line[6:11])
                serial_map[old_serial] = new_serial
                new_serial += 1

    raw_coords = []
    raw_labels = []
    raw_types = []
    raw_lig_idx = []
    raw_ring_idx = []
    raw_dist = []
    if os.path.exists(report_file):
        print("Report File Generated")
    else:
        return print("Report File Generation Failed. Cant Continue")
    
    with open(report_file, 'r') as file:
        lines = file.readlines()

    section = None
    table_started = False

    for line in lines:
        line = line.strip()

        if "**Hydrophobic Interactions**" in line:
            section = "hydrophobic"
            table_started = False
            continue
        elif "**Hydrogen Bonds**" in line:
            section = "hbond"
            table_started = False
            continue
        elif "**pi-Stacking**" in line:
            section = "pi-stacking"
            table_started = False
            continue
        elif "**Halogen Bonds**" in line:
            section = "halogens"
            table_started = False
            continue
        elif "**Salt Bridges**" in line:
            section = "salt-bridge"
            table_started = False
            continue
        elif "**pi-Cation Interactions**" in line:
            section = "pi-Cation"
            table_started = False
            continue


        if line.startswith("+") or line.startswith("| RESNR"):
            table_started = True
            continue

        if section and table_started and line.startswith("|"):
            fields = [x.strip() for x in line.strip('|').split('|')]
            try:
                if section == "hydrophobic":
                    protcoord = tuple(map(float, fields[10].split(',')))
                    lig_idxx = int(fields[7])
                    lig_idxx = serial_map.get(lig_idxx, -1)
                    distx = float(fields[6])
                    ring_indices = None

                elif section == "hbond":
                    prot_is_donor = fields[10].lower() == "true"
                    protcoord = tuple(map(float, fields[16].split(',')))
                    lig_idxx_raw = int(fields[13]) if prot_is_donor else int(fields[11])
                    lig_idxx = serial_map.get(lig_idxx_raw, -1)
                    distx = float(fields[7])
                    ring_indices = None

                elif section == "halogens":
                    protcoord = tuple(map(float, fields[15].split(',')))
                    lig_idxx_raw = int(fields[10])
                    lig_idxx = serial_map.get(lig_idxx_raw, -1)
                    distx = float(fields[7])
                    ring_indices = None

                elif section == "pi-stacking":
                    protcoord = tuple(map(float, fields[13].split(',')))
                    lig_ring_indices = list(map(int, fields[11].split(',')))
                    ring_indices = [serial_map.get(x, -1) for x in lig_ring_indices]
                    distx = float(fields[7])
                    lig_idxx = ring_indices[0]

                elif section == "pi-Cation":
                    protcoord = tuple(map(float, fields[13].split(',')))
                    lig_ring_indices = list(map(int, fields[11].split(',')))
                    ring_indices = [serial_map.get(x, -1) for x in lig_ring_indices]
                    distx = float(fields[7])
                    lig_idxx = ring_indices[0]

                elif section == "salt-bridge":
                    protcoord = tuple(map(float, fields[12].split(',')))  
                    raw_indices = list(map(int, fields[10].split(',')))   
                    lig_idxx_list = [serial_map.get(idx, -1) for idx in raw_indices]
                    distx = float(fields[7])
                    print(lig_idxx_list)
                    ring_indices = None

                    for lig_idxx in lig_idxx_list:
                        if lig_idxx == -1:
                            print(f"[WARNING] Ligand atom index not in mapping")
                            continue


                else:
                    continue

                if lig_idxx == -1:
                    print(f"[WARNING] Ligand atom index {lig_idxx_raw if 'lig_idxx_raw' in locals() else lig_idxx} not in mapping")
                    continue

                raw_coords.append(protcoord)
                raw_labels.append(f"{fields[1]}{fields[0]}")
                raw_types.append(section)
                print(raw_types)
                raw_lig_idx.append(lig_idxx)
                raw_ring_idx.append(ring_indices)
                raw_dist.append(distx)

            except (ValueError, IndexError) as e:
                print(f"[WARNING] Skipping row in section '{section}': {fields}")
                print("  → Error:", e)
                continue

    coord_dict = defaultdict(list)
    type_dict = {}
    ring_dict = {}
    dist_dict = {}

    for label, coord, lig_idx, typ, ring, dist in zip(
        raw_labels, raw_coords, raw_lig_idx, raw_types, raw_ring_idx, raw_dist):
        
        key = (label, lig_idx)
        coord_dict[key].append(coord)

        if key not in type_dict:
            type_dict[key] = typ
            ring_dict[key] = ring
            dist_dict[key] = dist

        coords, labels, lig_idxs, ring_idxs, types, dists = [], [], [], [], [], []

    for (label, lig_idx), coord_list in coord_dict.items():
        mean_coord = np.mean(coord_list, axis=0)
        coords.append(mean_coord)
        #NEW
        ring_atoms = ring_dict.get((label, lig_idx))
        if ring_atoms:
                valid_atoms = [np.array(coord_dict.get((label, idx), [])[0]) 
                               for idx in ring_atoms if (label, idx) in coord_dict]
                if valid_atoms:
                    mean_coord = np.mean(valid_atoms, axis=0)
        #END
        labels.append(label)
        lig_idxs.append(lig_idx)
        ring_idxs.append(ring_dict[(label, lig_idx)])
        types.append(type_dict[(label, lig_idx)])
        dists.append(dist_dict.get((label, lig_idx), -1))  # fallback in case of missing

    coords, labels, lig_idxs, ring_idxs, types, dists = (
        np.array(coords), labels, np.array(lig_idxs),  ring_idxs  , types, np.array(dists)
    )


    for i in range(len(coords)):
        print(f"[{types[i]}] {labels[i]}: lig_idx={lig_idxs[i]}, ring={ring_idxs[i]}, dist={dists[i]}")
    df = pd.DataFrame({
    'Interaction': types,
    'Residue': labels,
    'Distance': dists
})

    
    pathexcel = os.path.join(output_dir, r'interaction.xlsx')
    df.to_excel(pathexcel, index=False)
    filex = os.path.dirname(__file__)
    ligand_filename = f"{resname}.pdb"
    ligand_path = os.path.join(filex, ligand_filename)
    with open(pdb_path, 'r') as pdb_file:
        ligand_lines = [line for line in pdb_file if resname in line]
    with open(ligand_path, 'w') as out_file:
        out_file.writelines(ligand_lines)

    ligand_pdb_path = os.path.join(filex, f"{resname}.pdb")

    if os.path.exists(template):
        filenamex = os.path.basename(template)
        name, format = os.path.splitext(filenamex)
        print(f'FORMAT = {format}')
        
        if format.lower() == '.sdf':
            suppl = Chem.SDMolSupplier(template, sanitize=True)
            mol = suppl[0]
        elif format.lower() == '.mol2':
            mol = Chem.MolFromMol2File(template, sanitize=True)
        elif format.lower() == '.mol':
            mol = Chem.MolFromMolFile(template, sanitize=True)
    else:
        try:
            mol = Chem.MolFromSmiles(template)
            Chem.AddHs(mol)
            if mol:
                print("TEMPLATE = SMILES")
            else:
                print("INVALID TEMPLATE")
        except Exception as e:
            print(str(e))
    


    templatemol = Chem.MolFromPDBFile(ligand_pdb_path, removeHs=False, sanitize=False)
    rw_mol = Chem.RWMol(templatemol)

  
    h_indices = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetSymbol() == "H"]

    
    for idx in sorted(h_indices, reverse=True):
        rw_mol.RemoveAtom(idx)

   
    templatemol = rw_mol.GetMol()
    Chem.SanitizeMol(templatemol)



    
    if mol is None or templatemol is None:
        print("[ERROR] Could not load ligand or template molecule")
    mol = AllChem.AssignBondOrdersFromTemplate(mol, templatemol)

   
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)

    
    
    
    if mol is not None:
        
        atom_dict = {}
        for atom in mol.GetAtoms():
            og_idx = atom.GetIdx()
            transformed_idx = og_idx + 1
            atom_dict[(atom.GetSymbol(), transformed_idx)] = og_idx
            print(f"Atom symbol: {atom.GetSymbol()}, Index: {og_idx}, Transformed idx: {transformed_idx}")

        idx_to_labels = {}
        for label, typ, lig_idx in zip(labels, types, lig_idxs):
            idx_to_labels.setdefault(lig_idx, []).append((label, typ))

        AllChem.Compute2DCoords(mol)
        from rdkit.Chem.Draw import rdMolDraw2D
      
       
        rdDepictor.Compute2DCoords(mol)
        import io
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        opts = drawer.drawOptions()
        opts.padding = pad_setter
       
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png_data = drawer.GetDrawingText()

        ligand_img = Image.open(io.BytesIO(png_data))

    
        coord_map = {}
        for i in range(mol.GetNumAtoms()):
            pt = drawer.GetDrawCoords(i)
            coord_map[i] = np.array([pt.x, 500 - pt.y])  
        #NEW
        ring_centroids = compute_ring_centroids(mol, coord_map)
        #NEW
       
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(ligand_img, extent=[0, 500, 0, 500], zorder=0)

        color_map = {'hydrophobic': 'grey', 'hbond': 'blue', 'pi-stacking': 'green', 'halogens': 'cyan', 'salt-bridge': 'yellow', 'pi-Cation': 'orange'}
        style_map = {'hydrophobic': 'dashed', 'hbond': 'solid', 'pi-stacking': 'dashed', 'halogens': 'solid', 'salt-bridge': 'dotted', 'pi-cation': 'dashed' }
        legend_handles = {}

        # Optionally, plot atom positions as dots for debugging
        #for idx, (x, y) in coord_map.items():
           # ax.scatter(x, y, color='red', s=20, zorder=3, alpha=0.5)
        #NEW
        #for centroid in ring_centroids.values():
           # ax.scatter(centroid[0], centroid[1], color="purple", s=60, marker="x", zorder=4)
        #NEW
        
        residue_positions = []  
        center = np.array([250, 250])

        for (atom_symbol, transformed_idx), og_idx in atom_dict.items():
            if transformed_idx in idx_to_labels:
                aa_list = idx_to_labels[transformed_idx]

                
                anchor_point = coord_map[og_idx]

                
                for j, (label, typ) in enumerate(aa_list):
                    if typ.lower() in ["pi-stacking", "pi-cation"]:
                      
                        for ring_atoms, centroid in ring_centroids.items():
                            if og_idx in ring_atoms:
                                anchor_point = centroid
                                break

                    n = len(aa_list)
                    direction = anchor_point - center
                    direction = direction / (np.linalg.norm(direction) + 1e-8)

                    spread = np.pi if n > 1 else 0
                    angle_offset = (-spread / 2) + (spread * j / max(n - 1, 1)) if n > 1 else 0
                    rot_matrix = np.array([
                        [np.cos(angle_offset), -np.sin(angle_offset)],
                        [np.sin(angle_offset),  np.cos(angle_offset)]
                    ])
                    direction_rot = rot_matrix @ direction

                    dist_from_center = np.linalg.norm(anchor_point - center)
                    repulsion_strength = 100 + (250 - dist_from_center) * 1.2
                    repulsion_strength = np.clip(repulsion_strength, 120, 260)

                    aa_xy = anchor_point + direction_rot * repulsion_strength
                    aa_xy = np.clip(aa_xy, 0, 500)

                    residue_positions.append([aa_xy, label, typ, anchor_point])
                  
        min_dist = 40
        changed = True
        max_iter = 20
        iter_count = 0
       
        while changed and iter_count < max_iter:
            changed = False
            iter_count += 1
            for i in range(len(residue_positions)):
                for k in range(i+1, len(residue_positions)):
                    pos1, _, _, _ = residue_positions[i]
                    pos2, _, _, _ = residue_positions[k]
                    dist = np.linalg.norm(pos1 - pos2)
                    if dist < min_dist:
                        delta = (pos1 - pos2)
                        if np.linalg.norm(delta) < 1e-6:
                            delta = np.random.randn(2)
                        delta = delta / (np.linalg.norm(delta) + 1e-8)
                        move = (min_dist - dist) / 2
                        residue_positions[i][0] = np.clip(pos1 + delta * move, 20, 480)
                        residue_positions[k][0] = np.clip(pos2 - delta * move, 20, 480)
                        changed = True

      
        for aa_xy, label, typ, atom_xy in residue_positions:
            ax.plot([atom_xy[0], aa_xy[0]], [atom_xy[1], aa_xy[1]],
                    color=color_map.get(typ, 'black'),
                    linestyle=style_map.get(typ, 'solid'),
                    linewidth=2)

            circle = plt.Circle(aa_xy, 20, color='lightgreen', zorder=2)
            ax.add_patch(circle)
            ax.text(aa_xy[0], aa_xy[1], f"{label}", ha='center', va='center', fontsize=8, zorder=3)
         
            if typ not in legend_handles:
                legend_handles[typ] = mpatches.Patch(color=color_map.get(typ, 'black'), label=typ)

       
        ax.set_xlim(-50, 550)
        ax.set_ylim(-50, 550)
        ax.axis('off')
        plt.title('2D Protein Ligand Interactions')
        plt.legend(handles=list(legend_handles.values()), loc='best')

        
        plt.savefig(os.path.join(output_dir, f"2D-Interactions.png"), bbox_inches='tight', transparent=True, dpi=300)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        os.remove(ligand_pdb_path)
        if show_plot == True:
         plt.show()
        else:
          plt.close()

  except Exception as e:
       print(str(e))
