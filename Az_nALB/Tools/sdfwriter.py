from rdkit import Chem
from rdkit.Chem import AllChem

# Your SMILES string
smiles = "OC(=O)C(N)Cc1ccc(OC(=O)COCCOCCCCCCCl)cc1"  # Vanillyl alcohol

# Convert SMILES to RDKit molecule
mol = Chem.MolFromSmiles(smiles)

# Add explicit hydrogens
mol = Chem.AddHs(mol)

# Generate 3D coordinates
AllChem.EmbedMolecule(mol)

# Write to SDF
writer = Chem.SDWriter("/Users/nathanphan/Desktop/Projects/HT_Trial2/Parameters/CET.sdf")
writer.write(mol)
writer.close()
