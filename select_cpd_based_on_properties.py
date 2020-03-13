import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


def selectCompounds(NumRotatableBondsMin, NumRotatableBondsMax, MolWtMin, MolWtMax, MolLogPMin, MolLogPMax, sdfFile):
    counterSelected = 0
    suppl = Chem.SDMolSupplier(sdfFile)
    w = Chem.SDWriter(sdfFile.replace('.sdf','_selected.sdf'))
    for n,crude in enumerate(suppl):
        Name = crude.GetProp('_Name')
        all_smiles = Chem.MolToSmiles(crude)
        # need to get rid of salts etc.
        # assumption is that compound consists of longest smiles string
        # although there might be a neater version in rdkit
        smiles = max(all_smiles.split('.'), key=len)
        mol = Chem.MolFromSmiles(smiles)
        MolLogP = Descriptors.MolLogP(mol)
        MolWt = Descriptors.MolWt(mol)
        NumHAcceptors = Descriptors.NumHAcceptors(mol)
        NumHDonors = Descriptors.NumHDonors(mol)
        NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        if NumRotatableBonds >= NumRotatableBondsMin and NumRotatableBonds <= NumRotatableBondsMax:
            if MolWt > MolWtMin and MolWt < MolWtMax:
                if MolLogP > MolLogPMin and MolLogP < MolLogPMax:
                    try:
                        mol.SetProp("_Name", str(Name))
                        mol.SetProp("_MolLogP", str(MolLogP))
                        mol.SetProp("_MolWt", str(MolWt))
                        mol.SetProp("_NumRotatableBonds", str(NumRotatableBonds))
                        mol.SetProp("_NumHAcceptors", str(NumHAcceptors))
                        mol.SetProp("_NumHDonors", str(NumHDonors))
                        molH = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(molH)              # make 3D coordinates
                        AllChem.MMFFOptimizeMolecule(molH)        # optimize geometry
#                        AllChem.UFFOptimizeMolecule(molH)        # optimize geometry                        
                        w.write(molH)
                        counterSelected += 1
                    except ValueError:
                        pass
        sys.stdout.write("\r-> found %d ligands in SDF file; %d selected" %(n,counterSelected))
        sys.stdout.flush()

if __name__ == '__main__':
    NumRotatableBondsMax = 5
    NumRotatableBondsMin = 1
    MolWtMax = 500
    MolWtMin = 350
    MolLogPMax = 3
    MolLogPMin = -2
    sdfFile = 'chembl_26_0.sdf'
    selectCompounds(NumRotatableBondsMin,
                    NumRotatableBondsMax,
                    MolWtMin,
                    MolWtMax,
                    MolLogPMin,
                    MolLogPMax,
                    sdfFile)
