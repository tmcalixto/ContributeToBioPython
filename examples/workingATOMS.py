import sys
from FcfrpPDB import getStructure

def main():
    modelChoose = None
    id = sys.argv[1]
    directory = sys.argv[2] 
    structure = getStructure(id, directory, modelChoose)

    #Like biopython
    for model in structure:
        for chain in model:
            for residue in chain:
                print residue.get_resname()
                for atom in residue:
                    print atom.get_name(), atom.get_coord()

main()