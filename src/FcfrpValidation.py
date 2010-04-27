from scipy import *
from FcfrpAminoTopol import *
from FcfrpFile import FcfrpFile
import os
 
class FcfrpValidation():   
    
    def __init__(self, structure,pathFileName, chainId = None):
        self._structure = structure
        self._descricaoAmino = FcfrpAminoTopol()
        dir,name = os.path.split(pathFileName)
        self._file = FcfrpFile(dir,name)
        # Used to remove heteroatoms
        self._valid_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE" , "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        self._dicResidueFromSeqres = {}
        self._dicResidueFromAtom = {}
        self._dicMissingResidues = {}
        self._dicMissingAtoms = {}
        self._dicDuplicatedAtoms = {}
        self._listUnknowResidues = []
        self._erros = []
        self._chainId = chainId
    
    # Check if have missing residues and/or atoms on structure
    def checkStructure(self):
        self.setResiduesFromSeqres()
        self.setResiduesFromAtom()
        i = self.getDifferenceBetweenDictionaries(self.getResiduesFromSeqres(), self.getResiduesFromATOM())
        # HAVE DIFFERENCE
        # Find missing residues
        if i == 1:
            self.setPositionMissingResidues(self.getResiduesFromSeqres(), self.getResiduesFromATOM())        
        
        self.findMissingAtoms()
        self.findDuplicatedAtoms()

        
    # Build a dictionary with all residues per chain
    # The residues are come from SEQRES
    def setResiduesFromSeqres(self):
        residuos = []        
        for line in self._file.readLines():
            i = 0
            field = line.split()
            if field[0] == "SEQRES":
                if self.checkChain(field[2], self._chainId) == True:
                    for res in field:
                        if res in self._valid_residues:
                            cadeia = field[2]
                            if i > 3:
                                residuos.append(str(res))
                        else:
                            if i > 3:
                                self._listUnknowResidues.append(res) # Residue 
                                self._listUnknowResidues.append(field[2]) # chain
                        i = i + 1
                    
                    if self._dicResidueFromSeqres.has_key(cadeia) == 1:
                        tmp = []
                        tmp = self._dicResidueFromSeqres[cadeia]
                        for x in residuos:
                            tmp.append(x)
                        self._dicResidueFromSeqres[cadeia] = tmp
                        residuos = []
                    else:
                        self._dicResidueFromSeqres[cadeia] = residuos
                        residuos = []
        
        
    # Build a dictionary with all residues per chain
    # The residues are come from ATOM
    def setResiduesFromAtom(self):
        listResidues = []        
        for model in self._structure:
            for chain in model:
                if self.checkChain(chain.id, self._chainId) == True:
                    for residue in chain:
                        if residue.get_resname() in self._valid_residues:
                            listResidues.append(str(residue.get_resname()))                        
                    self._dicResidueFromAtom[chain.id] = listResidues
                    listResidues = []

    
    # Return residues dictionary from SEQRES                         
    def getResiduesFromSeqres(self):
        return self._dicResidueFromSeqres
    
    
    # Return residues dictionary from ATOM
    def getResiduesFromATOM(self):
        return self._dicResidueFromAtom
    
    
    # Return a dictionary with missing residues
    def getMissingResidues(self):
        return self._dicMissingResidues
    
    
    # Return a dictionary with missing atoms
    def getMissingAtoms(self):
        return self._dicMissingAtoms
    
    # Return a list with duplicated atoms
    def getDuplicatedAtoms(self):
        return self._dicDuplicatedAtoms
    
    # Return a list with unknown residues
    def getUnknownResidues(self):
        return self._listUnknowResidues
    
    
    # Return 1 if there are missing residues, and 0 if there aren't missing residues
    def haveMissingResidues(self):
        if self._dicMissingResidues.__len__() > 0:
            return 1
        else:
            return 0
        
    # Return 1 if there are missing atoms, and 0 if there aren't missing atoms
    def haveMissingAtoms(self):
        if self._dicMissingAtoms.__len__() > 0:
            return 1
        else:
            return 0
    
    # Return 1 if there are duplicated atoms, and 0 if there aren't duplicated atoms    
    def haveDuplicatedAtoms(self):
        if self._dicDuplicatedAtoms.__len__() > 0:
            return 1
        else:
            return 0
    
    # Return 1 if there are unknown residues, and 0 if there aren't unknown residues
    def haveUnknownResidues(self):
        if self._listUnknowResidues.__len__() > 0:
            return 1
        else:
            return 0
    
    
    def hasErrors(self):
        if self.haveMissingResidues() == 1:
            self._erros.append(1)
            
        if self.haveMissingAtoms() == 1:
            self._erros.append(2)
            
        if self.haveDuplicatedAtoms() == 1:
            self._erros.append(3)
            
        if self.haveUnknownResidues() == 1:
            self._erros.append(4)
            
        return self._erros
             
            
    # Return different length between two dictionaries (SEQRES and ATOM)
    def getDifferenceBetweenDictionaries(self, seqres, atom):
        listSeqres = []
        listAtom = []
        diff = 0
        for chain in seqres:
            listSeqres = seqres[chain]
            listAtom = atom[chain]
            if listSeqres.__len__() - listAtom.__len__() != 0:
                diff = 1
            listSeqres = []
            listAtom = []
        return diff 
    
    
    # Set position of missing residues
    def setPositionMissingResidues(self, dicResSeqres, dicResAtom):     
        for chain in dicResSeqres:
            listResFromSeqres = []
            listResFromAtom = []
            listMiss = []    
            
            listResFromSeqres = dicResSeqres[chain]
            listResFromAtom = dicResAtom[chain]
            
            i = 0
            j = 0   
                     
            while i < listResFromSeqres.__len__():
                if listResFromSeqres[i] == listResFromAtom[j]:
                    i = i + 1
                    if j < listResFromAtom.__len__()-1:
                        j = j + 1
                else:
                    listMiss.append(str(listResFromSeqres[i]))
                    listMiss.append(str(i+1))
                    self._dicMissingResidues[chain] = listMiss
                    i = i + 1

                            
    # Find for missing atoms in a structure
    def findMissingAtoms(self):
        result = []
        for model in self._structure:
            for chain in model:
                if self.checkChain(chain.id, self._chainId) == True:
                    for residue in chain:
                        if residue.get_resname() in self._valid_residues:                       
                            result = self.checkAtomsFromResidue(residue.get_resname(), self.getAtomsFromResidueInATOM(chain.id, residue.get_id()))
                            if result.__len__() > 0:
                                # result.append(chain.id)
                                # result.append(residue.get_id())
                                key_comp = residue.get_id()
                                self._dicMissingAtoms[chain.id, key_comp[1], residue.get_resname()] = result 
                            
            self.getAtomsFromMissingResidues()
    
    
    # Return a list of atom from residue - Search in ATOM
    def getAtomsFromResidueInATOM(self, chainRef, residueRef):
        atoms = []
        for model in self._structure:
            for chain in model:
                if chain.id == chainRef:
                    residue = chain[(residueRef[0],residueRef[1],residueRef[2])]
                    for atom in residue:
                        atoms.append(str(atom.get_name()))
        return atoms
    
    
    # check the amount atoms from residue
    def checkAtomsFromResidue(self, residue, atoms):
        dicRes = self._descricaoAmino.getDicionario(residue)
        atomRes = self.organizaAtomos(atoms)
        
        diff = {}
        missings = []
        
        for key in dicRes:
            if atomRes.has_key(key) == 1:
                if key != 'H':
                    if key == 'O':
                        diff[key] = str( int(dicRes[key])-1 - int(atomRes[key])) # Peptide Bond
                    else:
                        diff[key] = str(int(dicRes[key]) - int(atomRes[key]))    # Peptide Bond

                    if int(diff[key]) > 0: # It may to have duplicates
                        #          [0] Residue,    [1]amount,      [2]atom 
                        missings = [str(diff[key]), str(key)]            
        
        return missings
                    
        
    # Get all atoms from missing residues
    def getAtomsFromMissingResidues(self):
        missRes = []
        missResPosition = []
        atoms = []
        
        for chain in self._dicMissingResidues:
            missRes = []
            missResPosition = []

            i = 0
            res = self._dicMissingResidues[chain]
            while i < res.__len__():    
                missRes.append(res[i])
                missResPosition.append(res[i + 1])
                i = i + 2
                       
            j = 0
            for rs in missRes:
                dicRes = self._descricaoAmino.getDicionario(rs)
            
                # atoms.append(rs)
                for a in dicRes:
                    atoms.append(dicRes[a])
                    atoms.append(a)
                                    
                self._dicMissingAtoms[chain, str(missResPosition[j]), rs] = atoms
                atoms = []
                j = j + 1    
 
    # Build a dictionary with atoms and amount
    def organizaAtomos(self, atomos):
        dicAtomos = {}
        for atm in atomos:
            tmp = str(atm) # tmp and atom are use to get only the first character to each kind of atoms
            atom = tmp[0]
            if dicAtomos.has_key(str(atom)) == 1:
                valor = int(dicAtomos[str(atom)])
                valor = valor + 1
                dicAtomos[str(atom)] = str(valor)
            else:
                dicAtomos[str(atom)] = str(1)
                  
        return dicAtomos
    
    
    # Look for duplicated atoms
    def findDuplicatedAtoms(self):
        # [0]Residue [1]Chain [2]id from residue [3]Atom [4]Alternative Local
        result = [] 
        
        for model in self._structure:
            for chain in model:
                if self.checkChain(chain.id, self._chainId) == True:
                    for residue in chain:
                        for atom in residue:
                            result = self.checkDuplicatedAtomsInResidue(chain, residue, atom)
                            if result.__len__() > 0:
                                key_comp = residue.get_id()
                                key = chain.id, key_comp[1], residue.get_resname(), key_comp[2]
                            
                                if self._dicDuplicatedAtoms.has_key(key) == 1:
                                    tmp = []
                                    tmp = self._dicDuplicatedAtoms[key]
                                    for x in result:
                                        tmp.append(x)
                                        self._dicDuplicatedAtoms[key] = tmp
                                        result = []
                                else:
                                    self._dicDuplicatedAtoms[key] = result
                                    result = []
        
                            
    def checkDuplicatedAtomsInResidue(self, chain, residue, atom):
        x = ""
        duplicateds = [] 
        x = str(atom.get_altloc())
        if x != " ":
            duplicateds.append(str(atom.get_name()))
            duplicateds.append(str(atom.get_altloc()))
        
        return duplicateds

    def checkChain(self, chainIdFile, chainIdSelected):
        cont = False
        if chainIdSelected == None:
            cont = True
        elif chainIdFile == chainIdSelected:
            cont = True
        else:
            cont = False
        return cont 