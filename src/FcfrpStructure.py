#Based on Biopython project
from FcfrpEntity import FcfrpEntity

class FcfrpStructure(FcfrpEntity):
    
    def __init__(self, id):
        FcfrpEntity.__init__(self, id)

    def haveErrorsFound(self):
        return True #(self.get_StructureErrors().__len__() > 0)
     
    def get_SeqRes(self):
        return self._SeqRes
    
    def get_SSBonds(self):
        return self._SSBonds
    
    def get_StructureErrors(self):
        return self._structureErrors

    # The methods get_chains, get_residues and get_atoms pertain to biopython project
    # These methods were copied by us because this class hasn't inherited from biopython Structure class        
    def get_chains(self):
        for m in self:
            for c in m:
                yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c :
                yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a

        