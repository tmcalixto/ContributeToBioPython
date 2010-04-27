from Bio.PDB.StructureBuilder import StructureBuilder
from FcfrpStructure import FcfrpStructure
from FcfrpPDBSeqRes import FcfrpPDBSeqRes
from FcfrpPDBSSBonds import FcfrpPDBSSBonds
from FcfrpPDBErrors import FcfrpPDBErrors
from FcfrpPDBAtom import FcfrpPDBAtom

class FcfrpStructureBuilder(StructureBuilder):
    def __init__(self):
        StructureBuilder.__init__(self)
    
    def init_structure(self, structure_id):
        self.structure=FcfrpStructure(structure_id)
    
    def get_structure(self):
        return self.structure
    
    def init_SeqRes(self,id):
        self.structure._SeqRes = FcfrpPDBSeqRes(id)
    
    def init_SSBonds(self,id):
        self.structure._SSBonds = FcfrpPDBSSBonds(id)
        
    def init_ErrorStructure(self,id):
        self.structure._structureErrors = FcfrpPDBErrors(id)
    
    #def init_atom(self,idcode,model,name,x,y,z,bfactor,occupancy,altloc,fullname,serial_number,idamino,resseq,icode,chainid,segid,resName=None):
     #   self.structure._Atoms = FcfrpPDBAtom(idcode,model,name,x,y,z,bfactor,occupancy,altloc,fullname,serial_number,idamino,resseq,icode,chainid,segid,resName)
        