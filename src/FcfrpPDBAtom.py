from Bio.PDB.Atom import Atom  

import numpy

class FcfrpPDBAtom(Atom):
    
    def __init__(self,idcode,model,name,x,y,z,bfactor,occupancy,altloc,fullname,serial_number,idamino,resseq,icode,chainid,segid,resName=None):
        coord = coord=numpy.array((x, y, z), 'f')

        Atom.__init__(self, str(name), coord, bfactor, occupancy, str(altloc), str(fullname), int(serial_number))
        self._model = model
        self._idcode = idcode
        self._idamino = idamino
        self._resseq = resseq
        self._icode = str(icode)
        self._chainid = str(chainid)
        self._segid = segid
        #Used when we have the resname. 
        self._resname = resName
        
    
    def getModel(self):
        return self._model 
    
    def getIdCode(self):
        return self._idcode
    
    def getIdAmino(self):
        return self._idamino 
    
    def getResSeq(self):
        return self._resseq
    
    def getIcode(self):
        return self._icode
    
    def getChainId(self):
        return self._chainid
    
    def getSeqId(self):
        return self._segid
        
    def getResName(self): 
        return self._resname