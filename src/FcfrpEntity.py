from Bio.PDB.Entity import Entity

class FcfrpEntity(Entity):
    
    def __init__(self,id):
        Entity.__init__(self, id)
        self._SeqRes = None
        self._SSBonds = None
        self._structureErrors = None
    
    def add_SeqRes(self,dicSeqRes, dicResidue):
        if self._SeqRes == None:
            raise Exception("SeqRes is None in FcfrpEntity")
        self._SeqRes.setSeqResFromDic(dicSeqRes, dicResidue)
    
    def add_SSBonds(self,dicSSBonds):
        if self._SSBonds == None:
            raise Exception("SSBonds is None in FcfrpEntity")
        self._SSBonds.addSSBonds(dicSSBonds)

    def add_StructureErrors(self,dicErrors):
        if self._structureErrors == None:
            raise Exception("Structure Errors is None in FcfrpEntity")
        self._structureErrors.addStructureErrors(dicErrors)
