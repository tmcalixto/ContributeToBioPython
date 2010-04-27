#This class represents of the database Table t_pdb_seqres_resname or SEQRES section from PDB file.
# The database table stores idAmino and SeqLine values. Thus, to get another information, is necessary to work with 
# FcfrpDatabase class which its method getFullInformationAminoacid can bring them.
# When work with PDB file, the only information can bring is the of residue with 3 letters.
# Therefore, when need to work with Database, is necessary indicate FcfrpDatabase object because, its constructor
# calls some operations with database.  

class FcfrpResidue:
    def __init__(self,idAmino=None,seq=0,name=None,fcfrpDatabase=None,chainId=None): #name,l1,l3,
        
        if fcfrpDatabase != None: #Working with database
            database = fcfrpDatabase       
            v = database.getFullInformationAminoacid(idAmino)
            name = str(v[1]) 
            l1 = str(v[2])
            l3 = str(v[3])
            self.setValues(seq, idAmino, name, l1, l3,chainId)
        else: #Working with PDB File. Therefore, the number of information is low than Database
            self.setValues(None, seq, "", "", str(name),chainId)

    def setValues(self,seq,idAmino,name,l1,l3,chainId):
        self._seq = seq
        self._idAmino = idAmino        
        self._ResName = name
        self._idl1 = l1
        self._idl3 = l3
        self._chainId = chainId 
    
    def getChainId(self):
        return self._chainId
    
    def getResName(self):
        return self._ResName
    
    def getIdl1(self):
        return self._idl1    

    def getIdl3(self):
        return self._idl3
    
    def getSeq(self):
        return self._seq
    
    def getIdAmino(self):
        return self._idAmino
        