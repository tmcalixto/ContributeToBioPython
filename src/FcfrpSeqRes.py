#This class represents of the database Table t_pdb_seqres
class FcfrpSeqRes:
    def __init__(self):
        self._idcode = None
        self._seqline = None
        self._serNum = None
        self._chainid = None
        self._numRes = None
        
    def setIdCode(self,idCode):
        self._idcode = int(idCode)
    
    def setSeqLine(self,l):
        self._seqline = int(l)
    
    def setSerNum(self,S):
        self._serNum = int(S)
    
    def setChainId(self,C):
        self._chainid =str(C)
    
    def setNumRes(self,N):
        self._numRes = int(N)        
            
    def getIdcode(self):
        return self._idcode
    
    def getSeqLine(self):
        return self._seqline
    
    def getSerNum(self):
        return self._serNum
    
    def getChainId(self):
        return self._chainid
    
    def getNumRes(self):
        return self._numRes
    
    def copy(self,SeqRes):
        self._idcode = SeqRes.getIdcode()
        self._seqline = SeqRes.getSeqLine()
        self._serNum = SeqRes.getSerNum()
        self._chainid = SeqRes.getChainId()
        self._numRes = SeqRes.getNumRes()
         