class FcfrpSSBonds:
    def __init__(self):
        self._idcode = None
        self._seqline = None
        self._sernum = None
        self._chainid1 = None
        self._seqnum1 = None
        self._icode1 = None
        self._chainid2 = None
        self._seqnum2 = None
        self._icode2 = None
        self._sym1 = None
        self._sym2 = None
        self._lenght = None
    
    def getFullInformation(self):
        return (self._idcode, self._seqline, self._sernum, self._chainid1, self._seqnum1, self._icode1, self._chainid2, self._seqnum2, self._icode2, self._sym1, self._sym2, self._lenght)
    
    def setIdCode(self,idcode):
        self._idcode = int(idcode)
    
    def getIdCode(self):
        return self._idcode
    
    def setSeqLine(self,SeqLine):
        self._seqline = int(SeqLine)
    
    def getSeqLine(self):
        return self._seqline
    
    def setSerNum(self,SerNum):
        self._sernum = SerNum
    
    def getSerNum(self):
        return self._sernum
    
    def setChainId1(self,chainId1):
        self._chainid1 = chainId1

    def getChainId1(self):
        return self._chainid1

    def setSeqNum1(self,SeqNum1):    
        self._seqnum1 = SeqNum1
    
    def getSeqNum1(self):    
        return self._seqnum1

    def setIcode1(self,Icode1):
        self._icode1 = Icode1

    def getIcode1(self):
        return self._icode1

    def setChainId2(self,chainId2):
        self._chainid2 = chainId2
    
    def getChainId2(self):
        return self._chainid2

    def setSeqNum2(self,SeqNum2):    
        self._seqnum2 = SeqNum2
        
    def setLength(self, lenght):
        self._lenght = lenght
    
    def getLength(self):
        return self._lenght
    
    def getSeqNum2(self):    
        return self._seqnum2

    def setIcode2(self,Icode2):
        self._icode2 = Icode2

    def getIcode2(self):
        return self._icode2
    
    def setSym1(self,Sym1):
        self._sym1 = Sym1
        
    def getSym1(self):
        return self._sym1

    def setSym2(self,Sym2):
        self._sym2 = Sym2
        
    def getSym2(self):
        return self._sym2
    

