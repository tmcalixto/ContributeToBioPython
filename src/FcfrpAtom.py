class FcfrpAtom:
    def __init__(self):
        self._idatom = None
        self._name = None
        self._sigla = None
        
    #Public methods
    def getIdAtom(self):
        return self._idatom
    
    def getName(self):
        return self._name
    
    def getSigla(self):
        return self._sigla
    
    def setIdAtom(self,idA):
        self._idatom = idA
    
    def setName(self,n):
        self._name = n
    
    def setSigla(self,s):
        self._sigla = s