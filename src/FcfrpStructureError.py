class FcfrpStructureError:
    
    def __init__(self):
          self._idcode  = None
          self._iderror = None
          self._seq     = None
          self._message = None

    def getFullInformation(self):
        return (self._idcode, self._iderror, self._seq, self._message)

    def setIdCode(self,idcode):
        self._idcode  = idcode
    
    def getIdCode(self):
        return self._idcode
    
    def setIdError(self,e):
        self._iderror = e
    
    def getIdError(self):
        return self._iderror
    
    def setSeq(self,s):
        self._seq = s
    
    def getSeq(self):
        return self._seq
    
    def setMessage(self, m ):
        self._message = m
    
    def getMessage(self):
        return self._message
