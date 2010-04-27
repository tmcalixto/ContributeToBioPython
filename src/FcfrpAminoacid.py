class FcfrpAminoacid:

    def __init__(self):
        self._idAmino = None
        self._name  =  None
        self._idl1 = None
        self._idl3 = None
        self._tipo = None
        self._massMolar = None
        self._massMolarCalc = None
        self._radius = None
        self._volume = None
    
    def getFullInformation(self):
        #Return a Tuple with whole information about amino
        return (self._idAmino,self._name,self._idl1,self._idl3,self._tipo,self._massMolar,self._massMolarCalc, self._radius, self._volume)
    
    #Get methods
    def getIdAmino(self):
        return self._idAmino
    
    def getName(self):
        return self._name
    
    def getIdl1(self):
        return self._idl1    

    def getIdl3(self):
        return self._idl3
    
    def getTipo(self):
        self._tipo
    
    def getMassMolar(self):
        return self._massMolar
    
    def getMassMolarCalc(self):
        return self._massMolarCalc
    
    def getRaduis(self):
        return self._radius
    
    def getVolume(self):
        return self._volume

    # Set methods
    def setIdAmino(self,id):
        self._idAmino = int(id)
        
    def setName(self,n):
        self._name = str(n)
    
    def setIdl1(self,id1):
        self._idl1 = str(id1)
        
    def setIdl3(self,id3):
        self._idl3 = str(id3)
        
    def setTipo(self,tp):
        self._tipo = str(tp)
        
    def setMassMolar(self,mm):
        self._massMolar = mm
        
    def setMassMolarCalc(self,mmc):
        self._massMolarCalc = mmc
    
    def setRadius(self,r):
        self._radius = r
    
    def setVolume(self,v):
        self._volume = v
