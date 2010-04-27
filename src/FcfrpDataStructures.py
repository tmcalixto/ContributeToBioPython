#This File contains all Data Structures sharing the project 

class TorsionAngles ():
    def __init__(self, phi, psi):
        self._phi = phi
        self._psi = psi
    
    def getPhiPsi(self):
        return (self._phi, self._psi)
    
    def getPhi(self):
        return self._phi

    def getPsi(self):
        return self._psi

class Table_t_pdb_seqres_resname:
    def __init__(self):
        self.idcode = None
        self.seqline = None
        self.seq = None
        self.idamino = None
    
    def setValues(self, idcode, seqline, seq, idamino):
        self.setIdCode(idcode)
        self.setIdAmino(idamino)
        self.setSeq(seq)
        self.setSeqLine(seqline)
    
    def setIdCode(self,c):
        self.idcode = c
        
    def setSeqLine(self,l):
        self.seqline = l
    
    def setSeq(self,s):
        self.seq = s
    
    def setIdAmino(self, a):
        self.idamino = a
    
    def getIdCode(self):
        return self.idcode
    
    def getSeqLine(self):
        return self.seqline
    
    def getSeq(self):
        return self.seq
    
    def getIdAmino(self):
        return self.idamino

class Tablet_pdb_seqres_angles:
    def __init__(self):
        self.idcode = None
        self.seqline = None
        self.seq = None
        #Object Torsional Angles
        self.torsionalAngles = None
    
    def setValues(self, idcode, seqline, seq, phi, psi):
        self.setIdCode(idcode)
        self.setSeqLine(seqline)
        self.setSeq(seq)
        self.setTorsional(TorsionAngles(phi, psi))
        
    def setIdCode(self,c):
        self.idcode = c
        
    def setSeqLine(self,l):
        self.seqline = l
    
    def setSeq(self,s):
        self.seq = s
    
    def setTorsional(self, d):
        self.torsionalAngles = d

    def setPsi(self, p):
        self.psi = p
    
    def getIdCode(self):
        return self.idcode
    
    def getSeqLine(self):
        return self.seqline
    
    def getSeq(self):
        return self.seq

    def getTorsional(self):
        return self.torsionalAngles
        
    def getPhi(self):
        return self.torsionalAngles.getPhi()

    def getPsi(self):
        return self.torsionalAngles.getPsi()

class CHIAngles():
    def __init__(self):
        self.ch1 = None
        self.ch2 = None
        self.ch3 = None
        self.ch4 = None
        self.listAngles = []
    
    def setCHIAngles(self,ch1,ch2=None,ch3=None,ch4=None ):
        self.ch1 = ch1
        if ch2 != None:
            self.ch2 = ch2
        if ch3 != None:
            self.ch3 = ch3
        if ch4 != None:
            self.ch4 = ch4
    
    def setCHIAnglesFromList(self,listAngles):
        self.listAngles = listAngles            
    
    def getListAngles(self):
        return self.listAngles
    
    def getFullAngles(self):
        return (self.ch1,self.ch2,self.ch3,self.ch4)
    

        