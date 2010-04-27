from Bio.PDB import *
from FcfrpFile import *
from scipy import *
import os

class FcfrpUtils():  
     
     
    def __init__(self, dir=None, name="None"):
        if dir != None:
            self.pqrFile = FcfrpFile(dir,name+".pqr")
            
    def getRandomNumberFromPath(self, path):
        # Give a path, this method returns Path Random number.
        listPath = path.split("/")
        tam = listPath.__len__()-1
        i = 0
        while i < tam:
            n = listPath[tam-i]
            if n != "":
                 return n 
            i = i + 1
                    
    # Calculate the geometric center of protin
    def centroGeometricoProtein(self):
        quantidade = 0
        
        for line in self.pqrFile.getAsFile().readlines():
            x = 0.0
            y = 0.0
            z = 0.0
            
            coord = line.split()
            
            if coord[0] == "ATOM":
                quantidade = quantidade + 1
                x = x + float(coord[5])
                y = y + float(coord[6])
                z = z + float(coord[7])
        
        x = x/quantidade
        y = y/quantidade
        z = z/quantidade
        
        centroGeometrico = (x, y, z)
        self.pqrFile.close()
        
        return centroGeometrico
    
    
    # Calculate the ray of the protein
    # It is the greater distance from atom to geometric center of protein
    # cg -> centro geometrico
    def raioProtein(self, cg):
        maiorDistancia = 0.0
        distance = 0.0
        
        for line in self.pqrFile.getAsFile().readlines():
            fields = line.split()
            
            if fields[0] == "ATOM":
                coord = (fields[5], fields[6], fields[7])
                distance = self._getDistance(cg, coord)
                if distance > maiorDistancia:
                    maiorDistancia = distance
                distance = 0.0
            
        return maiorDistancia
    
        # Calculate the distance between two points in the tridimensional space
    def _getDistance(self, p1, p2):
        x1 = float(p1[0])
        y1 = float(p1[1])
        z1 = float(p1[2])
        
        x2 = float(p2[0])
        y2 = float(p2[1])
        z2 = float(p2[2])
        
        dx = float(x2 - x1)
        dy = float(y2 - y1)
        dz = float(z2 - z1)
        
        d = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2))
        return d
   
    def getVolumeRadiusProtein(self, dicProteinChainSeq, database):
        #dicProteinChainSeq is a dictionary which has its key chain and its value chain residues
        #Computes volume and radius of all chains of protein
        volRet = 0.0
        raioRet = 0.0
        for chain, sequence in dicProteinChainSeq.iteritems():
            vol, raio = self.getVolumeRadiusProteinChain(dicProteinChainSeq, database, chain)
            volRet = volRet + vol  
            raioRet = raioRet + raio
        return volRet, raioRet
    
    def getVolumeRadiusProteinChain(self, dicProteinChainSeq, database, chainId):
        #Computes volume and radius of chain of protein
        vol = 0.0
        raio = 0.0 
        sequence = dicProteinChainSeq[chainId]
        for s in sequence:
            #Computes the protein chain volume
            vol = vol + database.getFullInformationAminoacid(database.getAminoIdFromDatabase(s))[8] #Radius amino acid
        #Computes the protein chain radius
        raio = raio + math.pow(((3*vol)/12.56),1.0/3)
        return vol, raio
    
    def getMassProteinChain(self, dicProteinChainSeq, database, chainId):
        #Computes the mass of chain of protein
        mass = 0.0 
        sequence = dicProteinChainSeq[chainId]
        for s in sequence:
            #Computes the protein chain mass 
            mass = mass + database.getFullInformationAminoacid(database.getAminoIdFromDatabase(s))[6] #Mass of protein
        return mass
    
    def getChainStr(self,listchain):
        #Receives a list of chain (A;B;C) and returns a string such as 'A','B','C'
        # This result will be used in select command
        ChainsPBD = ""
        for c in listchain:
            c = "'" + c + "'"
            if ChainsPBD.__len__() == 0:
                ChainsPBD = ChainsPBD + str(c)
            else:
                ChainsPBD = ChainsPBD + "," + str(c)
        return ChainsPBD      