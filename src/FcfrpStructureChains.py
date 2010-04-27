import os
from FcfrpPDB import getStructureChains,saveStructure, getStructure
from FcfrpPDBSeqRes import FcfrpPDBSeqRes
from FcfrpFile import FcfrpFile
from FcfrpTitration import FcfrpTitration
from FcfrpCapacitance import FcfrpCapacitance

class FcfrpStructureChains():
    
    def __init__(self,path, fileNameSet):
        self._path = path
        self.fileSet = FcfrpFile(path, fileNameSet)
        
    def getPDBIdFromName(self,name):
        return str(name).split("_")[0]
        
    def splitPDBChains(self,id):
        return getStructureChains(id)
    
    def getDirectory(self,id):
        #Check if there isn't the path, it'll create and return it
        path = os.path.join(self._path, id) 
        if os.path.exists(path) == False:
            os.mkdir(path)
        return path
    
    def getPathFileName(self, name):
        id = self.getPDBIdFromName(name)
        path = self.getDirectory(id)
        FileName = name + ".pdb"
        return os.path.join(path,FileName)  
    
    def saveStructures(self,dicStructures):
        # Save each structure stored in the dicStructures dictionary.
        for k,structure in dicStructures.iteritems():
            pathFileName = self.getPathFileName(k)
            saveStructure(structure,pathFileName)
    
    def buildFileSetProteins(self,id,dicStructures):
        f = self.fileSet.getAsFile('a')
        f.write(";"+id+"\n")#The main pdb
        for k,v in dicStructures.iteritems():
             f.write(" "+k+"\n")
    
    def executeTitration(self,dicStructures, table):
        for k, structure in dicStructures.iteritems():
            id = self.getPDBIdFromName(k)
            path = self.getDirectory(id)
            FcfrpTitration(3,structure=structure,TableId=table,PDBid=k,path=path,FileName=None).execute(path)
            
    def executeCapacitance(self,dicStructures, table):
        for k, structure in dicStructures.iteritems():
            id = self.getPDBIdFromName(k)
            path = self.getDirectory(id)
            FcfrpCapacitance(3,structure=structure,TableId=table,PDBid=k,path=path,FileName=None).execute(path)
            
