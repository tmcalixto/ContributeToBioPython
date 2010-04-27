from FcfrpFile import FcfrpFile
from FcfrpConfig import FcfrpConfig
from Bio.PDB import *
import sys
import os

class FcfrpTopolPDB:
    def __init__(self, pdbName=None, pathPdb=None, pathOut=None):
        self._config = FcfrpConfig("GLOBAL")
        self._pdbName = pdbName
        
        if pathPdb == None:
            self._pathPdb = self._config.getParameter("Local_PDB_Files")
        else:
            self._pathPdb = pathPdb 
            
        if pathOut == None:
            self.pathOut = self._config.getParameter("Local_PDB_Files")
        else:
            self._pathOut = pathOut
            
            
    def checkPDB(self):
        qt = 0
        dell = []
        parser = PDBParser()
        name = ''
        if self._pdbName.__contains__('pdb'):
            name = self._pdbName
        else:
            name = self._pdbName + '.pdb'
            
        structure = parser.get_structure(self._pdbName, os.path.join(self._pathPdb, name)) 
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    qt = 0
                    for atom in residue:
                        qt += 1
                    if qt == 1:
                        dell.append(str(residue.id[1]))
        
        if dell.__len__() > 1:
            self.deleteResidue(dell, name)
        
            
    def deleteResidue(self, dell, name):
        pdbFile = open(os.path.join(self._pathPdb, name))
        str = ''
        strDell = ''
        for line in pdbFile.readlines():
            fields = line.split()
            
            if fields[0] == 'ATOM':
                num = fields[5]
                if num not in dell:
                    str = str + line
                    
                else:
                    strDell = strDell + line
            
            else:
                str = str + line
                    
        #print str
        print strDell
        self.saveFile(name, str)
        
    def saveFile(self, fileName, conteudo):  
        f = FcfrpFile(self._pathOut, fileName)
        f.getAsFile("w")
        f.write(conteudo)
        f.close()
                        
                