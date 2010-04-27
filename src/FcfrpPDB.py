import urllib
from Bio import File
from FcfrpFile import FcfrpFile
from FcfrpDatabase import FcfrpDatabase
import os
from Bio.PDB.PDBParser import PDBParser
from FcfrpPDBParser import FcfrpPDBParser
from FcfrpConfig import FcfrpConfig
from FcfrpPDBBio import FcfrpPDBBio 
from FcfrpValidation import FcfrpValidation
from FcfrpShowErrosDetails import FcfrpShowErrosDetails
 
def getAmountModel(structure):
    #Count the model number of structure. For NRM, this number is more than 1
    i = 0
    for model in structure:
        i = i + 1
    return i
    
def saveStructure(structure,fileName):
    #Save the a Structure in the filename
    w = FcfrpPDBBio(getAmountModel(structure))
    dir,file = os.path.split(fileName)
    w.save(structure, dir, file)

def loadStructure(id,fileName,PERMISSIVE=True):
    #Give id and fileName, return a biopython Structure
    #parser = PDBParser(PERMISSIVE)
    parser = FcfrpPDBParser(PERMISSIVE,LayoutDatabase=False)
    #return parser.get_structure(id, fileName)
    return parser.loadStructureFromFile(id, fileName) 
    
def checkPDBFile(id,fileName):
    #True structure is good. False, other option
    #return True
    parser = PDBParser(PERMISSIVE=True)
    structure = parser.get_structure(id, fileName)
    return checkStructure(structure,fileName)

def checkStructure(structure,fileName):
    result_CodError = []
    i = 0
    vs = FcfrpValidation(structure, fileName)
    vs.checkStructure()
    result_CodError = vs.hasErrors()    
    return result_CodError, vs

def removeHeteroAtoms(structure):
    #Remove whole hetero atoms from the Structure
    for model in structure:
        for chain in model:
            for residue in chain:
                id = residue.id
                if id[0] != ' ':
                    chain.detach_child(id)
            if len(chain) == 0:
                model.detach_child(chain.id)             
    return structure
        
def getPDBFromSite(id): 
    fullcgi = "http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&compression=None&pdbId=%s" % (id,id)
    handle = urllib.urlopen(fullcgi)
    uhandle = File.UndoHandle(handle)
    if not uhandle.peekline():
        raise IOError, "no results"
    return uhandle
    
def savePDBFile(pdb_file,pathFilename):
    #Generate a Pdb File from a File
    pdbfile2save = pdb_file.read()
    file2save = open(pathFilename, "w") 
    file2save.write(pdbfile2save) 
    file2save.close()

def isPDBFile(pathFileName):
    #This method check if the file is a PDB file or not. if the file contains ATOM, REMARK and SEQRES, this file is a PDB. Otherwise,
    # this file is not a PDB file.
    path,name = os.path.split(pathFileName)    
    fcfrpFile = FcfrpFile(path,name)
    #if fcfrpFile.find("ATOM") == -1 and fcfrpFile.find("REMARK")  == -1 and fcfrpFile.find("SEQRES") and -1:
    if fcfrpFile.find("ATOM") == -1 and fcfrpFile.find("SEQRES") == -1:
        mensage = "The %s file is not a PDB File. Please, check it." % pathFileName 
        raise Exception(mensage)

        
def getStructure(id,modelChoose=None,path):
    # This function returns a Structure that will be from PDB file.
    # If there isn't PDB file, it will be obtained from PDB site

    F = FcfrpFile(path,pdbid)
    if not F.existsFile():
        pdbFile = getPDBFromSite(id)
    savePDBFile(pdbFile, F._path)
    isPDBFile(F._path)        
    errors, vs = checkPDBFile(id, F.getPath())
    if len(errors) == 0:
       structure = loadStructure(id, F.getPath(), True)
       return structure
    else:
       s = FcfrpShowErrosDetails(vs, errors)
       s.showErrors()
       return loadStructure(id, F.getPath(), True)

def getStructureChainsFromFile(id,filename):
    #This method receives a PDBid and split it into its chains which will be stored in the dicStructures. It dictionary 
    # will be returned
    dicStructures = {}
    structure = loadStructure(id,filename)
    seqRes = structure.get_SeqRes()
    for chainId in seqRes.getAllChain():
        parser = FcfrpPDBParser(PERMISSIVE=True,LayoutDatabase=False)
        dicStructures[chainId] = parser.loadStructureFromFile(id, filename, chainId)
    return dicStructures