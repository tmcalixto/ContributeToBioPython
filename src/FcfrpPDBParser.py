#from Bio.PDB.StructureBuilder import StructureBuilder
import os
from FcfrpStructureBuilder import FcfrpStructureBuilder
from FcfrpDatabase import FcfrpDatabase
from FcfrpPDBDatabase import FcfrpPDBDatabase
from FcfrpFile import FcfrpFile
from FcfrpSeqRes import FcfrpSeqRes
from FcfrpResidue import FcfrpResidue
from FcfrpPDBAtom import FcfrpPDBAtom 
from FcfrpSSBonds import FcfrpSSBonds
from FcfrpValidation import FcfrpValidation
from FcfrpStructureError import FcfrpStructureError
from FcfrpStandardizationPDBErrorDetails import FcfrpStandardizationPDBErrorDetails

class PDBFileLayout:
    #This class works with PDB layout. It is used to load the file in a dictionary
    def __init__(self):
        #Store the chain has already ready. Because, The FcfrpPDBSeqRes class works with a SeqRes object by chain. This reason is
        # because with the unique SeqRes object by chain, we can obtain the whole information. And another reason is because the dictionary
        # which is read by FcfrpPDBSeqRes class needs a integer value in its key.
        #This dictionary has your key with chainid and your value is the key for self._dSeqRes    
        self._SeqResChain = {}
        #Store the SeqRes object from each chain
        self._dSeqRes = {}
        #Store the residue names only 
        self._dResName = {}
        self.chainIdRef = ""
        self._Seq = 0
        #Store whole biopython atom class 
        self._dAtom = {} 
        #Store whole SSBonds class
        self._dSSBonds = {}
        #Store The Remark
        self._dRemark = {}
        #Store structure errors. key = error and values is other dictionary
        self._dicErrors = {}
            
    def getSeqRes(self):
        return self._dSeqRes
    
    def getRes(self):
        return self._dResName
             
    def setSeqRes(self, record):
        seqRes = FcfrpSeqRes()
        seqRes.setSerNum(record[1])
        seqRes.setChainId(record[2])
        seqRes.setNumRes(record[3])
        if not self._SeqResChain.has_key(seqRes.getChainId()):
            ind = self._dSeqRes.__len__()
            self._SeqResChain[seqRes.getChainId()] = ind  
        i = self._SeqResChain[seqRes.getChainId()]
        self._dSeqRes[i] = seqRes

    def setRes(self, record):
        i = int(4)
        while i < len(record):
            r = str(record[i])
            res = FcfrpResidue(idAmino=None, seq=self._Seq, name=r, fcfrpDatabase=None,chainId=record[2])
            ind = self._dResName.__len__()
            self._dResName[ind] = res
            i += 1

    def setAtom(self, currentModel,line):
        #Here is read the line and store it in dAtom which will be read by getAtoms method
        fullname = line[12:16]
        name = str(fullname).strip()
        altloc = line[16:17]
        resname = line[17:20] 
        chainid = line[21:22] 
        if str(line[6:11]).strip() == "":
            serial_number = 0
        else:
            serial_number = int(line[6:11])
        resseq = line[22:26]
        icode = line[26:27]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occupancy = float(line[54:60])
        bfactor = float(line[60:66])
        segid = line[72:76]
        self._dAtom[self._dAtom.__len__()] = FcfrpPDBAtom(None, currentModel, name, x, y, z, bfactor, occupancy, altloc, fullname, serial_number, 0, resseq, icode, chainid, segid,resname)
        
    def getAtoms(self):
        return self._dAtom

    def setSSBonds(self,strLine):
        # Todos os indices iniciais esta configurados como n-1
        # onde n e a posicao indicada na documentacao do pdb
        # http://www.wwpdb.org/documentation/format32/sect6.html
        # campo      = strLine[0:6]
        sernum     = strLine[7:10]
        resi1      = strLine[10:14]
        chainid1   = strLine[14:16]
        seqnum1    = strLine[17:21]
        icode1     = strLine[21:22]
        resi2      = strLine[25:28]
        chainid2   = strLine[29:30]
        seqnum2    = strLine[31:35]
        icode2     = strLine[35:36]
        sym1       = strLine[59:65]
        sym2       = strLine[66:72]
        sym2       = strLine[66:72]
        length     = strLine[73:78]
        
        # print 'campo ' + campo
        # print 'sernum ' + sernum
        # print 'resi1 ' + resi1 
        # print 'chainid1 ' + chainid1 
        # print 'seqnum1 ' + seqnum1 
        # print 'icode1 ' + icode1 
        # print 'resi2 ' + resi2 
        # print 'chainid2 ' + chainid2 
        # print 'seqnum2 ' + seqnum2
        # print 'icode2 ' +  icode2
        # print 'sym1 ' +  sym1
        # print 'sym2 ' +  sym2
        # print 'length ' +  length
        # print '\n\n'        
        
        ssbond = FcfrpSSBonds()
        ssbond.setSerNum(sernum)
        ssbond.setChainId1(chainid1)
        ssbond.setSeqNum1(seqnum1)
        ssbond.setIcode1(icode1)
        ssbond.setChainId2(chainid2)
        ssbond.setSeqNum2(seqnum2)
        ssbond.setIcode2(icode2)
        ssbond.setSym1(sym1)
        ssbond.setSym2(sym2)
        ssbond.setLength(length)
        self._dSSBonds[self._dSSBonds.__len__()] = ssbond 
        
    def getSSBonds(self):
        return self._dSSBonds
    
    def setStructureErrors(self,errors,vs, standard = None):
        if standard == None:
            standard = FcfrpStandardizationPDBErrorDetails()
        dic = {}
        for e in errors:
            if e == 1:
                dic = standard.standardization(e,vs.getMissingResidues())
                self._setError(dic, e)
            elif e == 2:
                dic = standard.standardization(e,vs.getMissingAtoms())
                self._setError(dic, e)
            elif e == 3:
                dic = standard.standardization(e,vs.getDuplicatedAtoms())
                self._setError(dic, e)
            elif e == 4:
                dic = standard.standardization(e,vs.getUnknownResidues())
                self._setError(dic, e)
                
    def getStructureErrors(self):
        return self._dicErrors
    
    def _setError(self,dic,error):
        i = self._dicErrors.__len__()
        for seq, message in dic.iteritems():
            s = FcfrpStructureError()
            s.setIdCode( 0 )
            s.setIdError( error )
            s.setSeq( seq )
            s.setMessage( message )
            i = i + 1
            self._dicErrors[i] = s
                          
class FcfrpPDBParser():
    
    def __init__(self, PERMISSIVE=1, get_header=0, fcfrpDatabase=None, LayoutDatabase=False):
        self._PERMISSIVE = PERMISSIVE
        #If LayoutDatabase, I'll load the classes about database. Otherwise, I'll not load them.
        if LayoutDatabase == True:
            if fcfrpDatabase == None:
                self._fcfrpDatabase = FcfrpDatabase()
            else:
                self._fcfrpDatabase = fcfrpDatabase
            self._pdbDatabase = FcfrpPDBDatabase()
        else:
            self._pdbFileLayout = PDBFileLayout()
            self._dicResidue = {} 
            self._dicSeqRes = {}
                
        #self._structure_builder=StructureBuilder()
        self._structure_builder = FcfrpStructureBuilder()        
        
    def _loadAtomsfromDatabase(self, idcode, modelChoose, chainId = None):
        #Load the Dictionaries
        D_Atoms = self._pdbDatabase.readAtomsPDBFromResult(self._fcfrpDatabase.getAtomsPDB(idcode, modelChoose, chainId))
        i_a = 0 #atom index
        #Create reference variables
        reference_segid = None
        reference_chain = None
        reference_residue = None
        reference_resName = None
        reference_model = None
        inserted_model0 = False
        
        while (i_a < len(D_Atoms)):
            resname = str(self._fcfrpDatabase.getFullInformationAminoacid(D_Atoms[i_a].getIdAmino())[3])
            resseq = D_Atoms[i_a].getResSeq()
            icode = D_Atoms[i_a].getIcode()
            residueId = ("", resseq, icode)    
            #Model
            if reference_model != D_Atoms[i_a].getModel() or D_Atoms[i_a].getModel() == None:
                if D_Atoms[i_a].getModel() == None: #Cristalografy
                    reference_model = 0
                    if inserted_model0 == False: 
                        self._structure_builder.init_model(reference_model)
                    inserted_model0 = True    
                else:    
                    reference_model = D_Atoms[i_a].getModel()
                    self._structure_builder.init_model(reference_model)
                    reference_chain = None
                    reference_residue = None
                    reference_resName = None
            #Segid
            #if reference_segid != dicAtom[i_a].getSeqId():
            #    reference_segid = dicAtom[i_a].getSeqId()
            self._structure_builder.init_seg("")
            reference_segid = i_a
            #Chain and Residue
            if reference_chain != D_Atoms[i_a].getChainId():
                reference_chain = D_Atoms[i_a].getChainId()
                self._structure_builder.init_chain(reference_chain)
                #Because changed chain, it's necessary start a new residue in structure
                reference_residue = residueId
                reference_resName = resname
                self._structure_builder.init_residue(resname, "", resseq, icode)
            elif reference_resName != resname or reference_residue != residueId:
                reference_resName = resname
                reference_residue = residueId
                self._structure_builder.init_residue(resname, "", resseq, icode)
            #Atom
            #print (dicAtom[i_a].get_id(), dicAtom[i_a].get_coord(),  dicAtom[i_a].get_fullname(), dicAtom[i_a].get_serial_number(),dicAtom[i_a].getChainId(),resname,resseq)    
            self._structure_builder.init_atom(D_Atoms[i_a].get_id(), D_Atoms[i_a].get_coord(), D_Atoms[i_a].get_bfactor(), D_Atoms[i_a].get_occupancy() , D_Atoms[i_a].get_altloc(), D_Atoms[i_a].get_fullname(), D_Atoms[i_a].get_serial_number())            
            #self._structure_builder.init_atom(D_Atoms[i_a].get_id(), reference_model, reference_resName, D_Atoms[i_a].get_coord()[0], D_Atoms[i_a].get_coord()[1], D_Atoms[i_a].get_coord()[2], D_Atoms[i_a].get_bfactor(), D_Atoms[i_a].get_occupancy() , D_Atoms[i_a].get_altloc(), D_Atoms[i_a].get_fullname(), D_Atoms[i_a].get_serial_number(), reference_residue, resseq, reference_chain, reference_segid, reference_resName )
            i_a = i_a + 1
    
    def _loadSeqResfromDatabase(self, idcode, id, chainId = None):
        dicResidue = self._pdbDatabase.readResiduesFromResult(self._fcfrpDatabase.getResidues(idcode))
        dicSeqRes = self._pdbDatabase.readSeqResFromResult(self._fcfrpDatabase.getSeqRes(idcode,chainId))
        self._structure_builder.init_SeqRes(id)
        self._structure_builder.structure.add_SeqRes(dicSeqRes, dicResidue)
    
    def _loadSSBondsfromDatabase(self,idcode,id, chainId):
        dicSSBonds = self._pdbDatabase.readSSBondsPDBFromResult(self._fcfrpDatabase.getSSBonds(idcode, chainId)) 
        self._structure_builder.init_SSBonds(id)
        self._structure_builder.structure.add_SSBonds(dicSSBonds)
    
    def _loadErrorsFromDatabase(self,id):
        dicErrors = self._pdbDatabase.readPDBErrorsFromResult(self._fcfrpDatabase.getPDBErrors(id))
        self._structure_builder.init_ErrorStructure(id)
        self._structure_builder.structure.add_StructureErrors(dicErrors)
        
    def loadStructureFromDatabase(self, id, modelChoose=None,chainId=None):
        #Give an id, load of Structure from Database
        #if the Structure is a NRM, it's necessary inform its model which want to get
        idcode = self._fcfrpDatabase.getIdCode(id)
        #Start the Structure
        self._structure_builder.init_structure(id)
        #SeqRes
        self._loadSeqResfromDatabase(idcode, id, chainId)
        #SSbonds
        self._loadSSBondsfromDatabase(idcode, id, chainId)
        #Atoms
        self._loadAtomsfromDatabase(idcode, modelChoose, chainId)
        #Errors
        self._loadErrorsFromDatabase(id)
        return self._structure_builder.get_structure()
    
    def _loadAtomsfromFile(self, id, dicAtom):
        i_a = 0 #atom index
        #Create reference variables
        reference_chain = None
        reference_residue = None
        reference_resName = None
        reference_model = None
        inserted_model0 = False
        while i_a < dicAtom.__len__():
            resname = dicAtom[i_a].getResName()
            resseq = dicAtom[i_a].getResSeq()
            icode = dicAtom[i_a].getIcode()
            residueId = ("", resseq, icode)    
            #Model
            if reference_model != dicAtom[i_a].getModel() or dicAtom[i_a].getModel() == None:
                if dicAtom[i_a].getModel() == None: #Cristalografy
                    reference_model = 0
                    if inserted_model0 == False: 
                        self._structure_builder.init_model(reference_model)
                    inserted_model0 = True    
                else:    
                    reference_model = dicAtom[i_a].getModel()
                    self._structure_builder.init_model(reference_model)
                    reference_chain = None
                    reference_residue = None
                    reference_resName = None
            #Segid
            self._structure_builder.init_seg("")
            #Chain and Residue
            if reference_chain != dicAtom[i_a].getChainId():
                reference_chain = dicAtom[i_a].getChainId()
                self._structure_builder.init_chain(reference_chain)
                #Because changed chain, it's necessary start a new residue in structure
                reference_residue = residueId
                reference_resName = resname
                self._structure_builder.init_residue(resname, "", resseq, icode)
            elif reference_resName != resname or reference_residue != residueId:
                reference_resName = resname
                reference_residue = residueId
                self._structure_builder.init_residue(resname, "", resseq, icode)
            #Atom    
            self._structure_builder.init_atom(dicAtom[i_a].get_id(), dicAtom[i_a].get_coord(), dicAtom[i_a].get_bfactor(), dicAtom[i_a].get_occupancy() , dicAtom[i_a].get_altloc(), dicAtom[i_a].get_fullname(), dicAtom[i_a].get_serial_number())
            i_a = i_a + 1
        
    def _loadSeqResfromFile(self, id):
        self._structure_builder.init_SeqRes(id)
        self._structure_builder.structure.add_SeqRes(self._pdbFileLayout.getSeqRes(), self._pdbFileLayout.getRes())
        
    def _loadSSBondsfromfile(self,id):
        self._structure_builder.init_SSBonds(id)
        self._structure_builder.structure.add_SSBonds(self._pdbFileLayout.getSSBonds())

    def _loadErrorsfromFile(self, id, structure, pathfileName, chainId ):
        self._structure_builder.init_ErrorStructure(id)
        errors, vs = hasErrorFileStructure(structure,pathfileName , chainId)
        if len(errors) > 0:
            standard = FcfrpStandardizationPDBErrorDetails()
            self._pdbFileLayout.setStructureErrors(errors, vs, standard)
            self._structure_builder.structure.add_StructureErrors(self._pdbFileLayout.getStructureErrors()) 

    def loadStructureFromFile(self, id, pathfileName, chainId=None):
        dir, fileName = os.path.split(pathfileName)
        filePDB = FcfrpFile(dir, fileName)
        #Define Model value
        currentModel = None
        #Start the Structure
        self._structure_builder.init_structure(id)
        #Here the PDB file will read and separate in dictionaries. 
        for line in filePDB.readLines():
            record = line.split()
            if record[0] == "SEQRES":
                if self.checkChain(record[2], chainId) == True:
                    self._pdbFileLayout.setRes(record)
                    self._pdbFileLayout.setSeqRes(record)
            elif record[0] == "SSBOND":
                if self.checkChain(line[15:16], chainId) == True:
                    self._pdbFileLayout.setSSBonds(line)
            elif record[0] == "MODEL":
                #The MODEL line is MODEL        1                                                                  
                # So, when split command is used, the result'll be ['MODEL', '1']
                #Therefore, the currentModel represents the Model. If there isn't model, its value'll be always zero.
                currentModel = int(record[1]) 
            elif record[0] == "ATOM":
                if self.checkChain(line[21:22], chainId) == True:
                    self._pdbFileLayout.setAtom(currentModel,line)
        #Here the dictionaries have been created by for command above, they will be used in specific methods. These methods
        # represent each PDB file section.
        self._loadSeqResfromFile(id)
        self._loadSSBondsfromfile(id)
        self._loadAtomsfromFile(id, self._pdbFileLayout.getAtoms())
        structure = self._structure_builder.get_structure()
        
        self._loadErrorsfromFile(id, structure, pathfileName, chainId)
        #print "Comentado checagem de erro qdo arquvivo PDB"
        
        return structure
    
    def checkChain(self, chainIdFile, chainIdSelected):
        cont = False
        if chainIdSelected == None:
            cont = True
        elif chainIdFile == chainIdSelected:
            cont = True
        else:
            cont = False
        return cont 

    
def hasErrorFileStructure(structure,pathfileName, chainId):
    result_CodError = []
    i = 0
    vs = FcfrpValidation(structure, pathfileName, chainId)
    vs.checkStructure()
    result_CodError = vs.hasErrors()    
    return result_CodError, vs

