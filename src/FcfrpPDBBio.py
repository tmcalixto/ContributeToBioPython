from FcfrpFile import FcfrpFile
from FcfrpStandardizationPDBErrorDetails import FcfrpStandardizationPDBErrorDetails

class FcfrpPDBBio():
    def __init__(self,amountModels=0):
        #Based on PDBBio from biopython project
        self._ATOM_FORMAT="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"
        self._SEQRES_FORMAT="%s  %2i %s %4i  %s \n"
        #self._SSBOND_FORMAT = "%s  %2i %s %s   %3i   %s %s %s %4i                            %s %s %s\n"
        self._SSBOND_FORMAT = "%6s  %2i %3s%1s  %3i%3s %1s %1s %4i %1s                      %4s %4s %s\n"
        self._REMARK_FORMAT = "%s  1 %s \n"
        self._amountModels = amountModels
        #Default layouts
        self._MODEL = "Model %s\n"
        self._ATOM = "ATOM  "
        self._TER = "TER                                                                      \n"
        self._ENDMDL = "ENDMDL                                                               \n"
        self._SEQRES = "SEQRES"
        self._SSBOND = "SSBOND"
        self._REMARK = "REMARK"
        
    def _saveSeqRes(self,structure,filePDB):
        seqRes = structure.get_SeqRes()
        for chain in seqRes.getAllChain():
            residuesChain = seqRes.getResChain(chain)
            full = residuesChain[0].getFullInformation()
            serNum  = 0#full[0]
            chainid = full[1]
            NumRes  = full[2]
            i = 0
            residues = ""
            iRes = 0
            while iRes < len(residuesChain):
                residues = residues + seqRes.getResName(chain,iRes) + " "
                i += 1
                if i == 13:
                    serNum  += 1
                    values = (self._SEQRES,serNum,chainid,NumRes,residues)
                    i = 0
                    filePDB.write(self._SEQRES_FORMAT % values)
                    residues = ""
                iRes += 1
            if residues != "":
                serNum  += 1
                values = (self._SEQRES,serNum,chainid,NumRes,residues)
                filePDB.write(self._SEQRES_FORMAT % values)
                residues = "" 
        
    def _saveAtom(self,iAtom,atom,residue,chain,iRes,filePDB, resID):
        x, y, z=atom.get_coord()
        occupancy  = atom.get_occupancy()
        bfactor = atom.get_bfactor()
        segid = residue.get_segid()        
        icode = residue.id[2]
        #values = (self._ATOM,iAtom,' '+atom.get_fullname(),atom.get_altloc(),residue.get_resname(),chain.id,resID,icode, x,y,z,
        #                      occupancy,bfactor,segid," "," ")
        values = (self._ATOM,iAtom,atom.get_fullname(),atom.get_altloc(),residue.get_resname(),chain.id,resID,icode, x,y,z,
                              occupancy,bfactor,segid," "," ") 
        filePDB.write(self._ATOM_FORMAT % values)
        
    def _saveSSBonds(self,structure,filePDB):
        ssbonds = structure.get_SSBonds()
        i = 0
        while i < ssbonds.getamountSSBOnds():
            l = ssbonds.getLineFullInformation(i)
            sernum     = l[2]
            resi1      = "CYS"
            chainid1   = l[3]
            seqnum1    = l[4]
            icode1     = l[5]
            resi2      = "CYS"
            chainid2   = l[6]
            seqnum2    = l[7]
            icode2     = l[8]
            sym1       = l[9]
            sym2       = l[10]
            length     = l[11]
            
            if str(length).__contains__('None'):
                length = ''
                
            values = (self._SSBOND, int(sernum), resi1, chainid1, int(seqnum1), icode1, resi2, chainid2, int(seqnum2), icode2, sym1, sym2, length)
            filePDB.write(self._SSBOND_FORMAT % values)
            i += 1  
        
        
    def _saveRemark(self, structure, filePDB):
        #Structure Errors
        standard = FcfrpStandardizationPDBErrorDetails()
        for k,v in structure.get_StructureErrors().getStructureErrors().iteritems():
            m = standard.formatMessageError(v.getIdError(), v.getMessage())
            values = ( self._REMARK, m )
            filePDB.write(self._REMARK_FORMAT % values)

    def save(self,structure,dir,fileName):
        filePDB = FcfrpFile(dir, fileName)
        filePDB.getAsFile("w")
        #Create Remark Section
        self._saveRemark(structure, filePDB)
        #Create SEQRES Section
        self._saveSeqRes(structure, filePDB)
        #Create SSBOND Section
        self._saveSSBonds(structure, filePDB)
        #Create ATOM Section
        iModels = 0 #index for Models
        iAtom = 0 #index for atom
        iRes = 0 #index for Residue
        for model in structure:
            #Will write Model, if there are more than 1 models.
            if self._amountModels > 1:
                iModels += 1
                filePDB.write(self._MODEL % str(iModels)) 
            iAtom = 0
            iRes  = 0
            for chain in model:
                for residue in chain:
                    iRes  += 1
                    for atom in residue: #Atoms
                        iAtom += 1
                        self._saveAtom(iAtom,atom,residue,chain,iRes,filePDB, int(residue.id[1]))
            filePDB.write(self._TER)
            filePDB.write(self._ENDMDL)
        filePDB.close()