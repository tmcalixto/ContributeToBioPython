from Bio.PDB.Entity import Entity
from FcfrpSeqRes import FcfrpSeqRes
from FcfrpFasta import FastaChainSequence
    
#This class works with SEQRES and Residue objects. Because these objects are divide in Database.
# So, this class join them and it'll store in dictionary.  
class PDBSeqRes():
    def __init__(self, R=None,S=None):
        self._SeqRes = R
        self._Residue = S
        
    def setSeqRes(self,S):
        self._SeqRes = S
        
    def setResidue(self,R):
        self._Residue = R
    
    def getSeqRes(self):
        return self._SeqRes
    
    def getResidue(self):
        return self._Residue 

    def getChainId(self):
        return self._SeqRes.getChainId()

    def getIdl3(self):
        return self._Residue.getIdl3()
    
    def getSerNum(self):
        return self._SeqRes.getSerNum()
    
    def getNumRes(self):
        return self._SeqRes.getNumRes
                
    def getFullInformation(self):
        return (self._SeqRes.getSerNum(), self._SeqRes.getChainId(), self._SeqRes.getNumRes(), self._Residue.getIdl3())
    
class FcfrpPDBSeqRes(Entity):
    
    def __init__(self,id):
        Entity.__init__(self, id)
        #Store the SEQRES protein. Whole its SEQRES will be record it
        self._SeqRes={}
        self._iSeq = -1
        #Store chain and its residue number. Because this value is very used. Thus, this dictionary help to get it
        self._ChainNumRes = {}
            
    def setSeqResFromDic(self,dicSeqRes,dicResidue):
        # These dictionaries represent the database Tables t_pdb_seqres and t_pdb_seqres_resname, respectively.
        # Therefore, together these dictionaries represent the SEQRES information.
        
        #Obtain the chains of protein
        DicChains = {}
        i = 0
        #After obtained the protein chains, we must have to keep 
        # the chain position of original PDB. Therefore, we create   
        # a dictionary where its key is sortInserted
        sortInserted = 0
        while i < len(dicSeqRes):
            S = dicSeqRes[i]
            if DicChains.has_key(S.getChainId()) == False:
                DicChains[sortInserted] = [ str(S.getChainId()).strip() , S]
                sortInserted += 1
            i += 1
        dicSorted = DicChains.keys()
        dicSorted.sort()
        for sortInserted in dicSorted:
             listObj = DicChains[sortInserted]
             chain,seqRes = listObj 
             self.setChainNumRes(chain, seqRes.getNumRes())
             #Get Residues from chain
             i = 0
             while i < len(dicResidue):
                res = dicResidue[i]
                if (res.getChainId() == chain):
                    self.setSeqRes(seqRes, res)
                i += 1
            
    def setSeqRes(self,SeqRes,Residue):
        pdbSeqRes = PDBSeqRes(SeqRes,Residue)
        self._iSeq += 1
        self._SeqRes[self._iSeq] = pdbSeqRes
    
    def setChainNumRes(self,chainId,NumRes):
        self._ChainNumRes[chainId] = NumRes
            
    def getResChain(self,chainId):
        resChain = []
        i = 0
        while i < len(self._SeqRes):
            S = self._SeqRes[i]
            if (S.getChainId() == chainId) or (chainId == None):
                resChain.append(S)
            i += 1
        return resChain #PDBSeqRes Object 
    
    def getNumResChain(self,chainId):
        return self._ChainNumRes[chainId]
        
    def getResName(self,chainId,position):
        resChain = self.getResChain(chainId)
        if position > len(resChain):
            raise Exception("Position is more than residues number of chain %s" % chainId)  
        return resChain[position].getIdl3()
    
    def getAllRes(self):
        return self.getResChain(None)         
    
    def getAllChain(self):
        #Return a list of all chains which are present in protein.
        chain = []
        i = 0
        while i < len(self._SeqRes):
            S = self._SeqRes[i]
            if S.getChainId() not in chain:
                chain.append(S.getChainId())
            i += 1
        return chain
    
    def getResChainAsListSequence(self,chainId):
        #Return the protein chain residues as list
        listAux = []
        seqRes = self.getResChain(chainId)
        for r in seqRes:
            listAux.append(r.getIdl3())
        return listAux
            