from Bio.PDB.Entity import Entity

class FcfrpPDBSSBonds(Entity):
    def __init__(self,id):
        Entity.__init__(self, id)
        #Store the SSBONDS protein.
        self._SSBonds= {}
    
    def addSSBonds(self,dic):
        self._SSBonds = dic
        
    def getAllSSBonds(self):
        return self._SSBonds
    
    def getamountSSBOnds(self):
        return self._SSBonds.__len__()    
    
    def getLineFullInformation(self,line):
        if int(line) > self.getamountSSBOnds():
            raise Exception("The line number is more than SSBonds number")
        return self._SSBonds[line].getFullInformation()
    