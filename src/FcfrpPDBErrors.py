from Bio.PDB.Entity import Entity

class FcfrpPDBErrors(Entity):
    def __init__(self,id):
        Entity.__init__(self, id)
        #Store the structure errors.
        self._dicStructureErrors= {}

    def addStructureErrors(self,dic):
        self._dicStructureErrors = dic
        
    def getStructureErrors(self):
        return self._dicStructureErrors
    
    def getamountStructureErrors(self):
        return self._dicStructureErrors.__len__()    
    
    def getLineFullInformation(self,line):
        if int(line) > self.getamountStructureErrors():
            raise Exception("The line number is more than Error number")
        return self._dicStructureErrors[line].getFullInformation()
