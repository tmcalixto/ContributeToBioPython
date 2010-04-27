class FcfrpDescricaoAmino:
    
    dicionario = {}
    
    dicALA = {'C':'3', 'O':'2', 'N':'1', 'H':'7'}
    dicVAL = {'C':'5', 'O':'2', 'N':'1', 'H':'11'}
    dicLEU = {'C':'6', 'O':'2', 'N':'1', 'H':'13'}
    dicILE = {'C':'6', 'O':'2', 'N':'1', 'H':'13'}
    dicPRO = {'C':'5', 'O':'2', 'N':'1', 'H':'8'}
    dicMET = {'C':'5', 'O':'2', 'N':'1', 'H':'11', 'S':'1'}
    dicPHE = {'C':'9', 'O':'2', 'N':'1', 'H':'11'}
    dicTRP = {'C':'11', 'O':'2', 'N':'2', 'H':'13'}
    dicGLY = {'C':'2', 'O':'2', 'N':'1', 'H':'5'}
    dicSER = {'C':'3', 'O':'3', 'N':'1', 'H':'7'}
    dicTYR = {'C':'9', 'O':'3', 'N':'1', 'H':'11'}
    dicGLN = {'C':'5', 'O':'3', 'N':'2', 'H':'11'}
    dicASN = {'C':'4', 'O':'3', 'N':'2', 'H':'8'}
    dicCYS = {'C':'3', 'O':'2', 'N':'1', 'H':'7', 'S':'1'}
    dicTHR = {'C':'4', 'O':'3', 'N':'1', 'H':'9'}
    dicASP = {'C':'4', 'O':'4', 'N':'1', 'H':'6'}
    dicGLU = {'C':'5', 'O':'4', 'N':'1', 'H':'8'}
    dicLYS = {'C':'6', 'O':'2', 'N':'2', 'H':'15'}
    dicARG = {'C':'6', 'O':'2', 'N':'4', 'H':'15'}
    dicHIS = {'C':'6', 'O':'2', 'N':'3', 'H':'10'}
    
    def __init__(self):
        self.dicionario = {'ALA':self.dicALA, 'VAL':self.dicVAL, 'LEU':self.dicLEU, 'ILE':self.dicILE, 'PRO':self.dicPRO, 'MET':self.dicMET, 'PHE':self.dicPHE, 'TRP':self.dicTRP, 'GLY':self.dicGLY, 'SER':self.dicSER, 'TYR':self.dicTYR, 'GLN':self.dicGLN, 'ASN':self.dicASN, 'CYS':self.dicCYS, 'THR':self.dicTHR, 'ASP':self.dicASP, 'GLU':self.dicGLU, 'LYS':self.dicLYS, 'ARG':self.dicARG, 'HIS':self.dicHIS}
        
    def getDicionario(self, residue):
        return self.dicionario[str(residue)]