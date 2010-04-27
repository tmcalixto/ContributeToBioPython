class FcfrpShowErrosDetails:
    
    def __init__(self, s, errors):
        self.vs = s
        self.errors = errors
        
    def showErrors(self):
        for error in self.errors:
            # check for missing residues
            if error == 1:
                print "\nThere are missing residues"
                miss = self.vs.getMissingResidues()
                for residues in miss:
                    print residues
                    print miss[residues]
            # check for missing atoms        
            if error == 2:
                print "\nThere are missing atoms"
                miss = self.vs.getMissingAtoms()
                for atom in miss:
                    print atom, ": ", miss[atom]
            # check for duplicated atoms
            if error == 3:
                print "\nThere are duplicated atoms"
                dupli = self.vs.getDuplicatedAtoms()
                for atom in dupli:
                    print atom, ": ", dupli[atom]
            # check for unknown residues
            if error == 4:
                print "\nThere are unknown residues"
                unk = self.vs.getUnknownResidues()
                print unk