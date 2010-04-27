class FcfrpStandardizationPDBErrorDetails:
    
    def standardization(self,error,dicErrors):
        if error == 1:
            return self.detailsForMissingResidues(dicErrors)
        elif error == 2:
            return self.detailsForMissAtom(dicErrors)
        elif error == 3:
            return self.detailsForDupliAtom(dicErrors)
        elif error == 4:
            return self.detailsForUnkResidue(dicErrors)
    
    def formatMessageError(self,error,message):
        if error == 1:
            return self.formatMessageForMissingResidues(message)
        elif error == 2:
            return self.formatMessageForMissAtom(message)        
        elif error == 3:
            return self.formatMessageForDupliAtom(message)
        elif error == 4:
            return self.formatMessageForUnkResidue(message)

    def formatMessageForMissingResidues(self,message):
        chain, amino, position = str(message).split()
        m = "The amino %s on position %s and chain %s is missing " % (amino, position, chain) 
        return m
    
    def formatMessageForMissAtom(self, message):
        chain, number, amino, amount, atom = str(message).split()
        m = "The amino %s on position %s of chain %s has %s %s atoms missing " % (amino, number, chain, amount, atom)
        return m
    
    def formatMessageForDupliAtom(self,message):
        return message
    
    def formatMessageForUnkResidue(self,message):
        return message
        
    def detailsForMissingResidues(self,detail):
        dicRes = {}
        seq = 0
        for key in detail:
            listResult = []
            listResult = detail[key]
            indice = 0
            while indice < listResult.__len__():
                element = listResult
                seq += 1
                message = str(key) + " " + str(element[indice]) + " " + str(element[indice + 1])
                indice = indice + 2
                dicRes[seq] = message
        return dicRes
    
    def detailsForMissAtom(self,detail):
        dicRes = {}
        seq = 0
        for key in detail:           
            listResult = detail[key]
            indice = 0
            while indice < listResult.__len__():
                element = listResult
                seq += 1
                message = str(str(key[0]) + " " + str(key[1]) + " " + str(key[2]) + " " + str(element[indice]) + " " + str(element[indice + 1]))
                indice = indice + 2
                dicRes[seq] = message
        return dicRes
    
    def detailsForDupliAtom(self,detail):
        dicRes = {}
        seq = 0
        for key in detail:           
            listResult = detail[key]
            element = listResult
            indice = 0
            while indice < listResult.__len__():
                message = ''
                message = str(str(key[0]) + " " + str(key[1]) + " " + str(key[2]) + " " + str(key[3]))
                seq += 1
                message = message + " " + str(element[indice] + " " + str(element[indice + 1]))
                indice += 2
                dicRes[seq] = message
        return dicRes
                

    def detailsForUnkResidue(self,detail):
        dicRes = {}
        seq = 0
        indice = 0
        element = detail
        while indice < detail.__len__():
            seq += 1
            message = str(element[indice]) + " " + str(element[indice + 1])
            indice = indice + 2
            dicRes[seq] = message
        return dicRes
    