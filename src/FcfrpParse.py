from FcfrpFile import FcfrpFile
from FcfrpConfig import FcfrpConfig
from FcfrpExecuteProgram import FcfrpExecuteProgram
import os
import shutil
from FcfrpUtils import *

class FcfrpParse:
    def __init__(self, pqrName=None, pathPqr=None, pathOut=None):
        self._config = FcfrpConfig("GLOBAL")
        
        self._pqrName = pqrName
        
        #Path where are the input files
        if pathPqr == None:
            self._pathPqr = self._config.getParameter("Local_PQR_Files")
        else:
            self._pathPqr = pathPqr 
        
        if pathOut == None:
            self._pathOut = self._config.getParameter("Local_MEAD_Files")
        else:
            self._pathOut = pathOut
            
        #Path where are the MEAD Programs
        self._pathMEADPrograms = self._config.getParameter("MEAD_Program") 
        self._inicializated()     
    
    
    def _inicializated(self):             
        # Titrable residues to Sites file
        self.resIonic = ['ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'TYR']     
        #self.resIonic = ['ASP', 'GLU', 'LYS', 'ARG']
        self.positive = ['ARG', 'LYS', 'HIS', 'NT']
        self.negative = ['GLU', 'ASP', 'TYR', 'CT']     
        
        
     # Return pka's for individual aminoacids
    def getPkas (self, amino):
        if amino == 'LYS':
            return 10.4
        elif amino == 'ASP':
            return 4.0
        elif amino == 'GLU':
            return 4.4
        elif amino == 'CT':
            return 3.8
        elif amino == 'HIS':
            return 6.3
        elif amino == 'TYR':
            return 9.6
        elif amino == 'NT':
            return 7.5
        elif amino == 'ARG':
            return 12.0 
        
        
    # Create St files for all ionic residues and N and C Terminal    
    def createStFiles(self):
        sitesFile = open(os.path.join(self._pathPqr, self._pqrName[0:self._pqrName.__len__()-4]+'.sites'))
        st = []
        linhaCT =  0 # Reference line for the last residue : C-Terminal
        lineNT  = -1 # Reference line for the first residue: N-Terminal
        
        for line in sitesFile.readlines():
            # Build a list for ionic residues (not repeated) plus N and C-Terminal
            if not line.split()[1] in st:
                st.append(line.split()[1])
                amino  = line.split()[1]
                lineNT = line.split()[0]
                lineCT = line.split()[0]
                
                # It creates st files for each residue not repeated in sites file 
                # Residue as N-Terminal. pk = 7.5
                if amino[0:2] == 'NT': 
                    self.createStNT(amino, lineNT)
                
                # Residue as C-Terminal. pk = 3.8    
                elif amino[0:2] == 'CT':
                    self.createStCT(amino, lineCT)
                
                # Ionic Residue as N-Terminal. pk = pk_residue
                elif amino[3:5] == 'NT':
                    self.createStNT2(amino, lineNT)
                
                # Ionic Residue as C-Terminal. pk = pk_residue
                elif amino[3:5] == 'CT':
                    self.createStCT2(amino, lineCT)
                
        # Others residuos
        self.buildStsFIle()    
        
    # It buil N-Terminal st file        
    def createStNT(self, amino, lineNT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))  
        atoms  = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == 'ATOM': 
                if fields[4] == lineNT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]

                    if resName in self.resIonic:
                        if atom == 'N' or atom == 'H' or atom == 'H2' or atom == 'H3':
                            atoms[atom] = str(charge)
                    else:
                        atoms[atom] = str(charge)
                               
        self.mountStNT(amino, atoms)
        pqrFile.close()
    
    # It buil N-Terminal st file
    # N-term is charged positively
    def mountStNT(self, amino, atoms):
        s = ''
        pk = self.getPkas('NT')
        s = s + str(pk) + '\n'
        
        for a in atoms:
            s = s + amino[2:5] + '\t' + a + '\t' + atoms[a] + '\t' + '0.00' + '\n'
        
        self.saveFile(amino+'.st', s)
        
    # It buil N-Terminal st file
    # In this case the residue N-Terminal is an ionic residue    
    def createStNT2(self, amino, lineNT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))  
        atoms  = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == "ATOM":
                if fields[4] == lineNT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]
                    
                    if resName in self.resIonic:
                        if atom != 'N' and atom != 'H' and atom != 'H2' and atom != 'H3':
                            atoms[atom] = str(charge)
                               
        self.mountStNT2(amino, atoms)
        pqrFile.close()
        
    # It buil N-Terminal st file
    # N-term is charged positively       
    # In this case the residue N-Terminal is an ionic residue 
    def mountStNT2(self, amino, atoms):
        s = ''
        pk = self.getPkas(amino[0:3])
        s = s + str(pk) + '\n'
        
        for a in atoms:
            s = s + amino[0:3] + '\t' + a + '\t' + atoms[a] + '\t' + '0.00' + '\n'
        
        self.saveFile(amino+'.st', s)
    
        
    # It buil C-Terminal st file
    # C-term is charged negatively 
    def createStCT(self, amino, linhaCT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))
        atoms = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == "ATOM":
                if fields[4] == linhaCT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]
                
                    if resName in self.resIonic:
                        if atom == 'C' or atom == 'O' or atom == 'OH' or atom == 'OXT':
                            atoms[atom] = str(charge)
                    else:
                        atoms[atom] = str(charge)
        
        self.mountStCT(amino, atoms)
        pqrFile.close()
        
    # It buil C-Terminal st file
    def mountStCT(self, amino, atoms):
        s = ''
        pk = self.getPkas('CT')
        s = s + str(pk) + '\n'

        for a in atoms:
            s = s + amino[2:5] + '\t' + a + '\t' + '0.00' + '\t' + atoms[a] + '\n'
        
        self.saveFile(amino+'.st', s)
        
        
    # It buil C-Terminal st file
    # In this case the residue C-Terminal is an ionic residue   
    def createStCT2(self, amino, linhaCT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))
        atoms = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == "ATOM":
                if fields[4] == linhaCT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]
                
                    if resName in self.resIonic:
                        if atom != 'C' and atom != 'O' and atom != 'OH' and atom != 'OXT':
                            atoms[atom] = str(charge)
        
        self.mountStCT2(amino, atoms)
        pqrFile.close()
        
    # It buil C-Terminal st file
    # C-term is charged negatively       
    # In this case the residue C-Terminal is an ionic residue 
    def mountStCT2(self, amino, atoms):
        s = ''
        pk = self.getPkas(amino[0:3])
        s = s + str(pk) + '\n'

        for a in atoms:
            s = s + amino[0:3] + '\t' + a + '\t' + '0.00' + '\t' + atoms[a] + '\n'
        
        self.saveFile(amino+'.st', s)
        
    # It creates a site file    
    def createSitesFile(self, last):
        pqr_file = FcfrpFile(self._pathPqr, self._pqrName + '.renumbered')
        record = ''
        refLine = -1
        NTref = 1
        flag = 0
        lastResidue = last
                
        for line in pqr_file.readLines(): 
            fields = line.split()
            if fields[0] != "END_CHAIN":
                campo   = fields[0]
                numSeq  = fields[1]
                atom    = fields[2]
                resName = fields[3]
                numRes  = fields[4]
                px      = fields[5]
                py      = fields[6]
                pz      = fields[7]
                charge  = fields[8]
                radii   = fields[9]
                
                # It used to insert N Terminal
                if flag == 1:
                    NTref = int(numSeq)
                    flag = 0
                
                if numRes != refLine:
                    if int(numSeq) == NTref:
                        # Get the all (the first and others) N-Terminal
                        record = self.insertNterminal(numRes, resName, record)
                
                    elif resName in self.resIonic:   
                        # It does not insert the last ionic residue   
                        if int(numSeq) == lastResidue:
                            break   
                        else:   
                            record = record + str(numRes) + "  " + str(resName) + '\n'
                    
                    refLine = numRes
            
            if fields[0] == "END_CHAIN":
                # Get the others C-Terminal            
                record = self.insertCterminal(numRes, resName, record)
                flag = 1 # It indicates that the next residue is a N Terminal

        # Get the last C-Terminal       
        record = self.insertCterminal(numRes, resName, record)
        
        self.saveFile(self._pqrName[:self._pqrName.__len__()-3] + 'sites', record)


    # Insert C terminal residue in sites file
    def insertCterminal(self, numRes, resName, record):
        record = record + numRes + "  " + 'CT' + resName + '\n'
        if resName in self.resIonic:
            record = record + numRes + "  " + resName + 'CT' + '\n'
        return record

    
    # Insert N terminal residue in sites file    
    def insertNterminal(self, numRes, resName, record):
        record = record + numRes + "  " + 'NT' + resName + '\n'
        if resName in self.resIonic:
            record = record + numRes + "  " + resName + 'NT' + '\n'
        return record
            
        
    # Check the correct charges in N and C-Terminal st files    
    def checaSt(self):
        sites = open(os.path.join(self._pathPqr, self._pqrName[0:self._pqrName.__len__()-4]+'.sites'))
        st = []
        for line in sites:
            tmp = line.split()
            residuo = tmp[1]
            
            if residuo not in st:
                st.append(residuo)
                if residuo[0:2] == 'NT':
                    self.checaNT(residuo + '.st')
                elif residuo[3:5] == 'NT':
                    self.checaSiteNT(residuo + '.st')
                elif residuo[0:2] == 'CT':
                    self.checaCT(residuo + '.st')
                elif residuo[3:5] == 'CT':
                    self.checaSiteCT(residuo + '.st')
                elif residuo in self.resIonic:
                    self.checaSites(residuo + '.st')
    
        
    # Check the correct charges in st's file  
    def checaSites(self, fileName):
        st = open(os.path.join(self._pathPqr, fileName))
        atomsPROT = {}
        atomsDEPROT = {}
        s = ''
        
        for line in st:
            tmp = line.split()
            if tmp.__len__() > 1:
                atomsPROT[tmp[1]] = tmp[2]
                atomsDEPROT[tmp[1]] = tmp[3]
        
        protonates = 0.0
        deprotonates = 0.0
        sum = 0.0
        
        for a in atomsPROT:
            protonates = protonates + float(atomsPROT[a])
        for a in atomsDEPROT:
            deprotonates = deprotonates + float(atomsDEPROT[a])
        
        sum = protonates - deprotonates
        
        # Check if the sum of protonated states - the of deprotonated state = 1
        sum = float("%.2f" % sum)
        
        if sum != 1.00:
            # Check if residue is a positive residue
            if fileName[0:3] in self.positive:
                correcao = 1.0 - sum
                stCorrect = open(os.path.join(self._pathPqr, fileName))
                for line in stCorrect:
                    tmp = line.split()
                    if tmp.__len__() == 1:
                        s = s + tmp[0] + '\n'
                    else:
                        atm = self.getAtom4chargeInResidue('Residue', fileName)
                        if tmp[1] == atm:
                            c = float(tmp[3])+correcao
                            c = c*(-1.0)
                            s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + str("%.3f" % c) + '\n'
                        else:
                            s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + tmp[3] + '\n'
            # Check if residue is a negative residue                
            elif fileName[0:3] in self.negative:
                correcao = 1.0 + sum
                stCorrect = open(os.path.join(self._pathPqr, fileName))
                for line in stCorrect:
                    tmp = line.split()
                    if tmp.__len__() == 1:
                        s = s + tmp[0] + '\n'
                    else:
                        atm = self.getAtom4chargeInResidue('Residue', fileName)
                        if tmp[1] == atm:
                            c = float(tmp[2])+correcao
                            s = s + tmp[0] + '\t' + tmp[1] + '\t' + str("%.3f" % c) + '\t' + tmp[3] + '\n'
                        else:
                            s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + tmp[3] + '\n'
            
            self.saveFile(fileName, s)
            
            
    # Check the correct charges in N-Terminal st file
    # Residue is NTres    
    def checaNT(self, fileName):
        st = open(os.path.join(self._pathPqr, fileName))
        atoms = {}
        s = ''
        
        for line in st:
            tmp = line.split()
            if tmp.__len__() > 1:
                atoms[tmp[1]] = tmp[2]
                      
        sum = 0.0
        for a in atoms:
            sum = sum + float(atoms[a])
        
        # Check if the sum of protonated states - the of deprotonated state = 1
        sum = float("%.2f" % sum)
        if sum != 1.00:
            correcao = 1.0 - sum
            stCorrect = open(os.path.join(self._pathPqr, fileName))
            for line in stCorrect:
                tmp = line.split()
                if tmp.__len__() == 1:
                    s = s + tmp[0] + '\n'
                else:
                    atm = self.getAtom4chargeInResidue('NT', fileName)
                    if tmp[1] == atm:
                        c = float(tmp[3])+correcao
                        c = c*(-1.0)
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + str("%.3f" % c) + '\n'
                    else:
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + tmp[3] + '\n'
            
            self.saveFile(fileName, s)
        
    # Check the correct charges in N-Terminal st file
    # Residue is resNT. In this case the residue has the behavior 
    # as residue and not like N-Terminal
    def checaSiteNT(self, fileName):
        st = open(os.path.join(self._pathPqr, fileName))
        atomsPROT = {}
        atomsDEPROT = {}
        s = ''
        
        for line in st:
            tmp = line.split()
            if tmp.__len__() > 1:
                atomsPROT[tmp[1]] = tmp[2] # Charges on the first column
                atomsDEPROT[tmp[1]] = tmp[3] # Charges on the second column                    
                      
        sum = 0.0
        protonates = 0.0
        deprotonates = 0.0
        
        for a in atomsPROT:
            protonates = protonates + float(atomsPROT[a])
        for a in atomsDEPROT:
            deprotonates = deprotonates + float(atomsDEPROT[a])
            
        sum = protonates - deprotonates
        
        # Check if the sum of protonated states - the of deprotonated state = 1
        sum = float("%.2f" % sum)
        if sum != 1.00:
            stCorrect = open(os.path.join(self._pathPqr, fileName))
            for line in stCorrect:
                tmp = line.split()
                if tmp.__len__() == 1:
                    s = s + tmp[0] + '\n'
                else:
                    atm = self.getAtom4chargeInResidue('ResidueNT', fileName)
                    if tmp[1] == atm:
                        # Check if residue is positive
                        if fileName[0:3] in self.positive:
                            correcao = 1.0 - sum
                            c = float(tmp[3])+correcao
                            c = c*(-1.0)
                            s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + str("%.3f" % c) + '\n'
                        # Check if residue is negative
                        elif fileName[0:3] in self.negative:
                            correcao = 1.0 + sum
                            c = float(tmp[2])+correcao
                            s = s + tmp[0] + '\t' + tmp[1] + '\t' + str("%.3f" % c) + '\t' + tmp[3] + '\n'
                            
                    else:
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + tmp[3] + '\n'

            self.saveFile(fileName, s)
           
    # Check the correct charges in C-Terminal st file  
    # Residue is CTres  
    def checaCT(self, fileName):
        st = open(os.path.join(self._pathPqr, fileName))
        atoms = {}
        s = ''
        
        for line in st:
            tmp = line.split()
            if tmp.__len__() > 1:
                atoms[tmp[1]] = tmp[3]
                
        sum = 0.0
        for a in atoms:
            sum = sum + float(atoms[a])
        
        # Check if the sum of protonated states - the of deprotonated state = 1
        sum = float("%.2f" % sum)
        
        if sum != -1.00:
            correcao = 1.0 + sum
            stCorrect = open(os.path.join(self._pathPqr, fileName))
            for line in stCorrect:
                tmp = line.split()
                if tmp.__len__() == 1:
                    s = s + tmp[0] + '\n'
                else:
                    atm = self.getAtom4chargeInResidue('CT', fileName)
                    if tmp[1] == atm:
                        c = float(tmp[2])+correcao
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + str("%.3f" % c) + '\t' + tmp[3] + '\n'
                    else:
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + tmp[3] + '\n'
            
            self.saveFile(fileName, s)
            
    # Check the correct charges in C-Terminal st file  
    # Residue is resCT            
    def checaSiteCT(self, fileName):
        st = open(os.path.join(self._pathPqr, fileName))
        atoms = {}
        s = ''
        
        for line in st:
            tmp = line.split()
            if tmp.__len__() > 1:
                atoms[tmp[1]] = tmp[3]
                
        sum = 0.0
        for a in atoms:
            sum = sum + float(atoms[a])
        
        # Check if the sum of protonated states - the of deprotonated state = 1
        sum = float("%.2f" % sum)
        
        if sum != -1.00:
            correcao = 1.0 + sum
            stCorrect = open(os.path.join(self._pathPqr, fileName))
            for line in stCorrect:
                tmp = line.split()
                if tmp.__len__() == 1:
                    s = s + tmp[0] + '\n'
                else:
                    atm = self.getAtom4chargeInResidue('ResidueCT', fileName)
                    print atm
                    if tmp[1] == atm:
                        c = float(tmp[3])+correcao
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + str("%.3f" % c) + '\t' + tmp[3] + '\n'
                    else:
                        s = s + tmp[0] + '\t' + tmp[1] + '\t' + tmp[2] + '\t' + tmp[3] + '\n'
                
            self.saveFile(fileName, s)
        
            
    # Build the OGM file    
    def createOGM(self, p1=41, p2=8.0, p3=2.0, p4=0.25):
        #box = int(self.getRaioProtein(self._pqrName[0:self._pqrName.__len__()-4]))
        
        #if box % 2 == 0:
            #box += 1
        
        #p1 = box
        
        s = ''
        fileName = self._pqrName[0:self._pqrName.__len__()-4] + ".ogm"
        Origin = "ON_ORIGIN "+ str(p1) + " " + str(p2) 
        Center = "ON_CENT_OF_INTR "+ str(p1) + " " + str(p3)
        Center2 = "ON_CENT_OF_INTR "+ str(p1) + " " + str(p4)
        s = s + Origin + '\n'
        s = s + Center + '\n'
        s = s + Center2 + '\n'
        
        self.saveFile(fileName, s)
        
        
    # Build the MGM file    
    def createMGM(self, p1=41, p2=8.0, p3=2.0, p4=0.25):
        #box = int(self.getRaioProtein(self._pqrName[0:self._pqrName.__len__()-4]))
        
        #if box % 2 == 0:
            #box += 1
        
        #p1 = box
        
        s = ''
        fileName = self._pqrName[0:self._pqrName.__len__()-4] + ".mgm"
        Origin = "ON_GEOM_CENT "+ str(p1) + " " + str(p2) 
        Center = "ON_CENT_OF_INTR "+ str(p1) + " " + str(p3)
        Center2 = "ON_CENT_OF_INTR "+ str(p1) + " " + str(p4)
        
        s = s + Origin + '\n'
        s = s + Center + '\n'
        s = s + Center2 + '\n'
        
        self.saveFile(fileName, s)
        
         # Change the PQR file to permit various N and C Terminal
    # When the sites file 
    def changePQR(self):
        # Create a copy from genuine pqr file
        shutil.copy2(os.path.join(self._pathPqr, self._pqrName), os.path.join(self._pathPqr, self._pqrName + '.bkp'))
        self.cleanPQR()
        self.insertChainRemark()
        lastResidue = self.renumberPQR()
        return lastResidue
        
    # Remove from pqr file all fields that not is 'ATOM'
    def cleanPQR(self):
        pqr = open(os.path.join(self._pathPqr, self._pqrName)) 
        s = ''
        
        for line in pqr:
            tmp = line.split()
            if tmp[0] == 'ATOM':
                s = s + line
                    
        self.saveFile(self._pqrName + '.clean', s)
        
    # Insert the remark 'END_CHAIN' to indicate the end of each chain
    def insertChainRemark(self):
        pqr = open(os.path.join(self._pathPqr, self._pqrName + '.clean')) 
        s = ''
        endLine = False
        first = True
        
        firstLine = pqr.readlines()[0]
        back = int(firstLine.split()[4])
        pqr = open(os.path.join(self._pathPqr, self._pqrName + '.clean'))
        
        for line in pqr:
            tmp = line.split()
            if int(tmp[4]) < back:
                s = s + 'END_CHAIN' + '\n'
                s = s + line
                back = -1
                endLine = True
            elif int(tmp[4]) == 1 and endLine == False and first == False:
                s = s + 'END_CHAIN' + '\n'
                s = s + line
                endLine = True
                first = True
            else:
                s = s + line
                
            if int(tmp[4]) != 1:
                endLine = False
                first = False
        
        self.saveFile(self._pqrName + '.endChain', s)
        
    # Change the number order of numRes (if necessary)
    # The numRes numbering is always sequential.
    def renumberPQR(self):
        # Format field
        ATOM_FORMAT="%-8s %-4s %-5s %-3s %-4s %8.3f%8.3f%8.3f %7.3f%6.3f \n"
        
        pqr = open(os.path.join(self._pathPqr, self._pqrName + '.endChain')) 
        s = ''
        numResREF = 0
        i = 0
        lastResidue = -1 # numSeq of the last residue. These residue will be use in method to create sites file
        
        for line in pqr:
            tmp = line.split()
            
            if tmp[0] == 'END_CHAIN':
                s = s + line
            
            else:
                
                # New numbering for residues in pqr file
                # Start in 1 and go to N
                if numResREF != int(tmp[4]):
                    numResREF = int(tmp[4])
                    i += 1
                    lastResidue = int(tmp[1])
            
                campo   = tmp[0]
                numSeq  = tmp[1]
                atom    = tmp[2]
                resName = tmp[3]
                numRes  = str(i)
                px      = float(tmp[5])
                py      = float(tmp[6])
                pz      = float(tmp[7])
                charge  = float(tmp[8])
                radii   = float(tmp[9])
    
                values = (campo, numSeq, atom, resName, numRes, px, py, pz, charge, radii) 
                s = s + (ATOM_FORMAT % values)
        
        self.saveFile(self._pqrName + '.renumbered', s) 
        return lastResidue          


    # Return the atoms or list of atom witch will receive the charges
    # for a correct build of st's file 
    def getAtom4chargeInResidue(self, type, residue):
        # If residue is a N-Terminal as N-Terminal
        if type == 'NT':
            return 'N'
        # If residue is a C-Terminal as C-Terminal
        elif type == 'CT':
            return 'O'
        # If residue is a N-Terminal as Residue
        elif type == 'ResidueNT':
            if residue.__contains__('LYS'):
                return 'NZ'
            elif residue.__contains__('GLU'):
                return None
            elif residue.__contains__('HIS'):
                return None
            elif residue.__contains__('ARG'):
                return 'NE'
            elif residue.__contains__('ASP'):
                return 'O'
            elif residue.__contains__('CYS'):
                return None
            elif residue.__contains__('TYR'):
                return None
        # If residue is a C-Terminal as Residue    
        elif type == 'ResidueCT':
            if residue.__contains__('LYS'):
                return 'NZ'
            elif residue.__contains__('GLU'):
                return 'N'
            elif residue.__contains__('HIS'):
                return None
            elif residue.__contains__('ARG'):
                return 'NE'
            elif residue.__contains__('ASP'):
                return 'O'
            elif residue.__contains__('CYS'):
                return None
            elif residue.__contains__('TYR'):
                return 'N'
        # Just single residue
        elif type == 'Residue':
            if residue.__contains__('LYS'):
                return 'NZ'
            elif residue.__contains__('GLU'):
                return 'CD'
            elif residue.__contains__('HIS'):
                return 'NE2'
            elif residue.__contains__('ARG'):
                return 'NE'
            elif residue.__contains__('ASP'):
                return 'CG'
            elif residue.__contains__('CYS'):
                return None
            elif residue.__contains__('TYR'):
                return 'N'
    
    # Save pqr.renumbered file as pqr file for to be used by mead.
    # Save a pqr file without END_CHAIN remark
    def saveNewPQRfile(self):
        s = ''
        pqrRenum = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))
        
        for line in pqrRenum.readlines():
            if not line.__contains__('END_CHAIN'):
                s = s + line
        
        self.saveFile(self._pqrName, s)
        
        
    def buildStsFIle(self):
        self.buildARGst()
        self.buildASPst()
        self.buildGLUst()
        self.buildHISst()
        self.buildLYSst()
        self.buildTYRst()
        
        
    # It buil N-Terminal st file        
    def createStNT(self, amino, lineNT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))  
        atoms  = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == 'ATOM': 
                if fields[4] == lineNT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]

                    if resName in self.resIonic:
                        if atom == 'N' or atom == 'H' or atom == 'H2' or atom == 'H3':
                            atoms[atom] = str(charge)
                    else:
                        atoms[atom] = str(charge)
                               
        self.mountStNT(amino, atoms)
        pqrFile.close()
    
    # It buil N-Terminal st file
    # N-term is charged positively
    def mountStNT(self, amino, atoms):
        s = ''
        pk = self.getPkas('NT')
        s = s + str(pk) + '\n'
        
        for a in atoms:
            s = s + amino[2:5] + '\t' + a + '\t' + atoms[a] + '\t' + '0.00' + '\n'
        
        self.saveFile(amino+'.st', s)
        
    # It buil N-Terminal st file
    # In this case the residue N-Terminal is an ionic residue    
    def createStNT2(self, amino, lineNT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))  
        atoms  = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == "ATOM":
                if fields[4] == lineNT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]
                    
                    if resName in self.resIonic:
                        if atom != 'N' and atom != 'H' and atom != 'H2' and atom != 'H3':
                            atoms[atom] = str(charge)
                               
        self.mountStNT2(amino, atoms)
        pqrFile.close()
        
    # It buil N-Terminal st file
    # N-term is charged positively       
    # In this case the residue N-Terminal is an ionic residue 
    def mountStNT2(self, amino, atoms):
        s = ''
        pk = self.getPkas(amino[0:3])
        s = s + str(pk) + '\n'
        
        for a in atoms:
            s = s + amino[0:3] + '\t' + a + '\t' + atoms[a] + '\t' + '0.00' + '\n'
        
        self.saveFile(amino+'.st', s)
    
        
    # It buil C-Terminal st file
    # C-term is charged negatively 
    def createStCT(self, amino, linhaCT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))
        atoms = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == "ATOM":
                if fields[4] == linhaCT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]
                
                    if resName in self.resIonic:
                        if atom == 'C' or atom == 'O' or atom == 'OH' or atom == 'OXT':
                            atoms[atom] = str(charge)
                    else:
                        atoms[atom] = str(charge)
        
        self.mountStCT(amino, atoms)
        pqrFile.close()
        
    # It buil C-Terminal st file
    def mountStCT(self, amino, atoms):
        s = ''
        pk = self.getPkas('CT')
        s = s + str(pk) + '\n'

        for a in atoms:
            s = s + amino[2:5] + '\t' + a + '\t' + '0.00' + '\t' + atoms[a] + '\n'
        
        self.saveFile(amino+'.st', s)
        
        
    # It buil C-Terminal st file
    # In this case the residue C-Terminal is an ionic residue   
    def createStCT2(self, amino, linhaCT):
        pqrFile = open(os.path.join(self._pathPqr, self._pqrName + '.renumbered'))
        atoms = {}
       
        for line in pqrFile.readlines():
            fields = line.split()
            if fields[0] == "ATOM":
                if fields[4] == linhaCT:
                    campo   = fields[0]
                    numSeq  = fields[1]
                    atom    = fields[2]
                    resName = fields[3]
                    numRes  = fields[4]
                    px      = fields[5]
                    py      = fields[6]
                    pz      = fields[7]
                    charge  = fields[8]
                    radii   = fields[9]
                
                    if resName in self.resIonic:
                        if atom != 'C' and atom != 'O' and atom != 'OH' and atom != 'OXT':
                            atoms[atom] = str(charge)
        
        self.mountStCT2(amino, atoms)
        pqrFile.close()
        
    # It buil C-Terminal st file
    # C-term is charged negatively       
    # In this case the residue C-Terminal is an ionic residue 
    def mountStCT2(self, amino, atoms):
        s = ''
        pk = self.getPkas(amino[0:3])
        s = s + str(pk) + '\n'

        for a in atoms:
            s = s + amino[0:3] + '\t' + a + '\t' + '0.00' + '\t' + atoms[a] + '\n'
        
        self.saveFile(amino+'.st', s)
        
    # Build st file of residue ARG
    # ARG is a POSITIVE residue        
    def buildARGst(self):
        s = ''
        s = s + '12.0 \n'
        s = s + 'ARG N   -0.400 -0.400 \n'
        s = s + 'ARG CA  -0.000 -0.000 \n'
        s = s + 'ARG HN   0.400  0.400 \n'
        s = s + 'ARG HA   0.000  0.000 \n'
        s = s + 'ARG C    0.550  0.550 \n'
        s = s + 'ARG O   -0.550 -0.550 \n'
        s = s + 'ARG CB   0.000  0.000 \n'
        s = s + 'ARG HB1  0.000  0.000 \n'
        s = s + 'ARG HB2  0.000  0.000 \n'
        s = s + 'ARG CG   0.000  0.000 \n'
        s = s + 'ARG HG1  0.000  0.000 \n'
        s = s + 'ARG HG2  0.000  0.000 \n'
        s = s + 'ARG CD   0.350  0.280 \n'
        s = s + 'ARG HD1  0.000  0.000 \n'
        s = s + 'ARG HD2  0.000  0.000 \n'
        s = s + 'ARG NE  -0.350 -0.560 \n'
        s = s + 'ARG CZ   0.350  0.280 \n'
        s = s + 'ARG NH1 -0.700 -0.750 \n'
        s = s + 'ARG NH2 -0.700 -0.750 \n'
        s = s + 'ARG HE   0.450  0.000 \n'
        s = s + 'ARG HH1  0.400  0.375 \n'
        s = s + 'ARG HH2  0.400  0.375 \n'
        s = s + 'ARG HH11 0.400  0.375 \n'
        s = s + 'ARG HH12 0.400  0.375 \n'
        s = s + 'ARG HH21 0.400  0.375 \n'
        s = s + 'ARG HH22 0.400  0.375 \n'
        self.saveFile('ARG.st', s)
        
    # Build st file of residue GLU
    # GLU is a NEAGTIVE residue
    def buildGLUst(self):
        s = ''
        s = s + '4.4 \n'
        s = s + 'GLU N   -0.400  -0.400 \n'
        s = s + 'GLU CA  -0.000  -0.000 \n'
        s = s + 'GLU HN   0.400   0.400 \n'
        s = s + 'GLU HA   0.000   0.000 \n'
        s = s + 'GLU C    0.550   0.550 \n'
        s = s + 'GLU O   -0.550  -0.550 \n'
        s = s + 'GLU CB   0.000   0.000 \n'
        s = s + 'GLU HB1  0.000   0.000 \n'
        s = s + 'GLU HB2  0.000   0.000 \n'
        s = s + 'GLU CG   0.000  -0.000 \n'
        s = s + 'GLU HG1  0.000   0.000 \n'
        s = s + 'GLU HG2  0.000   0.000 \n'
        s = s + 'GLU CD   0.550   0.100 \n'
        s = s + 'GLU OE1 -0.495  -0.550 \n' 
        s = s + 'GLU OE2 -0.490  -0.550 \n'
        s = s + 'GLU HE2  0.435   0.000 \n'
        self.saveFile('GLU.st', s)
        
    # Build st file of residue HIS
    # HIS is a POSITIVE residue
    def buildHISst(self):
        s = ''
        s = s + '6.3 \n'
        s = s + 'HIS  N   -0.400  -0.400 \n'
        s = s + 'HIS  CA  -0.000  -0.000 \n'
        s = s + 'HIS  HN   0.400   0.400 \n'
        s = s + 'HIS  HA   0.000   0.000 \n'
        s = s + 'HIS  C    0.550   0.550 \n'
        s = s + 'HIS  O   -0.550  -0.550 \n'
        s = s + 'HIS  CB   0.125   0.125 \n'
        s = s + 'HIS  CG   0.142   0.155 \n'
        s = s + 'HIS  ND1 -0.350  -0.560 \n'
        s = s + 'HIS  CE1  0.141   0.155 \n'
        s = s + 'HIS  NE2 -0.350  -0.400 \n'
        s = s + 'HIS  CD2  0.142  -0.125 \n'
        s = s + 'HIS  HB1  0.000   0.000 \n'
        s = s + 'HIS  HB2  0.000   0.000  \n'
        s = s + 'HIS  HD1  0.450   0.000 \n'
        s = s + 'HIS  HE1  0.125   0.125 \n'
        s = s + 'HIS  HE2  0.450   0.400 \n'
        s = s + 'HIS  HD2  0.125   0.125 \n'
        self.saveFile('HIS.st', s)
        
    # Build st file of residue LYS
    # LYS is a POSITIVE residue
    def buildLYSst(self):
        s = ''
        s = '10.4 \n'
        s = s + 'LYS  N   -0.400  -0.400 \n'         
        s = s + 'LYS  CA  -0.000  -0.000 \n'
        s = s + 'LYS  HN   0.400   0.400 \n'
        s = s + 'LYS  HA   0.000   0.000 \n'
        s = s + 'LYS  C    0.550   0.550 \n'
        s = s + 'LYS  O   -0.550  -0.550 \n'
        s = s + 'LYS  CB   0.000   0.000 \n'
        s = s + 'LYS  HB1  0.000   0.000 \n'
        s = s + 'LYS  HB2  0.000   0.000 \n'
        s = s + 'LYS  CG   0.000   0.000 \n'
        s = s + 'LYS  HG1  0.000   0.000 \n'
        s = s + 'LYS  HG2  0.000   0.000 \n'
        s = s + 'LYS  CD   0.000   0.000 \n'
        s = s + 'LYS  HD1  0.000   0.000 \n'
        s = s + 'LYS  HD2  0.000   0.000 \n'
        s = s + 'LYS  CE   0.330   0.000 \n'
        s = s + 'LYS  HE1  0.000   0.000 \n'
        s = s + 'LYS  HE2  0.000   0.000 \n'
        s = s + 'LYS  NZ  -0.320  -0.780 \n'
        s = s + 'LYS  HZ1  0.330   0.390 \n'
        s = s + 'LYS  HZ2  0.330   0.390 \n'
        s = s + 'LYS  HZ3  0.330   0.000 \n'
        self.saveFile('LYS.st', s) 
    
    # Build st file of residue TYR
    # TYR is a NEAGTIVE residue    
    def buildTYRst(self):
        s = ''
        s = s + '9.6 \n'        
        s = s + 'TYR N   -0.400  -0.400 \n'
        s = s + 'TYR CA  -0.000  -0.000 \n'
        s = s + 'TYR HN   0.400   0.400 \n'
        s = s + 'TYR HA   0.000   0.000 \n'
        s = s + 'TYR C    0.550   0.550 \n'
        s = s + 'TYR O   -0.550  -0.550 \n'
        s = s + 'TYR CB   0.125   0.125 \n'
        s = s + 'TYR HB1  0.000   0.000 \n'
        s = s + 'TYR HB2  0.000   0.000 \n'
        s = s + 'TYR CG  -0.125  -0.195 \n'
        s = s + 'TYR CD1 -0.125  -0.195 \n'
        s = s + 'TYR HD1  0.125   0.125 \n'
        s = s + 'TYR CE1 -0.125  -0.195 \n'
        s = s + 'TYR HE1  0.125   0.125 \n'
        s = s + 'TYR CZ   0.055  -0.070 \n'
        s = s + 'TYR OH  -0.490  -0.580 \n'
        s = s + 'TYR HH   0.435   0.000 \n'
        s = s + 'TYR CE2 -0.125  -0.195 \n'
        s = s + 'TYR HE2  0.125   0.125 \n'
        s = s + 'TYR CD2 -0.125  -0.195 \n'
        s = s + 'TYR HD2  0.125   0.125 \n'
        self.saveFile('TYR.st', s)
        
    # Build st file of residue ASP
    # ASP is a NEAGTIVE residue
    def buildASPst(self):
        s = ''
        s = s + '4.0 \n'        
        s = s + 'ASP N   -0.400  -0.400 \n'   
        s = s + 'ASP CA  -0.000  -0.000 \n'
        s = s + 'ASP HN   0.400   0.400 \n'
        s = s + 'ASP HA   0.000   0.000 \n'
        s = s + 'ASP C    0.550   0.550 \n'
        s = s + 'ASP O   -0.550  -0.550 \n'
        s = s + 'ASP CB   0.000   0.000 \n'
        s = s + 'ASP HB1  0.000   0.000 \n'
        s = s + 'ASP HB2  0.000   0.000 \n'
        s = s + 'ASP CG   0.550   0.100 \n'
        s = s + 'ASP OD1 -0.495  -0.550 \n'
        s = s + 'ASP OD2 -0.490  -0.550 \n'
        s = s + 'ASP HD2  0.435   0.000 \n'
        self.saveFile('ASP.st', s)
        
        # Save some file in somewhere
    def saveFile(self, fileName, conteudo):       
        f = FcfrpFile(self._pathOut, fileName)
        f.getAsFile("w")
        f.write(conteudo)
        f.close()