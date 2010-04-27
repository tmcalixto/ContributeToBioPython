import os
import shutil 

class FcfrpFile:
    
    def __init__(self,dir, name):
        self._dir = dir
        self._name = name
        self._path = os.path.join(dir, name)
        self._asFile = None

    def _inicialize(self,dir,name,path=""):
        self._dir = dir
        self._name = name
        self._path = os.path.join(dir, name)
        self._asFile = None
            
    def existsFile(self):
        return os.path.exists(self._path)
    
    def getDir(self):
        return self._dir
    
    def getFileName(self):
        return self._name
    
    def getPath(self):
        return self._path
    
    def getAsFile(self,mode="r"):
        #The reference will be returned as a file and its mode is choose by user
        if self._asFile != None:
            raise "The file has already opened. Therefore, close it before open again"
        self._asFile = open(self._path,mode)
        return self._asFile

    def close(self):
        #if opened as File, it's necessary close it
        if self._asFile != None:
            self._asFile.close()
            self._asFile = None
            
    def read(self):
        #Read a file
        f = open(self._path,"r")
        t = f.read()
        f.close()
        return t
    
    def write(self,value):
        if self._asFile != None:
            self._asFile.write(value)
        else:
            raise Exception("File isn't open as write mode")

    def readLines(self):
        # reads all lines and returns each as part of a list
        f = open(self._path,"r")
        l = f.readlines()
        f.close()
        return l
    
    def readline(self):
        if self._asFile != None:
            return self._asFile.readline()
        else:
            raise Exception("File isn't open as read mode")     
    
    def remove(self):
        os.remove(self._path)
        self._inicialize("", "", "")
    
    def rename(self,new):
        pathNewFileName = os.path.join(self._dir,new)
        #os.rename(self._path, pathNewFileName)
        shutil.copy(self._path, pathNewFileName)
        os.remove(self._path)
        self._name = new
        self._inicialize(self._dir, self._name, pathNewFileName)
        
    def find(self,tokenSearch,start=0):
        #Search tokenSearch value in the file and returns if found or not.
        #if not found the return is -1. Otherwise, > -1
        text = self.read()
        return text.find(tokenSearch,start)
    