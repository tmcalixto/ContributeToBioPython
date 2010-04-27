from ConfigParser import *

class FcfrpConfig:
    # Stores the unique Singleton instance-
    _iInstance = None
    
    #fileConf = '/usr/share/tomcat5.5/webapps/prometheus/WEB-INF/classes/fcfrp/ConfigFCFRP.conf'
    fileConf = 'ConfigFCFRP.conf'
            
    class Singleton:
        def __init__(self):
            self.config_control = None
            
    def __init__(self,section,pathfilename=fileConf):
        if FcfrpConfig._iInstance is None:
            FcfrpConfig._iInstance = FcfrpConfig.Singleton()
            self._config = ConfigParser()
            #section is a recognize and its can use to divide the file in any parts 
            self._section = section
            self._pathFileName=pathfilename
            self._loadConfigure()
        # Store instance reference as the only member in the handle
        self.__dict__['_EventHandler_instance'] =  FcfrpConfig._iInstance   
    
    def __getattr__(self, aAttr):
        return getattr(self._iInstance, aAttr)
    
    def __setattr__(self, aAttr, aValue):
        return setattr(self._iInstance, aAttr, aValue)
        
    def _loadConfigure(self):
        self._config.read(self._pathFileName)

    def getsection(self):
        #As section it's possible any kind, it's necessary to know its value. Therefore, this method inform its value
        return self._section
        
    def getParameter(self,paramName):
        return self._config.get(self._section, paramName)

    def setParameter(self,paramName,value):
        self._config.set(self._section, paramName, value)
