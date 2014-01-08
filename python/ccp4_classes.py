#!/usr/bin/python
# A few classes needed to run CCP4 programs

## Modules and files needed
import os,sys,string
bhome=os.environ["BLEND_HOME"]                                            #  Directory where
if bhome[-1] == "/": sys.path.append(os.path.join(bhome,"src_python"))    #  all python code
if bhome[-1] != "/": sys.path.append(os.path.join(bhome,"/src_python"))   #  is stored
from j_baseClasses import CExeCode

## Classes

#########
class RunJob(CExeCode):

 def __init__(self,batch=None,log=None):
  CExeCode.__init__(self)
  self._cmdline=""

  ## Take care of batch and log file, if any
  if batch != None: self._batch_name=batch
  #if log != None: self._log_name=log
  if log != None: self.log_name=log

 def setCmdLineAndFiles(self,cmdline):
  self._cmdline=cmdline

#########
class RunGENERIC(RunJob):
 """ To run any CCP4 program and produce the corresponding log file """

 ## Constructor
 def __init__(self,program_name,types_list,files_list,batch_file,log_file):
  RunJob.__init__(self,batch_file,log_file)
  tmpcmdline=[program_name]
  for i in range(len(types_list)):
   tmpcmdline.append(types_list[i])
   tmpcmdline.append(files_list[i])
  tmpcmdline.append("<")
  tmpcmdline.append(batch_file)
  cmdline=string.join(tmpcmdline) 
  self.setCmdLineAndFiles(cmdline)

#########
