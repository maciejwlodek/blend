#!/usr/bin/python
#***********************************************************************************************************
#***********************************************************************************************************
#********* merge_clusters.py
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
#
# 
#

## Modules needed
import os,sys,shutil,string,math
#bhome=os.environ["BLEND_HOME"]                                            #  Directory where
#if bhome[-1] == "/": sys.path.append(os.path.join(bhome,"python"))    #  all python code
#if bhome[-1] != "/": sys.path.append(os.path.join(bhome,"/python"))   #  is stored
from ccp4_functions import generate_generic_batch_file
from ccp4_classes import RunGENERIC

## CCP4 programs used
program_mtzdump="mtzdump"
program_rebatch="rebatch"
program_pointless="pointless"

#########################################################
#######################  FUNCTIONS ######################
#########################################################
def next10(x):
 # This function returns the x next closest multiple of 10

 y=float(x)/10
 z=math.ceil(y)*10
 if int(z) == int(x):
  w=int(x+10)
 else:
  w=int(z)

 return w


#########################################################
#########################  MAIN #########################
#########################################################

## Input
li=sys.argv
filename=li[1]
file=open(filename,'r')
filecontents=file.readlines()
file.close()

# Directory where all merged files will be stored (it should already exists)
cl_dir=os.path.dirname(filename)

cluster_ini=[]
cl_number=[]
for i in range(len(filecontents)):
 linea=filecontents[i]
 if linea[0:3] == "***":
  wrds=linea[10:]
  cl_number.append(string.atoi(string.split(wrds)[1]))
  cluster_ini.append(i+1)
cluster_fin=[]
for i in range(len(cluster_ini)-1):
 cluster_fin.append(cluster_ini[i+1]-2)
cluster_fin.append(len(filecontents)-1)

## Check for CCP4 environment
is_ccp4=os.getenv("CCP4",0)
if is_ccp4 == 0: raise SystemExit("CCP4 has not been set up on this terminal. Please provide and start again.")

## Big loop over all clusters
for iclust in range(len(cluster_ini)):

 ## Suffix for merged file (cluster number)
 clnumber="%03d" % cl_number[iclust]
 stmp="Processing cluster %d:" % cl_number[iclust]
 print stmp

 ## Build selection of mtz files to be merged together
 selection=range(cluster_ini[iclust],cluster_fin[iclust])
 selection.append(cluster_fin[iclust])

 ## Build list with number of batches to be considered for each mtz file
 n_batches=[]
 for j in selection: n_batches.append(string.atoi(string.split(filecontents[j])[2]))

 ## Loop to obtain one merged file
 files_to_delete=[]
 files_for_pointless=[]
 ini_batch=1
 keywords_tmp=["ASCEND","H K L M/ISYM BATCH"]
 merge_list=["merged_from"]
 for j in range(len(selection)):
 
  ## Suffix for various files
  #anumber="%03d" % selection[j]
  anumber="%03d" % j
  merge_list.append(string.join(["_",anumber],""))
  
  ## Input file
  mtz_in=string.split(filecontents[selection[j]])[0]

  ## Extracting number of batches for this file using mtzdump
  batch_file=os.path.join(os.curdir,string.join(["mtzdump_",anumber,".com"],""))
  files_to_delete.append(batch_file)
  log_file=os.path.join(os.curdir,string.join(["mtzdump_",anumber,".log"],""))
  files_to_delete.append(log_file)
  keywords_list=["END"]
  generate_generic_batch_file(batch_file,keywords_list)
  types_list=["HKLIN"]
  files_list=[mtz_in]
  jobGENERIC=RunGENERIC(program_mtzdump,types_list,files_list,batch_file,log_file)
  jobGENERIC.execute()
  file=open(log_file,'r')
  fconts=file.readlines()
  file.close()
  nbatch=0
  blist=[]
  for iline in range(len(fconts)):
   line=fconts[iline]
   ltmp=string.split(line)
   if len(ltmp) != 0:
    if ltmp[0] == 'Batch':
     nbatch+=1
     blist.append(string.atoi(string.split(fconts[iline+1])[0]))
  print "MTZ file %s contains %3d batches. Number of batches used: %3d." % (mtz_in,nbatch,n_batches[j])

  ## Renumbering batches using rebatch
  mtz_out=os.path.join(os.curdir,string.join(["from_rebatch_",anumber,".mtz"],""))
  files_to_delete.append(mtz_out)
  keywords_tmp.append(string.join(['"',mtz_out,'"'],""))
  batch_file=os.path.join(os.curdir,string.join(["rebatch_",anumber,".com"],""))
  files_to_delete.append(batch_file)
  log_file=os.path.join(os.curdir,string.join(["rebatch_",anumber,".log"],""))
  files_to_delete.append(log_file)
  bnumber="%d" % ini_batch
  keywords_list=[string.join(["TITLE Rebatch file",mtz_in]),string.join(["BATCH ALL START",bnumber]),"END"]
  if len(n_batches) != 0:
   keywords_list=["TITLE Rebatch file for multicrystal merging",
                  string.join(["BATCH",repr(blist[0]+n_batches[j]),"TO",repr(blist[0]+nbatch-1),"REJECT"]),
                  string.join(["BATCH ALL START",bnumber]),"END"]
  generate_generic_batch_file(batch_file,keywords_list)
  files_list=[mtz_in,mtz_out]
  types_list=["HKLIN","HKLOUT"]
  jobGENERIC=RunGENERIC(program_rebatch,types_list,files_list,batch_file,log_file)
  jobGENERIC.execute()
  ini_batch+=nbatch-1
  ini_batch=next10(ini_batch)+1
 
  files_for_pointless.append(mtz_out)

 ## Pointless for alternate indexing, determining correct point group and final merging
 if len(selection) > 1:
  stmp="Merging data from cluster %d into single mtz file .........\n" % cl_number[iclust]
  print stmp
  batch_file=os.path.join(os.curdir,"pointless.com")
  files_to_delete.append(batch_file)
  shutil.copyfile(files_for_pointless[0],os.path.join(os.curdir,"reference.mtz"))
  files_to_delete.append(os.path.join(os.curdir,"reference.mtz"))
  keywords_list=[string.join(["HKLREF ",os.path.join(os.curdir,"reference.mtz")],"")]
  for a in files_for_pointless:
   keywords_list.append(string.join(["HKLIN ",a],""))
  keywords_list.append(string.join(["HKLOUT ",os.path.join(os.curdir,"merged_tmp.mtz")],""))
  keywords_list.append("END")
  generate_generic_batch_file(batch_file,keywords_list)
  log_file=os.path.join(os.curdir,"pointless.log")
  files_to_delete.append(log_file)
  mtz_out=os.path.join(os.curdir,"merged_tmp.mtz")
  files_to_delete.append(mtz_out)
  files_list=[""]
  types_list=[""]
  jobGENERIC=RunGENERIC(program_pointless,types_list,files_list,batch_file,log_file)
  jobGENERIC.execute()

  ## Final space group
  batch_file=os.path.join(os.curdir,"pointless2.com")
  files_to_delete.append(batch_file)
  keywords_list=[string.join(["HKLIN ",mtz_out],"")]
  mtz_out=os.path.join(cl_dir,string.join(["merged",clnumber,".mtz"],""))
  keywords_list.append(string.join(["HKLOUT ",mtz_out],""))
  keywords_list.append("END")
  generate_generic_batch_file(batch_file,keywords_list)
  log_file=os.path.join(os.curdir,"pointless2.log")
  files_to_delete.append(log_file)
  files_list=[""]
  types_list=[""]
  jobGENERIC=RunGENERIC(program_pointless,types_list,files_list,batch_file,log_file)
  jobGENERIC.execute()
 else:
  stmp="Merging data from cluster %d into single mtz file .........\n" % (iclust+1)
  print stmp

  ## Final space group
  batch_file=os.path.join(os.curdir,"pointless2.com")
  files_to_delete.append(batch_file)
  keywords_list=[string.join(["HKLIN ",mtz_out],"")]
  mtz_out=os.path.join(cl_dir,string.join(["merged",clnumber,".mtz"],""))
  keywords_list.append(string.join(["HKLOUT ",mtz_out],""))
  keywords_list.append("END")
  generate_generic_batch_file(batch_file,keywords_list)
  log_file=os.path.join(os.curdir,"pointless2.log")
  files_to_delete.append(log_file)
  files_list=[""]
  types_list=[""]
  jobGENERIC=RunGENERIC(program_pointless,types_list,files_list,batch_file,log_file)
  jobGENERIC.execute()
 
 ## Clean unnecessary files
 for a in files_to_delete: os.remove(a)
