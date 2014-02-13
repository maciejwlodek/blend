#!/usr/bin/python
#***********************************************************************************************************
#***********************************************************************************************************
#********* xds_to_mtz_list.py
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
# Simple script to turn xds files of an input list for BLEND into mtz, and update the list.
# This runs in the following way:
#
#          xds_to_mtz_list.py list.dat
#


## Modules needed
import os,sys,string, subprocess
#bhome=os.environ["BLEND_HOME"]                                            #  Directory where
#if bhome[-1] == "/": sys.path.append(os.path.join(bhome,"python"))    #  all python code
#if bhome[-1] != "/": sys.path.append(os.path.join(bhome,"/python"))   #  is stored
#from ccp4_functions import generate_generic_batch_file
#from ccp4_classes import RunGENERIC

## CCP4 programs used
#program_pointless="pointless"

########################
######### MAIN #########
########################

# Extract input from command line
li=sys.argv

# Stop if less or more than 2 words in command line
if len(li) > 3: raise SystemExit("Wrong command-line input")

# Read all file names from input file
file=open(li[1],'r')
filecontents=file.readlines()
file.close()
mtz_names=[]
for line in filecontents:
 if len(string.split(line,"\n")[0]) != 0: mtz_names.append(string.split(line,"\n")[0])

# Second word in command line is LAUEGROUP line
if len(li) == 3:
 lauegroup = li[2].strip()
 if lauegroup == "": lauegroup = None
else:
 lauegroup = None
if lauegroup is not None:
 ltmp = lauegroup.split()[1].upper()
 if ltmp == "AUTO": lauegroup = "AUTO"

# List with xds look up table information
xds_lu=[]

# Create file with list of mtz files
lines=[]
for iname in range(len(mtz_names)):
 name=mtz_names[iname]
 sfx="%03d.mtz" % (iname+1)
 newname=string.join(["dataset_",sfx],"")
 if os.path.splitext(name)[1] == ".mtz": lines.append(string.join([os.path.join(name),"\n"],""))
 if os.path.splitext(name)[1] == ".HKL" or os.path.splitext(name)[1] == ".hkl":
  # Name of subdirectory where all new mtz files are collected
  xds_files_dir=os.path.join(os.curdir,"xds_files")
  # First turn xds file into an mtz file (all new mtz files collected in directory xds_files)
  lu=string.join([newname,"     ",name,"\n"],"")
  xds_lu.append(lu)
  xdsin=string.join([os.path.join(name),"\n"],"")
  if not os.path.exists(xds_files_dir): os.mkdir(xds_files_dir)
  mtzout=os.path.join(xds_files_dir,newname)

  # Run POINTLESS to turn XDS files into MTZ files (this bit replace the following commented lines)
  if lauegroup == "AUTO":
   partline = "TOLERANCE  100\n" + "END\n"
  elif lauegroup is None:
   partline = None
  else:
   partline = "TOLERANCE  100\n" + lauegroup + "\n" + "END\n"
  if partline is None:
   cmdline = 'pointless -c xdsin ' + xdsin.rstrip("\n") + ' hklout ' + mtzout.rstrip("\n")
  else:
   cmdline = 'pointless xdsin ' + xdsin.rstrip("\n") + ' hklout ' + mtzout.rstrip("\n") + ' <<eof-pointless\n' + partline + 'eof-pointless'
  pout = string.join(["pointless_", "%03d" % (iname+1), ".log"],"")
  pout = os.path.join(xds_files_dir, pout)
  fileP = open(pout, 'w')
  try:
   ttt = subprocess.Popen(cmdline, shell = True, stdout = fileP)
   return_code = ttt.wait()   #!!! Important !!! Without wait() the shell launching this program can carry on without these results.
  except subprocess.CalledProcessError:
   raise SystemExit("Error has occurred while converting XDS file " + xdsin.rstrip("\n") + ".")
  fileP.close()

  # All temporary files will be deleted at the end of this script
  #files_to_delete=[]

  # Create keyword file for pointless
  #batch_file=os.path.join(os.curdir,"pointless.com")
  #files_to_delete.append(batch_file)
  #keywords_list=["NAME PROJECT From XDS to MTZ  CRYSTAL Xtal  DATASET Data"]
  #keywords_list.append(string.join(["XDSIN",xdsin]))
  #keywords_list.append(string.join(["HKLOUT",mtzout]))
  #keywords_list.append("COPY")
  #generate_generic_batch_file(batch_file,keywords_list)
  #log_file=os.path.join(os.curdir,"pointless.log")
  #files_to_delete.append(log_file)

  # Run pointless
  #files_list=[""]
  #types_list=[""]
  #jobGENERIC=RunGENERIC(program_pointless,types_list,files_list,batch_file,log_file)
  #jobGENERIC.execute()

  # Get rid of unwanted files
  #for a in files_to_delete:
  # os.remove(a)

  # Then list the new mtz file with the other
  lines.append(string.join([mtzout,"\n"],""))

# Open file to write all mtz - xds corespondences
if len(xds_lu) != 0:
 file=open("xds_lookup_table.txt",'w')
 file.writelines(xds_lu)
 file.close()

# Open file to write all mtz files to be used
file=open("mtz_names.dat",'w')
file.writelines(lines)
file.close()
