#!/usr/bin/python
# A few functions needed to run CCP4 programs

## Modules needed
import string

## Functions

#########
def generate_generic_batch_file(filename,keywords_list):
 """ Keywords file to serve as batch file for generic CCP4 program """

 file=open(filename,'w')
 filecontents=[]
 for linea in keywords_list:
  filecontents.append(string.join([linea,"\n"],""))
 file.writelines(filecontents)
 file.close()

#########
