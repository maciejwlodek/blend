#!/usr/bin/python
# Python script
#

# The date of last modification: 24/03/2005

import os
import sys
import glob
import shutil
import re
import math
import string
import time
import select, fcntl
import signal

#################################################   

class CExeCode :
    """ A generic abstract base class that is to be inheritted by other classes that warp 
        a executable code. Basically this class defines two main methods that characterize the
        procedures a job by an executable code, i.e.  forming a command line object and then
        run the code according to the command line object"""

    def __init__( self, path_obj = None ):

        # set pathes
        if path_obj :
            self._path_dict = path_obj
 
        self.exitCode = 0

        # test 
        self._setCmdLineAndFile()

    def _setCmdLineAndFile(self ):         # will be overrided by individual derived classes
        self._cmdline = 'ls '
        self.batch_name = None
        self.log_name = None
            
    def execute(self, mode = 0, err_str ="" ):

        if self.log_name :
            self.logMode(mode, err_str)
        else :
            self.interactiveMode()

    def logMode(self, mode = 0, err_str=""):
        """ execute self._cmdline and output to a log file (self.log_name)  
        in a nonblocking way"""

        logfile = open(self.log_name, 'w')
        logfile.write("\n============ PROCESS INFORMATION =============\n")

        # spawn a sub-process for a job and connect to its input/output(and error)
        # streams using pipes

        try : 
            import subprocess
            subProcess = subprocess.Popen(self._cmdline, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
            subProcess.stdin.close()
            outfile_sub  = subProcess.stdout
            outfd_sub    = outfile_sub.fileno()
        except ImportError:
            print "module 'subprocess' has not be found, you use a version of PYTHON below 2.4"
            import popen2
            subProcess = popen2.Popen4(self._cmdline)
            subProcess.tochild.close()
            outfile_sub  = subProcess.fromchild
            outfd_sub    = outfile_sub.fileno()
    
        self._pid  = subProcess.pid
        pid_c      = self._pid  + 1
        p_name = self._cmdline.strip().split()[0].split("/")[-1]
        t_log_name = self.log_name.strip().split("/")[-1]
        
       
        log_size = 0 
        log_max  = 100000000
        outfile_eof = 0
        if not mode:
            self._nonBlockingFile(outfd_sub)
            while not outfile_eof :
                ready_to_read, ready_to_write, in_error = \
                     select.select([outfd_sub],[],[],60.0)
                if outfd_sub in ready_to_read:
                    out_quantum = outfile_sub.read()
                    out_quantum_bt = len(out_quantum)
                    log_size = log_size + out_quantum_bt
                    if log_size > log_max:
                        logfile.write("\nTHE SIZE OF LOG FILE, %s, EXCEEDED THE LIMIT AND PROCESS STOPPED\n"%t_log_name)
                        logfile.write("The process running %s is killed\n"%p_name)
                        logfile.close()
                        os.kill(self._pid, signal.SIGKILL)
                        os.kill(pid_c, signal.SIGKILL)
                        print "The process running %s is killed " %p_name
                        time.sleep(1.0)
                        # be safe stop all process
                        sys.exit(1)
                    else :
                        if out_quantum == '':
                            outfile_eof = 1
                        logfile.write(out_quantum)
                        logfile.flush()
        else :
            self._nonBlockingFile(outfd_sub)
            while not outfile_eof :
                ready_to_read, ready_to_write, in_error = \
                     select.select([outfd_sub],[],[],10)
                if outfd_sub in ready_to_read:
                    out_quantum = outfile_sub.read()
                    if out_quantum == '':
                        outfile_eof = 1
                    logfile.write(out_quantum)
                    logfile.flush()
                    # allProcInfo.append(out_quantum)

        self.exitCode = subProcess.wait()
      
        # if mode:
        #    for item in allProcInfo:
        #        logfile.write(item)

        logfile.write("=========END OF PROCESS INFORMATION\n ==========\n")
        logfile.close()

        outfile_sub.close()

        if self.exitCode:   
            print "#-----------------------------------------------------------#"
            print "The process stoped because of a runtime error in code '%s'!"%p_name
            if err_str:
                print "%s"%err_str
            else:
                print "See the associated log file '%s'."%t_log_name
            print "#-----------------------------------------------------------#"
            if not err_str:
                # if normal runtime errors, stop the program(using log. no err_str)  
                # otherwise balbes will go further and let the next procedure to handle the error. 
                sys.exit(1)
                
        return True

    def _nonBlockingFile(self,fd):
        f_flag = fcntl.fcntl(fd, fcntl.F_GETFL,0)
        fcntl.fcntl(fd, fcntl.F_SETFL, f_flag | os.O_NONBLOCK)

    def interactiveMode(self):
        """ execute self._cmdline  just like what is done in a shell """
        
        os.system(self._cmdline)


# end
