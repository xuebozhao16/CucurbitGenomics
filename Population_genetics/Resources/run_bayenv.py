import sys, os
import commands
import time
import multiprocessing
import random

import numpy as np


# Process class to run the Bayenv test phase in paralell
class RunInProcess(multiprocessing.Process):
    
    def __init__(self, cmd, thread_id, testdata, testsize):
        multiprocessing.Process.__init__(self)
        self.cmd = cmd
        self.thread_id = thread_id
        self.testdata = testdata
        self.testsize = testsize

    def run(self):
        errors = open("errors.out", 'w')
        outfile = "results/bf_results_t" + str(self.thread_id)
        k = 0
        dataset = open(self.testdata, 'r')
        thread_id = str(self.thread_id)
        cmd = self.cmd
        data_size = self.testsize
        while dataset:
            k += 1
            line1 = dataset.readline()
            if line1 == '':
                break
            l = line1.split("\t")
            marker_file = l.pop(0) + "-" + str(self.thread_id)
            FILE = open(marker_file, 'w')
            line1 = "\t".join(l)
            line2 = dataset.readline()
            FILE.write(line1 + line2)
            FILE.close()
            print "BAYENV: process " + thread_id + " is processing " + marker_file + " (" + str(k) + ")...",
            start_test = time.time()
            sys.stdout.flush()
            failure, output = commands.getstatusoutput(cmd % (marker_file, outfile))
            elapsed = (time.time() - start_test)
            remaining = ((data_size-k)*elapsed)/60
            print "done. %f sec to complete. Estimated time remaining: %f minutes" % (elapsed, remaining)
        
            if failure:
                print output
                errors.write(output)
                error =  "Could not test locus: " + marker_file +"\n"\
                      + line1 + line2
                errors.write(error)
                os.remove(marker_file)
                os.remove(marker_file + ".freqs")
                continue
            os.remove(marker_file)
            os.remove(marker_file + ".freqs")
        errors.close()
        dataset.close()


############# Running BAYENV with multiprocessing ###############
def test_all_snps_multip(testdata, cmds, testsize):
    
    procs = []

    start = time.time()
    for i in range(len(cmds)):
        proc = RunInProcess(cmds[i], i, testdata, testsize)
        procs.append(proc)
        proc.start()
    
    for p in procs:
        p.join()
    print "Elapsed time: %s" % (time.time()-start)


#Calculating covariance matrices every 500 iterations - Bayenv 2.0
def compute_null_model_bayenv2(num_pops, iterations, snpfile):

    print "#############################################"
    
    time_usage = open("t_usage.out", 'w')

    rand_seed = str(int(random.uniform(1,99999)))
    print "Random seed = " + rand_seed

    cmd_str = open("covar-cmd.txt", "wb")

    cmd = "bayenv2 -i " + snpfile + " -p " + str(num_pops) + " -k " + str(iterations) \
        + " -r " + str(rand_seed) + " > covars.txt"

    print cmd
    cmd_str.write(cmd)
    cmd_str.close()
   
    print "BAYENV calculating the covariance matrices...",
    start_test = time.time()
    sys.stdout.flush()
    failure, output = commands.getstatusoutput(cmd)
    elapsed = (time.time() - start_test)
    usage = "done. %f sec to complete" % elapsed
    print usage
    time_usage.write(usage)
    time_usage.close()


#Separates the covariance matrices from bayenv2 output and 
#writes one  mean co-variance matrix
def write_mean_covar_bayenv2():

    covars = open("covars.txt", 'r')    

    #Removing the first 15 lines
    for i in range(0, 15):
        covars.readline()
    
    covar_lists = []
    cov = []

    matrix_counter = 0
    for line in covars:
        if ("VAR-COVAR" in line):
            matrix_counter += 1
            covar_lists.append(cov)
            cov = []
        elif line == "\n":
            continue
        else:
            line = line.strip("\t\n")
            cov.append(line.split("\t"))

    num_cov_matrix = np.array(covar_lists, np.float64)

    matrix_mean = np.average(num_cov_matrix, axis=0)

    covar_string = ""

    for list in matrix_mean:
        for elem in list:
            covar_string += str(elem) + "\t"
        covar_string += "\n"
    
    f = open("covar-cmd.txt", "a")
    f.write(covar_string)
    f.close()

    covar_file = open("mean_covar.txt", 'w')
    covar_file.write(covar_string)
    covar_file.close() 
    covars.close() 
