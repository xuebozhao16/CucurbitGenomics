#!/local/bin/python

from locus import *
from convert import *
from run_bayenv import *
from compute_average_bf import *
from bf2qvalues import *
from standardize import *

import sys, string, re,  os, commands, time, random

#from random import *

DEBUG = False
SKIP_COVAR = False #True if skipping building the null model
SKIP_TEST = False #True if tests are skipped
NULL_FILE = "" #Alternative Null file
DIFF_NULL = False #True if different null file

NUM_TESTS = 0 #Number of tests performed

NUM_POP = 0 #Number of populations in the file
NUM_ENV = 0 #Number of environmental variables

NULL_SIZE = 0
TEST_SIZE = 0
ITERATIONS = "5000" #Number of iterations for the tests
CUT_OFF = 0.0 #Default cut off value (no cut off)
EPSILON = 0.1 #Small constant

CLEAN = False

#Removing temporary files
def clean():
    
    print "Removing temporary files..."
    sys.stdout.flush()
    failure, output = commands.getstatusoutput("rm covar*")
    failure, output = commands.getstatusoutput("rm *.out")
    print "done."


############### Generate commands for BAYENV2 #################

def gen_commands(env_file, no_env_var):

    global NUM_TESTS
    global ITERATIONS
    global NUM_POP

    cmds = []
    num_tests = NUM_TESTS
    rand_seeds = [int(random.uniform(1,99999)) for i in range(1, num_tests+1)]

    additional_cmds = "" #Change to " -c" if computing rho

    for i, rand_seed in enumerate(rand_seeds):
        cmd = "bayenv2 -i %s -m mean_covar.txt -e " + env_file + \
            " -p " + str(NUM_POP) + \
            " -k " + str(ITERATIONS) + \
            " -n " + str(no_env_var) + \
            " -t -r " + str(rand_seed) + " -o %s" + additional_cmds

        f = open("test-cmd.txt", "a")
        f.write("Test " + str(i+1) + ":\n")
        f.write(cmd)
        f.write("\n")
        f.close()
        cmds.append(cmd)

    return cmds

    
####################### Main #############################
def main(in_file, var_file):

    global NULL_FILE
    global ITERATIONS
    global NUM_POP
    global NULL_SIZE
    global TEST_SIZE
    global CUT_OFF
    global NUM_ENV
    global NUM_TESTS
    global NULL_FILE

    max_diff = 0.5
    test = 0
    cut_off_percent = 99.5 #Percent to cut of from testing.

    #Standardizing the env variables
    env_file, no_env_var  = standardize_env(var_file) 

    #Temporary measures
    no_env_var = NUM_ENV #Number of environment variables to compute...

    #Generate the BAYENV command - IN USE!!!
    cmds = gen_commands(env_file, no_env_var)
    
    #Create the null modell
    if (SKIP_COVAR and SKIP_TEST):
        print "\nSkipping creating the covariance matrix and the test part of BAYENV"
        print "\nConverting to the BAYENV file format..."
        convert_fileformat(in_file, True)

    elif (SKIP_COVAR):
        print "Skipping null model..."
        out_file, new_out_file = convert_fileformat(in_file, True)
        test_all_snps_multip(new_out_file, cmds, TEST_SIZE)
    elif (SKIP_TEST):
        print "Skip the test part of BAYENV"
        try:
            if DIFF_NULL:
                print "Computing different null = " + NULL_FILE
                out_file, new_out_file = convert_fileformat(NULL_FILE)
                compute_null_model_bayenv2(num_pops, ITERATIONS, out_file)
                write_mean_covar_bayenv2()
            else:
                print "Computing same null = " + in_file
                out_file, new_out_file = convert_fileformat(in_file)
                compute_null_model_bayenv2(num_pops, ITERATIONS, out_file)
                write_mean_covar_bayenv2()
                
        except:
            print "Could not perform tests, no covariance matrix available!"

    else:
            #Split the output covar file into one for each 5000 iteration.
        if DIFF_NULL:
            print "Computing different null = " + NULL_FILE
            out_file, not_out_file = convert_fileformat(NULL_FILE)
            print out_file
            compute_null_model_bayenv2(NUM_POP, ITERATIONS, out_file)
            print "Co-variance matrices created created...1"
            write_mean_covar_bayenv2()
            #Run multiprocessing tests - THE ACTUAL TEST!!!
            out_file, new_out_file = convert_fileformat(in_file)
            #print new_out_file
            test_all_snps_multip(new_out_file, cmds, TEST_SIZE)
        else:
            out_file, new_out_file = convert_fileformat(in_file, True)
            print out_file
            compute_null_model_bayenv2(NUM_POP, ITERATIONS, out_file)
            print "Co-variance matrices created created...2"
            write_mean_covar_bayenv2() 
            #Run multiprocessing tests - THE ACTUAL TEST!!!
            test_all_snps_multip(new_out_file, cmds, TEST_SIZE)

    print ""

    print "Calculating average and median bf..."
    compute_average_bf(NUM_ENV, NUM_TESTS)
    print "Converting BF to q-values..."
    convert2qvalues("results/average_bf.txt", NUM_ENV)

  
if __name__ == '__main__':
    
    # Terminate if too few arguments
    if len(sys.argv) < 3:
        print 'usage: %s <infile> <varfile>' % sys.argv[0]
        sys.exit(-1)
    main(sys.argv[1], sys.argv[2])
