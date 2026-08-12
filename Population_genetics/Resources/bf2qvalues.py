#!/usr/bin/python

#Converts Bayes factors to q-values
#Use BAYENV results file as input and provide the number of env. variables as second arg.

import sys, string, re,  os, commands, time, math

import numpy as np


def convert2qvalues(filename, nvars):

    nvars = int(nvars)

    snp_dict = {}
    
    data = open(filename, "r")

    pi_1 = 0.1 #Prior prob.

    lines = data.readlines()
    for line in lines:
        cols = line.split("\t")
        snp_name = cols[0]
        snp_vals = []
        for i in range(1,nvars+1):
            bf = cols[i].strip()
            pp = (1/(1+((1-pi_1)/(float(bf)*pi_1))))
            snp_vals.append(pp)
        snp_dict[snp_name] = snp_vals

    snp_qvalues = {}
    
    
    for snp in snp_dict:
        snp_qvalues[snp] = []
        
    for i in range(0,nvars):
        for snp1 in snp_dict:
            pp1 = snp_dict[snp1][i]
            qval = 0
            number_of_significant = 0
            for snp2 in snp_dict:
                pp2 = snp_dict[snp2][i]
                if pp2 > pp1:
                    qval = qval + (1-pp2)
                    number_of_significant += 1
            if qval > 0.0:
                qval = qval/number_of_significant
            else:
                print snp1
            snp_qvalues[snp1].append(qval)
        print "finished converting variable " + str(i+1) + "..."
    
    lines = "snp"
    for i in range(1, nvars+1):
        lines += "\tvar" + str(i)
    lines += "\n"
 
    print "Total number of SNPs: " + str(len(snp_qvalues))
    for key in snp_qvalues:
        lines +=  key
        for q in snp_qvalues[key]:
            lines += "\t" + str(q)
        lines += "\n"

    print "Writing the q-values to file: " + filename + ".qvalues.txt"
    qval_file = open(filename + ".qvalues.txt", "w")
    qval_file.write(lines)
    qval_file.close()

#Algorithm for to convert to q-values
#For each SNP i:
#Do
#Set qval<-0
#Set number_of_significant<-0
#For each SNP j:
#Do
#If PP(j) is greater than PP(i)
#qval <- qval + (1-PP(j))
#number_of_significant++
#Done
#qval<-qval/number_of_significant
#Done

if __name__ == '__main__':
    # Terminate if too few arguments
    if len(sys.argv) < 3:
        print 'usage: %s <file name> <number of env. vars>' % sys.argv[0]
        sys.exit(-1)
    main(sys.argv[1], sys.argv[2])
