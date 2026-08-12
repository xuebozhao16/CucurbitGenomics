#!/usr/bin/python

import sys, string, re,  os, commands, time, math

#from scipy import stats

#import scipy as sp

import numpy as np

#import matplotlib as mpl

#from matplotlib import pyplot as plt

class SNP:
    
    def __init__(self, name, num_env, t):
        
        self.name = name
        self.num_env = [False] * num_env
        self.bf_list = [[0 for i in range(t)] for j in range(num_env)]
        self.rel_signal = []
        self.sum_signals = 0
        self.lg_info = []
        self.chr = 99
        self.lg = 99
        
    def get_name(self):
        return self.name
            
    def get_num_env(self):
        return self.num_env

    def set_num_env(self, n):
        self.num_env[n] = True  
      
    def add_to_list(self, bf, k, i):
        self.bf_list[k][i] = bf

    def set_signal(self, gamma):
        self.rel_signal.append(gamma)
        self.sum_signals += gamma #Add to the total of signals

    #Return the bf signal in variable k
    def get_signal(self, k):
        return self.rel_signal[k]

    #Return the bf signal list
    def get_signals(self):
        return self.rel_signal

    def get_sum_signals(self):
        return self.sum_signals
    
    def print_env(self):
        print self.num_env   

    def get_median_bf(self, k):
        #print self.bf_list[k]
        bfs = np.array(self.bf_list[k])
        median = np.median(bfs)
        return median

    def get_avg_bf(self, k):
        #print self.bf_list[k]
        bfs = np.array(self.bf_list[k])
        avg = np.average(bfs)
        return avg

    def add_bf(self, bf):
        self.sum_bf += bf
        
    def get_sum_bf(self):
        return self.sum_bf
    
    def get_num_runs(self):
        return self.num_runs

    def get_bf_list(self):
        return self.bf_list

    def get_bf_list(self):
        return self.bf_list

    def set_lg_info(self, info):
        self.lg_info.append(info)


    def get_lg_info(self):
        return self.lg_info

    def set_chr(self, ch):
        self.chr = ch
    
    def get_chr(self):
        return self.chr

    def set_linkage_group(self, lg):
        self.lg = lg

    def get_linkage_group(self):
        return self.lg
  

def compute_average_bf(num_var, num_tests):

    N = int(num_var)
    t = int(num_tests)
    
    snp_dict = {}

    for i in range (0, t):
        filename = "results/bf_results_t" + str(i) + ".bf"
        data = open( filename, "r")
        print filename
        lines = data.readlines()
        for line in lines:
            cols = line.split("\t")
            snp_name = cols[0][0:-2]
            if i > 9:
                snp_name = snp_name[0:-1]                
            if snp_name in snp_dict:
                snp = snp_dict[snp_name]
                for k in range(0, N):
                    snp.add_to_list(float(cols[k+1]), k, i)
            else:
                snp = SNP(snp_name, N, t)
                snp_dict[snp_name] = snp
                for k in range(0, N):
                    snp.add_to_list(float(cols[k+1]), k, i)
                
        data.close()

    print "################LENGTH:" + str(len(snp_dict))

    FILE1 = open("results/median_bf.txt", "w")
    FILE2 = open("results/average_bf.txt", "w")
    
    #bf_median = "marker\tsal1\tsal2\ttemp1\ttemp2\tox1\tox2\n"
    #bf_avg = "marker\tsal1\tsal2\ttemp1\ttemp2\tox1\tox2\n"    
    bf_median = ""
    bf_avg = ""
    for key in snp_dict:
        snp = snp_dict[key]
        bf_avg += snp.get_name()
        bf_median += snp.get_name()
        for k in range(0, N):
            bf_a = snp.get_avg_bf(k)
            bf_m = snp.get_median_bf(k)
            bf_avg += "\t" + str(bf_a)
            bf_median += "\t" + str(bf_m)
        bf_avg += "\n"
        bf_median += "\n"
    FILE1.write(bf_median)
    FILE2.write(bf_avg)
    FILE1.close()
    FILE2.close()

if __name__ == '__main__':
    # Terminate if too few arguments
    if len(sys.argv) < 3:
        print 'usage: %s <number of vars> <num tests>' % sys.argv[0]
        sys.exit(-1)
    main(sys.argv[1], sys.argv[2])
