import sys, string, re,  os, commands, time

from scipy import stats

import scipy as sp

import numpy as np



################## Class Env_var ######################
class Env_var:

    
    def __init__(self, v):
        self.std_env_list = []
        self.env_list = []
        v = v.replace('\r\n', '')
        v = v.replace(',', '.') 
        data = v.split('\t')
        self.name = data[0]
        
        variables = data[1::]
        for var in variables:
            self.env_list.append(float(var))
            
        self.env_data = np.array(self.env_list)
        
        mean = np.mean(self.env_data)
        std = np.std(self.env_data)
        
        self.std_env_data = np.ones((len(self.env_data), 1))
        for i in range(0, len(self.env_data)):
            std_var = (self.env_data[i] - mean)/std
            self.std_env_data[i] = std_var
        

    def print_env_data(self):
        print self.name
        print self.env_data

    def get_std_env(self):
        self.std_env_list = self.std_env_data.tolist()
        return self.std_env_list
    

        
        #self.pop = [[0 for x in xrange(2)] for x in xrange(pops)]



####################### Main #############################
def standardize_env(in_file):

    #num_env = "4"
    out_file = "std_" + in_file


    env_list = [] #List of locus names
    
    env_data = open(in_file, 'r')
    lines = env_data.readlines()

    num_vars =  len(lines)

    for line in lines:
        env = Env_var(line)
        env_list.append(env)

    env_lines = ""
    for i in range(0, len(env_list)):
        var = env_list[i].get_std_env()
        for v in var:
            env_lines += str(v) + '\t'
        env_lines += '\n'
    
    env_lines = env_lines.replace('[', '')
    env_lines = env_lines.replace(']', '')
    env_out = open(out_file, 'w')
    env_out.write(env_lines)
    env_out.close()

    return out_file, num_vars

  
if __name__ == '__main__':
    
    # Terminate if too few arguments
    if len(sys.argv) < 2:
        print 'usage: %s <infile>' % sys.argv[0]
        sys.exit(-1)
    main(sys.argv[1])
