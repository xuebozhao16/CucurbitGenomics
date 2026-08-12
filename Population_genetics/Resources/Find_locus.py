

HOMOZYGOTE_1 = 1
HOMOZYGOTE_2 = 2
HETEROZYGOTE = 0
MISSING_DATA = -1

#Returning the allele frequency of each sample
def get_allele_type(alleles):
    if (alleles[1] == alleles[3]): #Homozygot sample
        if (alleles[1] == '1'): #Homoz. for allele 1 (0101)
            return HOMOZYGOTE_1
        if (alleles[1] == '0'): #Data missing (0000)
            return MISSING_DATA
        return HOMOZYGOTE_2 #Homoz. for allele 2 (0202 or 0303)
    else:
        return HETEROZYGOTE #Heteroz. (0103 or 0102)

def get_allele_freq(allele1, allele2):
    if (allele1 == 0 and allele2 == 0):
        #print "#######################"
        #print allele1
        #print allele2
        #print "#######################"
        return 99
    else:
        return float(allele1) / (allele1 + allele2)

################## Class Locus ######################
class Locus:
    
    def __init__(self, name, pops):
        
        self.name = name.strip()
        self.pop = [[0 for x in xrange(2)] for x in xrange(pops)]
        self.freqs = [[0 for x in xrange(1)] for x in xrange(pops)]
        self.is_monomorphic = False #The locus is monomorphic for one allele.
        self.is_consensus = False #True if freq differences are pops < cut off.
        self.freq_diff = 0 #Difference in allele freq among pop.

    #Updates the allele frequencies for each population.
    def update_freqs(self, al_type, n):
        if (al_type == HOMOZYGOTE_1):
            self.pop[n][0] += 2
        elif (al_type == HOMOZYGOTE_2):
            self.pop[n][1] += 2
        elif (al_type == MISSING_DATA):
            pass
        else:
            self.pop[n][0] += 1
            self.pop[n][1] += 1
        #print self.pop

    def get_name(self):
        return self.name

    #Prints the name and allele frequency for each locus
    def to_string(self):
        if (self.is_monomorphic_locus()):
            print "Removing monomorphic locus " + self.name + " with the following allele counts:"
            for i in range(0, len(self.pop)):
                curr_pop = self.pop[i]
                print curr_pop

    #Adding the allele freqs to bayenv format (null model)   
    def freqs_to_lines(self):
        line_1 = ""
        line_2 = ""
        if (not self.is_monomorphic_locus()):
            for i in range(0, len(self.pop)):
                line_1 += str(self.pop[i][0]) + "\t"
                line_2 += str(self.pop[i][1]) + "\t"
            line_1 += "\n" + line_2 + "\n"
        return line_1

        #Adding the allele freqs to bayenv format for test purposes
    def freqs_to_lines2(self):
        line_1 = ""
        line_2 = ""
        if (not self.is_monomorphic_locus()):
            line_1 += self.name.replace(":","-") + "\t"
            for i in range(0, len(self.pop)):
                line_1 += str(self.pop[i][0]) + "\t"
                line_2 += str(self.pop[i][1]) + "\t"
            line_1 += "\n" + line_2 + "\n"
        return line_1
    
##      def to_string(self):
##         if (not self.equal_freqs_among_pops()):
##             print self.name
##             for i in range(0, len(self.pop)):
##                 curr_pop = self.pop[i]
##                 print curr_pop


    #Returns true if all pops are zero in one of the alleles.
    def is_monomorphic_locus(self):
        prev_freq = None
        if (self.pop[0][0] == 0): #First allele is zero
            for i in range(0, len(self.pop)):
                curr_freq = self.pop[i][0]
                if (prev_freq is not None and curr_freq != 0):
                    return self.is_monomorphic
                prev_freq = curr_freq
            self.is_monomorphic = True
        elif (self.pop[0][1] == 0): #Second allele is zero
            for i in range(0, len(self.pop)):
                curr_freq = self.pop[i][1]
                if (prev_freq is not None and curr_freq != 0):
                    return self.is_monomorphic
                prev_freq = curr_freq
            self.is_monomorphic = True
        else:
            pass
        return self.is_monomorphic

    #Returns true if the allele frequecy is less than x.
    def equal_freqs_among_pops(self):
        min_diff = 0.01
        prev_freq = None
        is_equal_freq = True
        for i in range(0, len(self.pop)):
            curr_pop = self.pop[i]
            curr_freq = get_allele_freq(curr_pop[0], curr_pop[1])
            if (prev_freq != None):
                diff = abs(curr_freq-prev_freq)
            else:
                diff = 0
            if (prev_freq is not None and diff > min_diff):
                is_equal_freq = False
                break
            prev_freq = curr_freq
        return is_equal_freq
    
    #Calculates the allele freqs for the locus
    def set_freqs(self):
        for i in range(0, len(self.pop)):
            curr_pop = self.pop[i]
            #print self.name
            ####
            self.freqs[i] = get_allele_freq(curr_pop[0], curr_pop[1])
            if self.freqs[i] == 99: #Marking the SNP for deletion
                self.is_consensus = True
                self.is_monomorphic = True
                
            
    #Returns the frequency array for the population
    def get_freqs(self):
        return self.freqs


    #Returns true if the allele frequecy is less than x.
    def is_consensus_locus(self, min_diff):
        is_consensus = True
        for i in range(1, len(self.pop)):
            #print self.freqs
            diff = abs(self.freqs[0]-self.freqs[i])
            if (diff > min_diff):
                is_consensus = False
                break
        return is_consensus

    #Sets the MAFD
    def set_max_freq_diff(self,):
        max_diff = 0
        k = 0
        while (k < len(self.pop)-1):
            for i in range(k+1, len(self.pop)):
                diff = abs(self.freqs[k]-self.freqs[i])
                if (diff > max_diff):
                    max_diff = diff
                    self.freq_diff = diff
            k += 1
        #print "Maximum difference is: " + str(self.freq_diff)

    def get_max_freq_diff(self):
        return self.freq_diff


