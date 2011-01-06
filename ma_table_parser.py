#    Copyright 2010 Frederico G. C. Arnoldi <fgcarnoldi /at/ gmail /dot/ com>
#
#    This file is part of pyBioSig.
#
#    pyBioSig is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    pyBioSig is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with pyBioSig.  If not, see <http://www.gnu.org/licenses/>.


#import numpy
#import time
from rpy import *
import threading 

class MAGE_TAB_data_matrix(threading.Thread):
    """
    A class for modeling a MAGE_TAB matrix
    Input: 
        matrix_filename
        group_filename
        output filename
        method for filtering data: Kruskal-Wallis or Anova
        p-value threshold
    Functions:
        get_table_header
        load_groups
        get_R_factors
        run: filter informative genes using the especified method and p-value threshold
    """

    def __init__(self,matrix_filename, group_filename, output = "filtered_data.txt", method = "kruskal-wallis", p_value_threshold = 0.05):
        threading.Thread.__init__(self)    
        self.matrix_filename = matrix_filename
        self.header = self.get_table_header(self.matrix_filename)
        self.groups_filename = group_filename
        self.output = output
        print "Output at: ", self.output
        self.filter_method = method
        self.p_value_threshold = p_value_threshold
        print "Inespecific filtering method: ", self.filter_method, " with P < ", str(self.p_value_threshold)
        self.groups = {}
        print "loading groups"
        self.load_groups()
        print self.groups, len(self.groups)
        print "setting factors"
        self.r_factors = self.get_R_factors()
        print self.r_factors, len(self.r_factors)

    def get_table_header(self, filename):
        """Given a file with a microarray data matrix, return its header names"""
        print "getting table's header"
        raw_data = open(filename).readlines()
        data_temp = [x.split("\t") for x in raw_data]
        if len(data_temp[0]) == 1:
            data_temp[0] = data_temp[0].split(" ")
        labels = data_temp[0]
        labels[-1] = labels[-1].replace("\r\n", "") 
        return labels

    
    def load_groups(self):
        """
        Given a file with each line containing 'Sample name, group it belongs' (class file passed to __init__),
        update self.groups as following:
        {Group_1: [Sample_1, Sample_2], Group_n:[Sample_x, Sample_y]}
        """
        raw_data = open(self.groups_filename).readlines()
        data = [line.strip().split(",") for line in raw_data]
        for datum in data:
            if datum[1].strip('"') in self.groups.keys():
               self.groups[datum[1]].append(datum[0])
            else:
                self.groups[datum[1]] = [datum[0]]


    def get_R_factors(self):
        """Return a string to use with r() to set 'factors'"""
        labels = self.header
        print "labels"
        print labels

        r_factors = "groups <- c("
        #adiciona primeira amostra, sem virgula antes
        print "Assigning samples to groups: "
        for group_name, members in self.groups.iteritems():
             if labels[1] in members:
                    r_factors += "'" + group_name + "'" 
                    print labels[1], group_name, "1"
        #adiciona demais amostras, com virgula antes
        i = 2
        for label in labels[2:]:      
            label = label.strip("\n")
            for group_name, members in self.groups.iteritems():
                if label in members:
                    r_factors += ", '" + group_name + "'" 
                    print label, group_name, i
                    i += 1
                    
                else:
                    pass

        r_factors += ")"
        print "\n************ R factors **************"
        print r_factors
        return r_factors


    def run(self):
        """Run unspecific filtering""" 
        r('data <- read.table("%s", dec=".", header=TRUE, row.names=1)'%self.matrix_filename)
        r('data <- as.matrix(data)')
        r('%s' %self.r_factors)
        r('groups <- factor(groups)')
        print "Checking loaded data before filtering"
        print "\nSample's group, recognized by R: ", r('groups')
        print "\nNumber of groups recognized by R: ", r('length(groups)')
        print "\nData dimensions loaded in R: ", r('dim(data)')
        print "\nFirst line of loaded data (5 cols): ", r('data[1,1:5]')
        print "\nColumn names: ", r('colnames(data)')
        print "\nNumber of columns: ", r('length(colnames(data))')
        print "\n\nUnspecific filtering can take some minutes. Please, wait for the ending message."
        filter_function = ""
        if self.filter_method == "kruskal-Wallis":
            filter_function = 'filter_function <- function(x){kruskal.test(x~groups)$p.value}' 
        elif self.filter_method == "anova":
            filter_function = 'filter_function <- function(x){anova(lm(x~groups))[[5]][1]}'

        r('%s' %filter_function)
        r('raw_pvalues <- apply(data, 1, filter_function)')
        if self.filter_method == "kruskal-Wallis":
            r('raw_pvalues[is.na(raw_pvalues)] <- 1') 

        r('positives <- raw_pvalues < %s' %self.p_value_threshold)
        r('row_names <- rownames(data)[positives]')
        r('col_names <- c(c("Sonde"), colnames(data))')
        r('print(raw_pvalues)')
        data_dim = r('dim(data)')
        selected_dim = r('length(row_names)')
        print "Positive genes: %s of %s => %s" %(str(selected_dim), str(data_dim[0]), str(float(selected_dim)/data_dim[0])) 
        r('write.table(data[positives,], "%s", row.names=row_names, col.names=TRUE, quote=FALSE, sep="\t")' %self.output)
        #The code starting here ....
        mt_content = "Sondes\t" + open(self.output).read()
        table_output = open(self.output, "w")
        table_output.write(mt_content)
        table_output.close()
        #...and end here, can be optimized. However, I am not aware how to do it directly in R
        print "Unspecific filtering completed!"
        

if __name__ == "__main__":
    print "Nothing yet"
    
