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

from rpy import *
from hclust_parser import *
from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Util
from pyevolve import Selectors
from pyevolve import Crossovers
import sys
import threading 

class Temporary_output:
    def __init__(self):
        """This class is just implemented to receive temporarily the sys.stdout"""
        self.content = []
    def write(self, string):
        self.content.append(string)
    def flush(self):
        pass

class BioSig_console:
    def __init__(self, datafile ="filtered_data.txt", classes_file = "classe.txt", 
                 output_file = "filtered_data", generations = 100, population_size = 5, 
                 chromossome_lenght = 4039, ga_selector = "GRouletteWheel", optimization_method = "distance",  
                 mutation_ratio = 0.02, crossover_ratio = 0.8,
                 method_dist = "euclidean", method_hclust="average", n_bootstrap = 100, size_weight = 2, 
                 groups_weight = 1, bootstrap_weight = 2, distance_ratio_weight = 9, variation_weight = 1, 
                 parallel = False, biosignature_max_size = 50, biosignature_min_size = 10, 
                 biosignature_min_acceptable_size = 5, buffer_penalty = 1, display_output = 1,
                 F_factor = 0.2, ga_repetition = 1, n_extremes_values = 1):
        """
            BioSig is a program developed to assist Biomarker discovery from microarray or proteomic data.
            Given samples belonging to different groups, find those variables (expressed mRNAs or proteins, 
            for transcriptomic or proteomic analysis respectively) that, by a hierarchical clustering, 
            could recovery those groups. 
            The solution is optimized by a genetic algorithm procedure, implemented by pyevolve. 
            Inputs: 
                normalized and filtered data file name (argument: datafile); 
                file containing the groups of each sample (classes_file); 
                output filename (output_file);
                number of generations for the G.A. (generations); default: 100
                population size for the G.A. (population_size); default: 5
                chromossome lenght for the G.A. (chromossome_lenght); 
                selector (ga_selector)
                optimization method (optimization_method): 'distance' or 'clustering'; default: clustering; 
                mutation ratio (mutation_ratio); default: 0.8
                crossover ratio(crossover_ratio); default: 0.02
                distance method for hierachical clustering and distance calculation (method_dist);
                - for Hierarchical clustering
                    clustering method for hierachical clustering (method_hclust); default: average;
                    number of bootstrap pseudo replications for hierachical clustering (n_bootstrap); default: 100;
                weight for chromossome size ratio (size_weight), i.e., number of genes included (standardized);
                weight for bootstrap support (bootstrap_weight);
                weight for distance ratio (distance_ratio_weight);
                weight for F factor (F_factor), default: 0.2;
                - for Biosignature size               
                    biosignature_max_size, default: 50;
                    biosignature_min_size, default: 10;
                    biosignature_min_acceptable_size, default: 5;
                    buffer_penalty, default: 1;

                display_output: 1 if GA progression should be displayed
                n_extremes_values = 1
            Output:
                those mRNA/proteins that gave best solution
        """
        #Define optimization parameters
        self.display_output = display_output
        self._independent_best_solution = -1000
        self._independent_best_chromossome = [] 
        self.optimization_method = optimization_method
        self.method_dist = method_dist
        self.method_hclust = method_hclust
        self.n_bootstrap = n_bootstrap
        self.parallel = parallel
        #Define I/O
        self.classes_file = classes_file
        self.data_file = datafile
        self.output_file = output_file
        self.output = open(self.output_file + ".txt", "w")
        self.plot_output_filename = output_file[:-4] + ".png" 
        self.groups = self.parse_groups(self.classes_file)
        print "\nGroups and their Samples"
        for group, samples in self.groups.iteritems():
            print group, samples
        #Define GA parameter
        self.chromossome_lenght = chromossome_lenght
        self.population_size = population_size
        self.generations = generations
        self.ga_selector = ga_selector
        self.mutation_ratio = mutation_ratio
        self.crossover_ratio = crossover_ratio
        #Define weights
        self.size_weight = size_weight 
        self.groups_weight = groups_weight
        self.bootstrap_weight = bootstrap_weight
        self.distance_ratio_weight = distance_ratio_weight
        self.variation_weight = variation_weight
        self.F_factor = F_factor
        self.n_extremes_values = n_extremes_values
        #Define Bisignature size parameter
        self.biosigbiosignature_max_size = biosignature_max_size
        self.biosigbiosignature_min_size = biosignature_min_size
        self.biosigbiosignature_min_acceptable_size = biosignature_min_acceptable_size
        self.buffer_penalty = buffer_penalty
    	#The next lines record analysis's parameters
        self.output.write("I/O parameters\n")
        self.output.write("\tInput data: %s\n" %datafile)
        self.output.write("\tOutput data: %s\n" %output_file)
        self.output.write("\tGroups data: %s\n" %classes_file)
        self.output.write("Genetic Algorithms parameters\n")
        self.output.write("\tChromossome Lenght: %s\n" %str(chromossome_lenght))
        self.output.write("\tPopulation size: %s\n" %str(population_size))
        self.output.write("\tGenerations: %s\n" %str(generations))
        self.output.write("\tGA selector: %s\n" %str(ga_selector))
        self.output.write("\tMutation Ratio: %s\n" %str(mutation_ratio))
        self.output.write("\tCrossOver Ratio: %s\n" %str(crossover_ratio))
        self.output.write("Weights parameters\n")
        self.output.write("\tSize weight: %s\n" %str(size_weight))
        self.output.write("\tGroups weight: %s\n" %str(groups_weight))
        self.output.write("\tBootstrap weight: %s\n" %str(bootstrap_weight))
        self.output.write("\tDistance ratio weight: %s\n" %str(distance_ratio_weight))
        self.output.write("\tVariation weight: %s\n" %str(variation_weight))
        self.output.write("Biosignature size parameters\n")
        self.output.write("\tmax size: %s\n" %str(biosignature_max_size))
        self.output.write("\tmin size: %s\n" %str(biosignature_min_size))
        self.output.write("\tmin acceptable size: %s\n" %str(biosignature_min_acceptable_size))
        self.output.write("\tBuffer penalty: %s\n" %str(buffer_penalty))
        self.output.write("\tObjective Function parameters:")
        self.output.write("\tF-factor: %s\n" %str(F_factor))
        self.output.write("\tNumber of extreme values included: %s\n" %self.n_extremes_values)
        self.output.write("\tDistance measure: %s\n" %self.method_dist)

        # load data from microarrays using RPy
        r('library("pvclust")')
        r('data <- read.table("%s", header=TRUE, row.names=1)' %self.data_file)
        r('data_sd <- sd(t(data))')
        #Getting and setting samples indexes and colors
        self.matrix_labels = r('colnames(data)')  
        self.groups_indexes = {}
        self.samples_colors = []
        self.get_sample_indexes_and_colors()
        #Creating group_median_matrix
        #this matrix will contain the median of each sonde for each group 
        r('group_median_matrix <- matrix(nrow=nrow(data), ncol=%s)' %(str(len(self.groups.keys()))))
        r('dimnames(group_median_matrix)[[1]] <- labels(data)[[1]]')
        r('dimnames(group_median_matrix)[[2]] <- c(%s)'%str(self.groups.keys()).replace(" ","")[1:-1])
        #Creating within_group_dispersion_matrix
        #this matrix will contain the a dispersion measure of each sonde for each group 
        r('within_group_dispersion_matrix <- matrix(nrow=nrow(data), ncol=%s)' %(str(len(self.groups.keys()))))
        r('dimnames(within_group_dispersion_matrix)[[1]] <- labels(data)[[1]]')
        r('dimnames(within_group_dispersion_matrix)[[2]] <- c(%s)'%str(self.groups.keys()).replace(" ","")[1:-1])
        #Creating between_groups_dispersion_matrix
        #this matrix will contain the a dispersion measure of each sonde among all groups 
        r('between_groups_dispersion_matrix <- matrix(nrow=nrow(data), ncol=1)')
        r('dimnames(between_groups_dispersion_matrix)[[1]] <- labels(data)[[1]]')
        #Calculate group means and standard deviation
        print "\nGroups and their indexes"
        for group, indexes in self.groups_indexes.iteritems():         
            print group, indexes
            r('group_median_matrix[,"%s"] <- apply(data[,c(%s)],1,median)' %(group.strip(), str(indexes).replace(" ","")[1:-1]))
            r('within_group_dispersion_matrix[,"%s"] <- apply(data[,c(%s)],1,IQR)' %(group.strip(), str(indexes).replace(" ","")[1:-1]))            

        r('between_groups_dispersion_matrix <- apply(group_median_matrix,1,IQR)')            

        #Starting 
        print "Dimensions of loaded data: ", r('dim(data)')
        print "Microarray data loaded"
        print "Setting G.A. parameters"
        self.genome = G1DList.G1DList(self.chromossome_lenght)
        if self.optimization_method == "clusters":
            self.genome.evaluator.set(self.check_clusters)
            self.output.write("Clustering: %s; Bootstrap: %s; Chromossome lenght: %s; Population size: %s; Generations: %s \n" \
		                  %(method_hclust, str(n_bootstrap), str(chromossome_lenght), str(population_size), str(generations)))      
            self.output.write("Groups weight: %s; Size weight: %s; Bootstrap weight: %s\n" \
		                  %(str(self.groups_weight), str(self.size_weight), str(self.bootstrap_weight)))      

        elif self.optimization_method == "distance":
            self.genome.evaluator.set(self.check_distance)
            self.output.write("Chromossome lenght: %s; Population size: %s; Generations: %s \n" \
		                  %(str(chromossome_lenght), str(population_size), str(generations)))      
            self.output.write("Distance weight: %s; Size weight: %s\n" \
		                  %(str(self.distance_ratio_weight), str(self.size_weight)))      

        self.genome.setParams(rangemin=0.0, rangemax=1.0)
        self.genome.crossover.set(Crossovers.G1DBinaryStringXTwoPoint)
        self.ga = GSimpleGA.GSimpleGA(self.genome, interactiveMode=False)
        self.ga.setPopulationSize(self.population_size)
        self.ga.setGenerations(self.generations)
        self.ga.selector.set(Selectors.GRouletteWheel)
        self.ga.setMutationRate(self.mutation_ratio)
        self.ga.setCrossoverRate(self.crossover_ratio)
    ####################################### END #################################

    def get_sample_indexes_and_colors(self):
        groups_indexes = {}
        groups_colors = self.set_group_colors()     
        print "\nGetting indexes for samples"
        print self.matrix_labels
       
        for group, samples in self.groups.iteritems(): 
        #these loops create a dictionary with group names as keys and a list with the samples indexes 
        #in the distance matrix as values and a dictionary relating sample and group's color
            sample_indexes = []
            for label_index in range(len(self.matrix_labels)):
                if self.matrix_labels[label_index] in samples:
                    sample_indexes.append(label_index+1)
                    self.samples_colors.append(groups_colors[group])
                else:
                    pass
            
            self.groups_indexes[group] = sample_indexes
            print group, samples, groups_indexes

        return groups_indexes
    ####################################### END #################################

    def set_group_colors(self):
        """Set a color for each group for ploting the results"""
        colors = ["blue", "black", "brown", "magenta", "green", "orange", "red", 
                  "yellow", "green4", "tomato", "violet", "wheat", "red4", "purple4", 
                  "cadetblue1", "blue4", "chocalate"]
        groups = self.groups.keys()
        group_colors = {}
        for index in range(len(groups)): 
            group_colors[groups[index]] = colors[index]

        return group_colors
    ####################################### END #################################

    def run(self):
        print "Evolving the solution... this can take a long time"
        self.output.write("Following the output of BioSig\n")

        if self.optimization_method == "clusters":
            self.output.write("Expected groups, recognized groups, genes included, bootstrap values, objective function\n") 
            print "Expected groups, recognized groups, genes included, bootstrap values, objective function\n"

        elif self.optimization_method == "distance":
            self.output.write("Max within group, Min between group, Distance ratio, genes included, chromossome penalty, dispersion of groups (median of IQR), objective function\n") 
            print "Max within group, Min between group, Distance ratio, genes included, chromossome penalty, dispersion of groups (median of IQR), objective function\n"


        self.ga.evolve(freq_stats=1)
        self.best_solution = self.ga.bestIndividual()  
        self.output.write("The best chromossome:\n")
        self.output.write(self.parse_chromossome(self.best_solution))
        self.output.close()

        #Creating signature file of the best solution in the last generation
        data = open(self.data_file).readlines()
        signature_output = open(self.output_file + "_signature.csv", "w")
        signature_output.write(data[0])
        for linha in self.parse_chromossome(self.best_solution).split(","):
            signature_output.write(data[int(linha)])

        signature_output.close()

        #Creating signature file of the best solution found in all generations
        print "Best solution: ", self._independent_best_solution
        print "Best chromossome: ", self._independent_best_chromossome
        best_signature_output = open(self.output_file + "_best_signature.csv", "w")
        best_signature_output.write(data[0])
        for linha in self._independent_best_chromossome.split(","):
            best_signature_output.write(data[int(linha)])

        best_signature_output.close()

        print "Ploting best signature found in all generations"
        r('print(rownames(data[c(%s),]))'%self._independent_best_chromossome)
        r('png("%s_best_signature.png", unit="cm", width=30, height=30, res=200)' %self.output_file)
        try:
            r('result <- pvclust(data[c(%s),], method.dist="%s", method.hclust="%s", nboot=%s)' %(self._independent_best_chromossome, self.method_dist, self.method_hclust, str(500)))
            r('plot(result)')
            r('dev.off()')
        except RException:
            print RException
        
        print "Performing Multidimensional Scaling of the best signature"
        r.assign('mds_colors',self.samples_colors)
        r('dist_matrix <- dist(t(data[c(%s),]))'%self._independent_best_chromossome)
        r('fit <- cmdscale(dist_matrix,eig=TRUE, k=2)')
        r('png("%s_best_mds.png", unit="cm", width=40, height=20, res=200)' %self.output_file)
        r('plot(fit$points[,1], fit$points[,2], main="Multidimensional scaling of the best signature", col=mds_colors)')
        r('dev.off()')

        print "Ploting best signature found at the last generation"
        #Clustering signature
        Util.set_normal_term() 
        try:
            self.clustering(self.best_solution, True)
        except RException:
            print RException

        print "Performing Multidimensional Scaling of the GA signature"
        r('dist_matrix <- dist(t(data[c(%s),]))'%self.parse_chromossome(self.best_solution))
        r('fit <- cmdscale(dist_matrix,eig=TRUE, k=2)')
        r('png("%s_ga_mds.png", unit="cm", width=40, height=20, res=200)' %self.output_file)
        r('plot(fit$points[,1], fit$points[,2], main="Multidimensional scaling of the GA signature", col=mds_colors)')
        r('dev.off()')



        print "Optimization is finished"
        r('rm(list = ls())')

    ####################################### END ##################################

    def parse_chromossome(self,chromossome):
        """
        input: a binary chromossome in pyevolve format 
        eturns: a list with indexes from '1', 
        example: input [1,1,0,0,1], returns [0,1,4]
        """
        genes = ""
        for bit_index in range(len(chromossome)):
            if bit_index == len(chromossome)-1:
                if chromossome[bit_index] == 1:
                    genes += str(bit_index+1)
                else:
                    pass
                    genes = genes[:-1]
            else:
                if chromossome[bit_index] == 1:
                    genes += str(bit_index+1) + ","
                else:
                    pass

        return genes      
    ####################################### END #################################

    def multiply(self,list_to_be_multiplied):
        """Given a list of numbers, return the product of them"""
        value = 1
        for x in list_to_be_multiplied:
            value *= x
    
        return value
    ####################################### END #################################

    def parse_groups(self,groups_filename):
        """
        Given a file with each line containing 'Sample name, group it belongs', return a dictionary with like:
        {Group_1: [Sample_1, Sample_2], Group_n:[Sample_x, Sample_y]}
        """
        raw_data_file = open(groups_filename) 
        raw_data = raw_data_file.readlines()
        raw_data_file.close()
        data = [line.strip().split(",") for line in raw_data]
        groups = {}
        for datum in data:
            if datum[1] in groups.keys():
                groups[datum[1]].append(datum[0])
            else:
                groups[datum[1]] = [datum[0]]

        return groups
    ####################################### END #################################

    def clustering(self,chromossome, plot = False, hclust= False):
        """
        Inputs: chromossome in pyevolve format, 
                plot: True ou False for ploting the result
                hclust: True ou False to use hclust instead of pvclust
        
        return: return the cluster tree in pvclust format. If passed True, plot it graphically
        """   
        parsed_chromossome = self.parse_chromossome(chromossome)
        teste = parsed_chromossome.strip().split(",")
        list_parsed_chromossome = [int(x) for x in teste]
        r.assign('chromossome',list_parsed_chromossome)
        r('print(rownames(data[chromossome,]))')
        temp_output = Temporary_output()
        sys.stdout = temp_output

        if hclust:
            r('result <- hclust(dist(t(data[chromossome,])))')
            result_rpy = r("result")
            sys.stdout = sys.__stdout__
            return result_rpy

        if self.parallel:
            r('result <- parPvclust(data[chromossome,], method.dist="%s", method.hclust="%s", nboot=%s)' %(self.method_dist, self.method_hclust, str(self.n_bootstrap)))
        else:
            if plot:
                r('png("%s_ga_signature.png", unit="cm", width=30, height=30, res=200)' %self.output_file)
                r('result <- pvclust(data[chromossome,], method.dist="%s", method.hclust="%s", nboot=%s)' %(self.method_dist, self.method_hclust, str(500)))
                r('plot(result)')
                r('dev.off()')
            else:
                r('result <- pvclust(data[chromossome,], method.dist="%s", method.hclust="%s", nboot=%s)' %(self.method_dist, self.method_hclust, str(self.n_bootstrap)))

        result_rpy = r("result")
        sys.stdout = sys.__stdout__
        return result_rpy
    ####################################### END #################################

    def check_chromossome_validity(self, parsed_chromossome):
        testing_chromossome = parsed_chromossome.strip().split(",")
        try:
            list_parsed_chromossome = [int(x) for x in testing_chromossome]
            return list_parsed_chromossome
        except ValueError:
            print "ValueError"
            message_console = "wm:0, bm:0, wd:0, bd:0, dr:0, chr:0(0), sd:0, of:0\n"   
            message_file_output = "0, 0, 0, 0, 0, 0, 0\n" 
            if self.display_output:
                print message_console
            else:
                pass 

            self.output.write(message_file_output)
            return 0
    ####################################### END #################################


    def calculate_chromossome_lenght_penalty(self, chromossome):
        """Given a chromossome, calculates its size penalty"""
        chromossome_len = float(len(chromossome))
        chromossome_sum = float(sum(list(chromossome)))
        number_of_groups = float(len(self.groups_indexes.keys()))
        chromossome_ratio = 0.0

        if chromossome_sum < self.biosigbiosignature_min_acceptable_size:
            chromossome_ratio = 1 
            if self.display_output:
                print "unacceptable chromossome size"

        elif chromossome_sum >= self.biosigbiosignature_min_acceptable_size and chromossome_sum < self.biosigbiosignature_min_size:
            chromossome_ratio = (self.buffer_penalty*((self.biosigbiosignature_min_size+1)- chromossome_sum))/chromossome_len
            if self.display_output:
                print "chromossome size too small(buffer zone)"

        elif chromossome_sum >= self.biosigbiosignature_min_size and chromossome_sum <= self.biosigbiosignature_max_size:
            chromossome_ratio = 0#0.01
            if self.display_output:
                print "ideal chromosssome size"

        else:
            chromossome_ratio = chromossome_sum/chromossome_len

        return chromossome_ratio
    ####################################### END #################################

    def check_clusters(self,chromossome):
        """
        input: chromossome in pyevolve format
        return: the value of the objetive function that is.... in development
        """
        ############## Calculate chromossome size penalties ############
        chromossome_len = float(len(chromossome))
        chromossome_sum = float(sum(list(chromossome)))
        chromossome_ratio = self.calculate_chromossome_lenght_penalty(chromossome)
        ###########

        results = []
        bootstrap_values = []
        #here, execute pvclust with data from the "chromossome"
        cluster = pvclust_tree(self.clustering(chromossome))
        #check if it clustered the expected groups
        #if clustered, save its bootstrap support, if not, save 0
        for group, samples in self.groups.iteritems():
            result = cluster.tree.check_node_existance(samples)
            if result[0]:
                self.output.write("Found group %s \n" %group)
                print "Found group %s" %group
                results.append(True)
                bootstrap_values.append(result[2])
            else:
                results.append(False)
                bootstrap_values.append(0.0)

        # calculate the objective function
   
        # com poucos grupos e possivel otimizar os genes apenas utilizando aquelas topologia que apresentem todos os agrupamentos
        # entretanto, conforme o numero de grupos aumenta, a probabilidade de obter uma topologia com todos os grupo diminui
        # assim eh preciso considerar aquelas topologia que apresentem alguns grupos, mas nao todos... e entao ir otimizando
        # metodologia utilzada foi:
        # se todos os agrupamentos esperados forem achados, multiplica o valor de bootstrap de cada um deles, e entao multiplica pelo numero de clados achados
        # essa funcao objetiva otimiza duas variaveis: o numero de grupos recuperado e o valor de suporte deles
        # a funcao entao e a soma do numero de grupos recuperados com a multiplicacao do suporte de cada grupo
        bootstraps_to_be_multiplied = []            
        for cluster_index in range(len(results)):
           if results[cluster_index]:
               bootstraps_to_be_multiplied.append(bootstrap_values[cluster_index])
           else:
               pass
        message = "%s, %s, %s, %s, %s\n" %(str(len(results)), 
                  str(len(bootstraps_to_be_multiplied)), 
                  str(chromossome_sum), self.multiply(bootstraps_to_be_multiplied), str(self.bootstrap_weight*self.multiply(bootstraps_to_be_multiplied) + self.groups_weight*len(bootstraps_to_be_multiplied) - self.size_weight*chromossome_ratio))
        print message
        self.output.write(message)

        return self.groups_weight*len(bootstraps_to_be_multiplied) + self.bootstrap_weight*self.multiply(bootstraps_to_be_multiplied) - self.size_weight*chromossome_ratio
    ####################################### END #################################

    def check_distance(self,chromossome):                             
        """This function tries to minimize the distances among members of the same group and the number of genes included
           while try to maximize the distances between members of different groups and the mean variation among them"""   

        ############## Calculate chromossome size penalties ############
        chromossome_len = float(len(chromossome))
        chromossome_sum = float(sum(list(chromossome)))
        chromossome_ratio = self.calculate_chromossome_lenght_penalty(chromossome)

        ############## Creating distances matrix with R ############
        parsed_chromossome = self.parse_chromossome(chromossome)
        list_parsed_chromossome = self.check_chromossome_validity(parsed_chromossome)
        #This 'if' prevents checking a null signature
        #and prevents an unexpected interruption 
        if not list_parsed_chromossome:
            return 0

        r.assign('chromossome',list_parsed_chromossome)
        temp_output = Temporary_output()
        sys.stdout = temp_output
        matrix_labels = []
        if self.method_dist == "correlation":
            r('distance <- 1-cor(data[chromossome,])')
            r('dist_matrix <- as.matrix(distance)')
            matrix_labels = r("labels(distance)")[1]
        else:
            r('distance <- dist(t(data[chromossome,]))')
            r('dist_matrix <- as.matrix(distance)')
            matrix_labels = r("labels(distance)")


        sys.stdout = sys.__stdout__
        #Creating matrix with different distances
        groups_indexes = self.groups_indexes
        r('within_group_distances_matrix <- dist_matrix')
        r('between_group_distances_matrix <- dist_matrix')
        samples = 0
        within_measures = 0
        for group, indexes in groups_indexes.iteritems():
            samples += len(indexes)
            within_measures += (len(indexes)*(len(indexes)-1))/2
            r.assign('indexes', indexes)
            r('between_group_distances_matrix[indexes,indexes] <- NA')            
            r('indexes_r <- setdiff(c(1:dim(dist_matrix)[2]), indexes)')
            r('within_group_distances_matrix[indexes,indexes_r] <- NA')            

        between_measures = (samples*(samples-1)/2) - within_measures
        within_sum = r('sum(within_group_distances_matrix, na.rm = TRUE)')
        between_sum = r('sum(between_group_distances_matrix, na.rm = TRUE)')
        between_distances_mean = between_sum/between_measures
        within_distances_mean = within_sum/within_measures
        within_max = r('max(within_group_distances_matrix, na.rm = TRUE)')
        between_min = r('min(between_group_distances_matrix, na.rm = TRUE)')
        distance_ratio = 0
        if self.n_extremes_values == 1:
            distance_ratio =  -(within_max + self.F_factor*within_distances_mean)/(between_min + self.F_factor*between_distances_mean)
        else:
            within_highest_values = r('sum(sort(as.vector(within_group_distances_matrix), decreasing=TRUE, na.last=NA)[1:%s])' %str(self.n_extremes_values))
            between_lowest_values = r('sum(sort(as.vector(between_group_distances_matrix), decreasing=FALSE, na.last=NA)[1:%s])'%str(self.n_extremes_values))
            distance_ratio =  -(within_highest_values + self.F_factor*within_distances_mean)/(between_lowest_values + self.F_factor*between_distances_mean)
        
        #distance_ratio = -(within_max/between_min) - 0.1*(within_distances_mean/between_distances_mean)
        #distance_ratio = -(within_max/between_min)

        ################ Calculating standard deviation ###########################        
        standard_deviation = r('median(between_groups_dispersion_matrix[chromossome])')
        ################ Calculating objective_function ###########################
        objective_function = (1000 + self.distance_ratio_weight*distance_ratio + self.variation_weight*standard_deviation - self.size_weight*chromossome_ratio)
        ################ Reporting result ###########################
        message_console = "wm:%s, bm:%s, wd:%s, bd:%s, dr:%s, chr:%s(%s), sd:%s, of:%s\n" \
                                               %(str(within_max),str(between_min), 
                                               str(within_distances_mean), str(between_distances_mean),
                                               str(distance_ratio), str(chromossome_sum), str(chromossome_ratio), 
                                               str(standard_deviation), str(objective_function))

        message_file_output = "%s, %s, %s, %s, %s, %s, %s\n" %(str(within_max),str(between_min), 
                                                           str(distance_ratio), str(chromossome_sum), 
                                                           str(chromossome_ratio), str(standard_deviation), 
                                                           str(objective_function))
        if self.display_output:
            print message_console
        
        self.output.write(message_file_output)

        #GA returns the best solution in the last generation. 
        #In order to have the best solution of all optimization,  keep and update de the best solution 
        #found in all generations
        if distance_ratio > self._independent_best_solution:
            self._independent_best_solution = distance_ratio
            self._independent_best_chromossome = parsed_chromossome
        else:
            pass

        #Assert that objective function is positive
        if objective_function < 0:
            print "Some error have ocurred. Objective function can not be negative."
            return 0
        else:
            return objective_function

class BioSig(BioSig_console, threading.Thread):
#        threading.Thread.__init__(self)    

    def __init__(self, datafile ="filtered_data.txt", classes_file = "classe.txt", 
                 output_file = "filtered_data", generations = 100, population_size = 5, 
                 chromossome_lenght = 4039, ga_selector = "GRouletteWheel", optimization_method = "distance",  
                 mutation_ratio = 0.02, crossover_ratio = 0.8,
                 method_dist = "euclidean", method_hclust="average", n_bootstrap = 100, size_weight = 2, 
                 groups_weight = 1, bootstrap_weight = 2, distance_ratio_weight = 9, variation_weight = 1, 
                 parallel = False, biosignature_max_size = 50, biosignature_min_size = 10, 
                 biosignature_min_acceptable_size = 5, buffer_penalty = 1, display_output = 1,
                 F_factor = 0.2, ga_repetition = 1, n_extremes_values = 1):
        """
            BioSig is a program developed to assist Biomarker discovery from microarray or proteomic data.
            Given samples belonging to different groups, find those variables (expressed mRNAs or proteins, 
            for transcriptomic or proteomic analysis respectively) that, by a hierarchical clustering, 
            could recovery those groups. 
            The solution is optimized by a genetic algorithm procedure, implemented by pyevolve. 
            Inputs: 
                normalized and filtered data file name (argument: datafile); 
                file containing the groups of each sample (classes_file); 
                output filename (output_file);
                number of generations for the G.A. (generations); default: 100
                population size for the G.A. (population_size); default: 5
                chromossome lenght for the G.A. (chromossome_lenght); 
                selector (ga_selector)
                optimization method (optimization_method): 'distance' or 'clustering'; default: clustering; 
                mutation ratio (mutation_ratio); default: 0.8
                crossover ratio(crossover_ratio); default: 0.02
                distance method for hierachical clustering and distance calculation (method_dist);
                - for Hierarchical clustering
                    clustering method for hierachical clustering (method_hclust); default: average;
                    number of bootstrap pseudo replications for hierachical clustering (n_bootstrap); default: 100;
                weight for chromossome size ratio (size_weight), i.e., number of genes included (standardized);
                weight for bootstrap support (bootstrap_weight);
                weight for distance ratio (distance_ratio_weight);
                weight for F factor (F_factor), default: 0.2;
                - for Biosignature size               
                    biosignature_max_size, default: 50;
                    biosignature_min_size, default: 10;
                    biosignature_min_acceptable_size, default: 5;
                    buffer_penalty, default: 1;

                display_output: 1 if GA progression should be displayed

            Output:
                those mRNA/proteins that gave best solution
        """
        threading.Thread.__init__(self)    
    
        #Define optimization parameters
        self.display_output = display_output
        self._independent_best_solution = -1000
        self._independent_best_chromossome = [] 
        self.optimization_method = optimization_method
        self.method_dist = method_dist
        self.method_hclust = method_hclust
        self.n_bootstrap = n_bootstrap
        self.parallel = parallel
        #Define I/O
        self.classes_file = classes_file
        self.data_file = datafile
        self.output_file = output_file
        self.output = open(self.output_file + ".txt", "w")
        self.plot_output_filename = output_file[:-4] + ".png" 
        self.groups = self.parse_groups(self.classes_file)
        print "\nGroups and their Samples"
        for group, samples in self.groups.iteritems():
            print group, samples
        #Define GA parameter
        self.chromossome_lenght = chromossome_lenght
        self.population_size = population_size
        self.generations = generations
        self.ga_selector = ga_selector
        self.mutation_ratio = mutation_ratio
        self.crossover_ratio = crossover_ratio
        #Define weights
        self.size_weight = size_weight 
        self.groups_weight = groups_weight
        self.bootstrap_weight = bootstrap_weight
        self.distance_ratio_weight = distance_ratio_weight
        self.variation_weight = variation_weight
        self.F_factor = F_factor
        self.n_extremes_values = n_extremes_values
        #Define Bisignature size parameter
        self.biosigbiosignature_max_size = biosignature_max_size
        self.biosigbiosignature_min_size = biosignature_min_size
        self.biosigbiosignature_min_acceptable_size = biosignature_min_acceptable_size
        self.buffer_penalty = buffer_penalty
    	#The next lines record analysis's parameters
        self.output.write("I/O parameters\n")
        self.output.write("\tInput data: %s\n" %datafile)
        self.output.write("\tOutput data: %s\n" %output_file)
        self.output.write("\tGroups data: %s\n" %classes_file)
        self.output.write("Genetic Algorithms parameters\n")
        self.output.write("\tChromossome Lenght: %s\n" %str(chromossome_lenght))
        self.output.write("\tPopulation size: %s\n" %str(population_size))
        self.output.write("\tGenerations: %s\n" %str(generations))
        self.output.write("\tGA selector: %s\n" %str(ga_selector))
        self.output.write("\tMutation Ratio: %s\n" %str(mutation_ratio))
        self.output.write("\tCrossOver Ratio: %s\n" %str(crossover_ratio))
        self.output.write("Weights parameters\n")
        self.output.write("\tSize weight: %s\n" %str(size_weight))
        self.output.write("\tGroups weight: %s\n" %str(groups_weight))
        self.output.write("\tBootstrap weight: %s\n" %str(bootstrap_weight))
        self.output.write("\tDistance ratio weight: %s\n" %str(distance_ratio_weight))
        self.output.write("\tVariation weight: %s\n" %str(variation_weight))
        self.output.write("Biosignature size parameters\n")
        self.output.write("\tmax size: %s\n" %str(biosignature_max_size))
        self.output.write("\tmin size: %s\n" %str(biosignature_min_size))
        self.output.write("\tmin acceptable size: %s\n" %str(biosignature_min_acceptable_size))
        self.output.write("\tBuffer penalty: %s\n" %str(buffer_penalty))
        self.output.write("\tObjective Function parameters:")
        self.output.write("\tF-factor: %s\n" %str(F_factor))
        self.output.write("\tNumber of extreme values included: %s\n" %self.n_extremes_values)
        self.output.write("\tDistance measure: %s\n" %self.method_dist)
        # load data from microarrays using RPy
        r('library("pvclust")')
        r('data <- read.table("%s", header=TRUE, row.names=1)' %self.data_file)
        r('data_sd <- sd(t(data))')
        #Getting and setting samples indexes and colors
        self.matrix_labels = r('colnames(data)')  
        self.groups_indexes = {}
        self.samples_colors = []
        self.get_sample_indexes_and_colors()
        #Creating group_median_matrix
        #this matrix will contain the median of each sonde for each group 
        r('group_median_matrix <- matrix(nrow=nrow(data), ncol=%s)' %(str(len(self.groups.keys()))))
        r('dimnames(group_median_matrix)[[1]] <- labels(data)[[1]]')
        r('dimnames(group_median_matrix)[[2]] <- c(%s)'%str(self.groups.keys()).replace(" ","")[1:-1])
        #Creating within_group_dispersion_matrix
        #this matrix will contain the a dispersion measure of each sonde for each group 
        r('within_group_dispersion_matrix <- matrix(nrow=nrow(data), ncol=%s)' %(str(len(self.groups.keys()))))
        r('dimnames(within_group_dispersion_matrix)[[1]] <- labels(data)[[1]]')
        r('dimnames(within_group_dispersion_matrix)[[2]] <- c(%s)'%str(self.groups.keys()).replace(" ","")[1:-1])
        #Creating between_groups_dispersion_matrix
        #this matrix will contain the a dispersion measure of each sonde among all groups 
        r('between_groups_dispersion_matrix <- matrix(nrow=nrow(data), ncol=1)')
        r('dimnames(between_groups_dispersion_matrix)[[1]] <- labels(data)[[1]]')
        #Calculate group means and standard deviation
        print "\nGroups and their indexes"
        for group, indexes in self.groups_indexes.iteritems():         
            print group, indexes
            r('group_median_matrix[,"%s"] <- apply(data[,c(%s)],1,median)' %(group.strip(), str(indexes).replace(" ","")[1:-1]))
            r('within_group_dispersion_matrix[,"%s"] <- apply(data[,c(%s)],1,IQR)' %(group.strip(), str(indexes).replace(" ","")[1:-1]))            

        r('between_groups_dispersion_matrix <- apply(group_median_matrix,1,IQR)')            

        #Starting 
        print "Dimensions of loaded data: ", r('dim(data)')
        print "Microarray data loaded"
        print "Setting G.A. parameters"
        self.genome = G1DList.G1DList(self.chromossome_lenght)
        if self.optimization_method == "clusters":
            self.genome.evaluator.set(self.check_clusters)
            self.output.write("Clustering: %s; Bootstrap: %s; Chromossome lenght: %s; Population size: %s; Generations: %s \n" \
		                  %(method_hclust, str(n_bootstrap), str(chromossome_lenght), str(population_size), str(generations)))      
            self.output.write("Groups weight: %s; Size weight: %s; Bootstrap weight: %s\n" \
		                  %(str(self.groups_weight), str(self.size_weight), str(self.bootstrap_weight)))      

        elif self.optimization_method == "distance":
            self.genome.evaluator.set(self.check_distance)
            self.output.write("Chromossome lenght: %s; Population size: %s; Generations: %s \n" \
		                  %(str(chromossome_lenght), str(population_size), str(generations)))      
            self.output.write("Distance weight: %s; Size weight: %s\n" \
		                  %(str(self.distance_ratio_weight), str(self.size_weight)))      

        self.genome.setParams(rangemin=0.0, rangemax=1.0)
        self.genome.crossover.set(Crossovers.G1DBinaryStringXTwoPoint)
        self.ga = GSimpleGA.GSimpleGA(self.genome, interactiveMode=False)
        self.ga.setPopulationSize(self.population_size)
        self.ga.setGenerations(self.generations)
        self.ga.selector.set(Selectors.GRouletteWheel)
        self.ga.setMutationRate(self.mutation_ratio)
        self.ga.setCrossoverRate(self.crossover_ratio)
    ####################################### END #################################



        
if __name__ == "__main__":
    a = BioSig()
           
