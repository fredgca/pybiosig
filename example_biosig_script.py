# This file is part of pyBioSig.
# It is an example of how to use pyBioSig using python scripts for batch analyses


from BioSig import *
import time
class batch:
    def __init__(self, data_file, classes_file, output, n_analyses):
        self.output = output
        self.n_analyses = n_analyses
        print "Initializing parameters....."
        self.parameters = {}
        ###I/O data
        self.parameters["datafile"] = data_file
        self.parameters["classes_file"] = classes_file
        self.PATH = self.parameters["datafile"][:self.parameters["datafile"].rfind("/")+1]
        #Evolve soluion using distance 
        self.parameters["optimization_method"] = "distance"
        ###Biosignature size parameters
        self.parameters["buffer_penalty"] = 10
        self.parameters["biosignature_max_size"] = 10
        self.parameters["biosignature_min_size"] = 6
        self.parameters["biosignature_min_acceptable_size"] = 4
        ###Weight values
        self.parameters["size_weight"] = 20
        self.parameters["distance_ratio_weight"] = 10
        self.parameters["variation_weight"] = 0
        ###Distance parmeters
        self.parameters["method_dist"] = "euclidean"
        ###GA parameters
        self.parameters["population_size"] = 4
        self.parameters["generations"] = 25000
        self.parameters["chromossome_lenght"] = 432
        self.parameters["ga_selector"] = "GRouletteWheel"
        self.parameters["display_output"] = True
        #self.parameters["from_interface"] = 0

    def run(self):
        ## Send them to the pyBiomarkerFinder  
        output_time_file = self.output + "time.txt"
        output_time = open(output_time_file, "w")
        for x in range(1,self.n_analyses):
            print "\n\n*********** Round %s ***********\n" %str(x)
            self.parameters["output_file"] = self.PATH + self.output + str(x)
            time_start = time.time()
            analysis = BioSig_console(**self.parameters)
            analysis.run()
            time_end = time.time()
            total_time = time_end - time_start
            output_time.write("Round %s took %s secs\n" %(str(x), str(total_time)))
    
        output_time.close()


if __name__ == "__main__":
    #1
    artificial_analysis_1 = batch("example_artificial_data_anova_001.csv", "example_classes_artificial.txt", 
                                  "output_example_artifical_analysis_1", 3)
    artificial_analysis_1.parameters["chromossome_lenght"] = 388
    artificial_analysis_1.parameters["population_size"] = 40
    artificial_analysis_1.parameters["generations"] = 15000
    artificial_analysis_1.parameters["biosignature_max_size"] = 12
    artificial_analysis_1.parameters["biosignature_min_size"] = 8
    artificial_analysis_1.parameters["biosignature_min_acceptable_size"] = 6
    artificial_analysis_1.parameters["size_weight"] = 30
    artificial_analysis_1.parameters["distance_ratio_weight"] = 10
    artificial_analysis_1.parameters["variation_weight"] = 0
    artificial_analysis_1.parameters["buffer_penalty"] = 50
    artificial_analysis_1.run()

    #2
    artificial_analysis_2 = batch("example_artificial_data_anova_001.csv", "example_classes_artificial.txt", 
                                  "output_example_artifical_analysis_2", 3)
    artificial_analysis_2.parameters["chromossome_lenght"] = 388
    artificial_analysis_2.parameters["population_size"] = 40
    artificial_analysis_2.parameters["generations"] = 15000
    artificial_analysis_2.parameters["biosignature_max_size"] = 12
    artificial_analysis_2.parameters["biosignature_min_size"] = 6
    artificial_analysis_2.parameters["biosignature_min_acceptable_size"] = 5
    artificial_analysis_2.parameters["size_weight"] = 30
    artificial_analysis_2.parameters["distance_ratio_weight"] = 10
    artificial_analysis_2.parameters["variation_weight"] = 0
    artificial_analysis_2.parameters["buffer_penalty"] = 50
    artificial_analysis_2.parameters["n_extremes_values"] = 3    
    artificial_analysis_2.run()
    
