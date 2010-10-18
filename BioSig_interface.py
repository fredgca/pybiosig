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


import pygtk, gtk.glade
from BioSig import *
import gtk
import gobject
from ma_table_parser import *
gtk.gdk.threads_init()

class Interface:
    ###funcoes da gui####
    def get_glade_widgets(self):
        """Load all widgets from the glade xml file that are listed at 'widget_names'"""
        widget_names = ["main_window", "distance_measure_combo", "clustering_method_combo", "bootstrap_spinbutton",
                        "chromossome_lenght_spinbutton", "generations_spinbutton", "population_size_spinbutton", 
                        "data_matrix_filechooserbutton", "classes_filechooserbutton", "run_button", 
                        "parallel_computing_checkbutton", "filter_ma_data_radiobutton", "distance_measure_4dist_combo",
                        "filter_method_combo", "filter_p_value_combo", "filter_frame", "ga_frame",
                        "clusters_found_weight_spinbutton", "bootstrap_weight_spinbutton", "size_weight_spinbutton", 
                        "distance_weight_spinbutton", "evolve_solution_distance_radiobutton", "distance_parameters_frame",
                        "evolve_solution_cluster_radiobutton", "clustering_frame", "output_entry", 
                        "variation_among_groups_weight_spinbutton", "biosignature_max_size_spinbutton",
                        "biosignature_min_size_spinbutton", "min_acceptable_biosignature_lenght_spinbutton",
                        "buffer_penalty_spinbutton", "ga_selector_combo", "mutation_rate_spinbutton",
                        "crossover_rate_spinbutton", "f_factor_spinbutton", "main_notebook", "n_extreme_values_spinbutton"]

        for widget_name in widget_names:
            setattr(self, "_" + widget_name, self.xml_glade.get_widget(widget_name))

    def on_close_button_clicked(self, *args):
        gtk.main_quit()

    def on_run_button_clicked(self, *args):
        if self._evolve_solution_distance_radiobutton.get_active():
            self.evolve_solution_distance()            
        elif self._evolve_solution_cluster_radiobutton.get_active():
            self.evolve_solution_cluster()
        elif self._filter_ma_data_radiobutton.get_active():
            self.filter_ma_data()


    def evolve_solution_cluster(self):
        parameters = {}
        #Evolve soluion using distance 
        parameters["optimization_method"] = "clusters"
        ## Get all parameter from the interface
        ###Biosignature size parameters
        parameters["biosignature_max_size"] = self._biosignature_max_size_spinbutton.get_value_as_int()
        parameters["biosignature_min_size"] = self._biosignature_min_size_spinbutton.get_value_as_int()
        parameters["biosignature_min_acceptable_size"] = self._min_acceptable_biosignature_lenght_spinbutton.get_value_as_int()
        parameters["buffer_penalty"] = self._buffer_penalty_spinbutton.get_value_as_int()
        ###Weight values
        parameters["groups_weight"] = self._clusters_found_weight_spinbutton.get_value_as_int()
        parameters["bootstrap_weight"] = self._bootstrap_weight_spinbutton.get_value_as_int()
        parameters["size_weight"] = self._size_weight_spinbutton.get_value_as_int()
        parameters["variation_weight"] = self._variation_among_groups_weight_spinbutton.get_value_as_int()#**
        ###Clustering parameter
        parameters["method_dist"] = self._distance_measure_combo.get_active_text()
        parameters["method_hclust"] = self._clustering_method_combo.get_active_text()
        parameters["n_bootstrap"] = self._bootstrap_spinbutton.get_value_as_int()
        ###GA parameters
        parameters["population_size"] = self._population_size_spinbutton.get_value_as_int()
        parameters["generations"] = self._generations_spinbutton.get_value_as_int()
        parameters["chromossome_lenght"] = self._chromossome_lenght_spinbutton.get_value_as_int()
        parameters["ga_selector"] = self._ga_selector_combo.get_active_text()
        parameters["mutation_ratio"] = self._mutation_ratio_spinbutton.get_value()
        parameters["crossover_ratio"] = self._mutation_ratio_spinbutton.get_value()
        parameters["ga_repetition"] = self.ga_repetition_spinbutton.get_value()
        ###I/O data
        parameters["datafile"] = self._data_matrix_filechooserbutton.get_filename()
        parameters["classes_file"] = self._classes_filechooserbutton.get_filename()
        PATH = parameters["datafile"][:parameters["datafile"].rfind("/")+1]
        parameters["output_file"] = PATH + self._output_entry.get_text()
        parallel_computing = self._parallel_computing_checkbutton.get_active()
        ## Send them to the BioSig
        analysis = BioSig(**parameters)
        ## And evolve the solution in a thread
        analysis.start()

    def evolve_solution_distance(self):
        parameters = {}
        #Evolve soluion using distance 
        parameters["optimization_method"] = "distance"
        ## Get all parameter from the interface
        ###Biosignature size parameters
        parameters["biosignature_max_size"] = self._biosignature_max_size_spinbutton.get_value_as_int()
        parameters["biosignature_min_size"] = self._biosignature_min_size_spinbutton.get_value_as_int()
        parameters["biosignature_min_acceptable_size"] = self._min_acceptable_biosignature_lenght_spinbutton.get_value_as_int()
        parameters["buffer_penalty"] = self._buffer_penalty_spinbutton.get_value_as_int()
        ###Weight values
        parameters["size_weight"] = self._size_weight_spinbutton.get_value_as_int()
        parameters["distance_ratio_weight"] = self._distance_weight_spinbutton.get_value_as_int()
        parameters["variation_weight"] = self._variation_among_groups_weight_spinbutton.get_value_as_int()#**
        ###Distance parmeters
        parameters["method_dist"] = self._distance_measure_4dist_combo.get_active_text()
        parameters["F_factor"] = self._f_factor_spinbutton.get_value()
        parameters["n_extremes_values"] = self._n_extreme_values_spinbutton.get_value_as_int()
        ###GA parameters
        parameters["population_size"] = self._population_size_spinbutton.get_value_as_int()
        parameters["generations"] = self._generations_spinbutton.get_value_as_int()
        parameters["chromossome_lenght"] = self._chromossome_lenght_spinbutton.get_value_as_int()
        parameters["ga_selector"] = self._ga_selector_combo.get_active_text()
        parameters["mutation_ratio"] = self._mutation_rate_spinbutton.get_value()
        parameters["crossover_ratio"] = self._crossover_rate_spinbutton.get_value()
        ###I/O data
        parameters["datafile"] = self._data_matrix_filechooserbutton.get_filename()
        parameters["classes_file"] = self._classes_filechooserbutton.get_filename()
        PATH = parameters["datafile"][:parameters["datafile"].rfind("/")+1]
        parameters["output_file"] = PATH + self._output_entry.get_text()
        ## Send them to the BioSig
        analysis = BioSig(**parameters)
        ## And evolve the solution in a thread
        analysis.start()


    def filter_ma_data(self):
        data_matrix_file = self._data_matrix_filechooserbutton.get_filename()
        classes_file = self._classes_filechooserbutton.get_filename()
        PATH = data_matrix_file[:data_matrix_file.rfind("/")+1]
        output_file = PATH + self._output_entry.get_text()
        filter_method = self._filter_method_combo.get_active_text()
        p_value = self._filter_p_value_combo.get_active_text()
        print data_matrix_file, classes_file, output_file, filter_method, p_value
        matrix = MAGE_TAB_data_matrix(data_matrix_file, classes_file, output_file, filter_method, float(p_value))
        matrix.start()


    def on_action_radiobutton_changed(self, *args):
        cluster = self._evolve_solution_cluster_radiobutton.get_active()
        distance = self._evolve_solution_distance_radiobutton.get_active()
      
        if cluster == True and distance == False:
            #Optimizing using clustering criterium
            self._filter_frame.set_sensitive(False)
            self._ga_frame.set_sensitive(True)
            self._clustering_frame.set_sensitive(True)
            self._distance_parameters_frame.set_sensitive(False)
            self._clusters_found_weight_spinbutton.set_sensitive(True)
            self._bootstrap_weight_spinbutton.set_sensitive(True)
            self._size_weight_spinbutton.set_sensitive(True)
            self._distance_weight_spinbutton.set_sensitive(False)
            self._biosignature_max_size_spinbutton.set_sensitive(True)
            self._biosignature_min_size_spinbutton.set_sensitive(True)
            self._min_acceptable_biosignature_lenght_spinbutton.set_sensitive(True)
            self._buffer_penalty_spinbutton.set_sensitive(True)
            self._variation_among_groups_weight_spinbutton.set_sensitive(True)

        elif cluster == False and distance == True:
            #Optimizing using distance criterium
            self._filter_frame.set_sensitive(False)
            self._ga_frame.set_sensitive(True)
            self._clustering_frame.set_sensitive(False)
            self._distance_parameters_frame.set_sensitive(True)
            self._clusters_found_weight_spinbutton.set_sensitive(False)
            self._bootstrap_weight_spinbutton.set_sensitive(False)
            self._size_weight_spinbutton.set_sensitive(True)
            self._distance_weight_spinbutton.set_sensitive(True)
            self._biosignature_max_size_spinbutton.set_sensitive(True)
            self._biosignature_min_size_spinbutton.set_sensitive(True)
            self._min_acceptable_biosignature_lenght_spinbutton.set_sensitive(True)
            self._buffer_penalty_spinbutton.set_sensitive(True)
            self._variation_among_groups_weight_spinbutton.set_sensitive(True)

        elif cluster == False and distance == False:
            #Filtering
            self._filter_frame.set_sensitive(True)
            self._ga_frame.set_sensitive(False)
            self._clustering_frame.set_sensitive(False)
            self._distance_parameters_frame.set_sensitive(False)
            self._clusters_found_weight_spinbutton.set_sensitive(False)
            self._bootstrap_weight_spinbutton.set_sensitive(False)
            self._size_weight_spinbutton.set_sensitive(False)
            self._distance_weight_spinbutton.set_sensitive(False)
            self._biosignature_max_size_spinbutton.set_sensitive(False)
            self._biosignature_min_size_spinbutton.set_sensitive(False)
            self._min_acceptable_biosignature_lenght_spinbutton.set_sensitive(False)
            self._buffer_penalty_spinbutton.set_sensitive(False)
            self._variation_among_groups_weight_spinbutton.set_sensitive(False)


    def __init__ (self):
        """BioSig Interface"""
        self.xml_glade= gtk.glade.XML("gui/interface.glade")
        # --- Dicionario com as funcoes callback ---
        funcoes_callback = {"on_close_button_clicked":self.on_close_button_clicked,
                            "on_run_button_clicked":self.on_run_button_clicked,
                            "on_action_radiobutton_changed" : self.on_action_radiobutton_changed

                           }
        self.get_glade_widgets()
        self._main_window.show()
        self.xml_glade.signal_autoconnect(funcoes_callback)
        self._main_window.connect("delete-event", self.on_close_button_clicked)

        self._distance_measure_combo.set_active(0)
        self._distance_measure_4dist_combo.set_active(0)
        self._clustering_method_combo.set_active(0)
#        self._filter_ma_data_radiobutton.set_group(self._evolve_solution_cluster_radiobutton)
#        self._filter_ma_data_radiobutton.connect("group-changed", self.on_action_radiobutton_changed, "Button 1")
#        self._evolve_solution_cluster_radiobutton.connect("group-changed", self.on_action_radiobutton_changed, "Button 1")
        self._main_notebook.set_tab_label_text(self._main_notebook.get_nth_page(0),"Genetic Algorithms parameters")
        self._main_notebook.set_tab_label_text(self._main_notebook.get_nth_page(1),"Objective Function parameters")
        self._main_notebook.set_tab_label_text(self._main_notebook.get_nth_page(2),"Unspecific Filtering parameters")

        #GA parameter
        self._population_size_spinbutton.set_value(40)
        self._generations_spinbutton.set_range(10,999999)
        self._generations_spinbutton.set_value(5000)
        self._bootstrap_spinbutton.set_value(15)
        self._bootstrap_spinbutton.set_range(0,100000)
        self._filter_method_combo.set_active(1)
        self._filter_p_value_combo.set_active(1)
        self._chromossome_lenght_spinbutton.set_range(10,50000)
        self._chromossome_lenght_spinbutton.set_value(324)
        self._ga_selector_combo.set_active(1)
        self._mutation_rate_spinbutton.set_value(0.02)
        self._mutation_rate_spinbutton.set_increments(0.01, 0.1)
        self._mutation_rate_spinbutton.set_range(0.0,1.0)
        self._crossover_rate_spinbutton.set_value(0.8)
        self._crossover_rate_spinbutton.set_increments(0.05, 0.1)
        self._crossover_rate_spinbutton.set_range(0.0,1.0)

        #weights
        self._clusters_found_weight_spinbutton.set_range(0,100000)
        self._clusters_found_weight_spinbutton.set_value(1)
        self._clusters_found_weight_spinbutton.set_sensitive(False)
        self._bootstrap_weight_spinbutton.set_range(0,100000)
        self._bootstrap_weight_spinbutton.set_value(10)
        self._bootstrap_weight_spinbutton.set_sensitive(False)
        self._size_weight_spinbutton.set_value(10)
        self._distance_weight_spinbutton.set_value(10)
        self._size_weight_spinbutton.set_range(0,100000)
        self._distance_weight_spinbutton.set_range(0,100000)
        self._variation_among_groups_weight_spinbutton.set_range(0,100)
        self._variation_among_groups_weight_spinbutton.set_value(0)
        self._parallel_computing_checkbutton.set_active(False)
        self._parallel_computing_checkbutton.hide()
        self._evolve_solution_distance_radiobutton.set_active(True)
        self._filter_frame.set_sensitive(False)
        self._f_factor_spinbutton.set_range(0,10)
        self._f_factor_spinbutton.set_value(0.2)
        self._f_factor_spinbutton.set_increments(0.1, 0.1)

        #Biosignature size parameters
        self._biosignature_max_size_spinbutton.set_range(5,1000)
        self._biosignature_max_size_spinbutton.set_value(16)
        self._biosignature_min_size_spinbutton.set_range(3,500)
        self._biosignature_min_size_spinbutton.set_value(8)
        self._min_acceptable_biosignature_lenght_spinbutton.set_range(2,450)
        self._min_acceptable_biosignature_lenght_spinbutton.set_value(4)
        self._buffer_penalty_spinbutton.set_range(1,50)
        self._buffer_penalty_spinbutton.set_value(1)

        #Distance optimization parameters
        self._n_extreme_values_spinbutton.set_range(1,100)
        self._n_extreme_values_spinbutton.set_value(3)

if __name__ == "__main__":              
    Interface()
    gtk.gdk.threads_enter()
    gtk.main()
    gtk.gdk.threads_leave()                  
