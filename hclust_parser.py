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

class Sample:
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name

class Node:
    def __init__(self,leaf_1, leaf_2, node_id, bootstrap=0):
        """Just a simple class that models a tree node"""
        self.leaves = [leaf_1, leaf_2]
        self.node_id = node_id
        self.bootstrap = bootstrap


    def __getitem__(self):
        return self.leaves

    def __str__(self): 
        return str(self.leaves) + str(node_id) + str(self.bootstrap)


class Tree:
    def __init__(self):
        """Just a simple class that models a Tree"""
        self.root = None
        self.nodes = {}
        self.nodes_alphabetically_ordered = {}
        

    def add_node(self,node, node_id):
        self.nodes[node_id] = node


    def find_node(self, node_id):
        return self.nodes["node_id"]


    def get_newick(self, node):
        """Given a node, return a newick format of it it"""
        nodes = []
        if isinstance(node.leaves[0], Sample):       
            nodes.append(node.leaves[0].name)
        else:
            nodes.append(self.get_labels(node.leaves[0]))


        if isinstance(node.leaves[1], Sample):
            nodes.append(node.leaves[1].name)
        else:
            nodes.append(self.get_labels(node.leaves[1]))


        return tuple(nodes)


    def get_labels(self, node):
        """Given a node, return all names clustered in it"""
        if isinstance(node.leaves[0], Sample) and isinstance(node.leaves[1], Sample):       
            return node.leaves[0].name + "," + node.leaves[1].name

        elif isinstance(node.leaves[0], Sample) and isinstance(node.leaves[1], Node):       
            return node.leaves[0].name + "," + self.get_labels(node.leaves[1])        

        elif isinstance(node.leaves[0], Node) and isinstance(node.leaves[1], Sample):       
            return self.get_labels(node.leaves[0]) + "," + node.leaves[1].names

        elif isinstance(node.leaves[0], Node) and isinstance(node.leaves[1], Node):       
            return self.get_labels(node.leaves[0]) + "," + self.get_labels(node.leaves[1])


    def get_labels_for_each_node(self):
        self.nodes_alphabetically_ordered = {}
        for node_id, node in self.nodes.iteritems():
            node_labels = self.get_labels(node).split(",")
            node_labels.sort()
            self.nodes_alphabetically_ordered[node_id] = node_labels

    def check_node_existance(self, node_to_be_checked, check_bootstrap=True):
        node_to_be_checked.sort()
        if check_bootstrap:
            for node_id, node in self.nodes_alphabetically_ordered.iteritems():
                if node == node_to_be_checked:
                    return True, node_id, self.get_bootstrap_support(node_id)
                else:
                    pass

            return False, 0

        else:
            for node_id, node in self.nodes_alphabetically_ordered.iteritems():
                if node == node_to_be_checked:
                    return True, node_id
                else:
                    pass

            return False



    def get_bootstrap_support(self, node_id):
        return self.nodes[node_id].bootstrap

    def __str__(self):
        return str(self.get_newick(self.root))

      

class pvclust_tree:
    def __init__(self, pvclust_object):
        """This class models/parse an output from pvclust library in R"""
        self.pvclust = pvclust_object
        self.hclust = self.pvclust["hclust"]
        self.hclust_merge = self.hclust["merge"]
        self.tree = Tree()
        self.reconstruct_tree()
        self.tree.get_labels_for_each_node()
              

    def get_label(self,label_index):
        """Given an index to hclust$labels, returnt its string, that should be a sample name"""
        #the -1 is to adapt R's pvclust index range to its correspondent RPy object
        return self.hclust["labels"][label_index-1]


    def reconstruct_tree(self):
        """Given a hclust$merge table, create a Tree object"""
        node_id = 1
        #here somethings to be aware:
        #node or label identification in pvclust object (R) start at number 1
        #however, when they are imported by rpy their index start at number 0
        #node_id is set to R range
        #so if we want to obtain node or label 10, we need to refer to index 9 to its respective list in rpy object
        #
        #in pvclust merge table, if a leaf is a sample, its index is negative
        #but, if it is a internal node, its index is positive
        for node in self.hclust_merge:
           
            if node[0] < 0 and node[1] < 0: #if all leaves are tips/samples
                label_0, label_1 = self.hclust["labels"][int(abs(node[0])-1)], self.hclust["labels"][int(abs(node[1])-1)]
                #print label_0, label_1, node
                node_obj = Node(Sample(label_0), Sample(label_1), node_id, self.get_bootstrap_value(node_id-1))     
                self.tree.add_node(node_obj, node_id)     
                self.tree.root = node_obj

            elif node[0] < 0  and node[1] > 0: #if leaf 0 is a tip and leaf 1 is a internal node
                label_0 = self.hclust["labels"][int(abs(node[0])-1)]
                #print label_0, node
                node_obj = Node(Sample(label_0), self.tree.nodes[abs(node[1])], node_id, self.get_bootstrap_value(node_id-1))
                self.tree.add_node(node_obj, node_id)
                self.tree.root = node_obj

            elif node[0] > 0 and node[1] < 0:
                label_1 = self.hclust["labels"][int(abs(node[1])-1)]
                #print label_1, node
                node_obj = Node(self.tree.nodes[abs(node[0])], Sample(label_1), node_id, self.get_bootstrap_value(node_id-1))
                self.tree.add_node(node_obj, node_id)                 
                self.tree.root = node_obj

            else:
                #print "Only internal nodes", node
                node_obj = Node(self.tree.nodes[abs(node[0])], self.tree.nodes[abs(node[1])], node_id, self.get_bootstrap_value(node_id-1))
                self.tree.add_node(node_obj, node_id)
                self.tree.root = node_obj

            node_id += 1


    def get_bootstrap_value(self, node_id):
        return self.pvclust["edges"]["bp"][node_id]


    def __str__(self):
        return
"""
#from hclust_parser import *
if __name__ == "__main__": 
    from rpy import *
    # load data from microarrays using RPy
    print "Loading Microarray data"
    r('library("pvclust")')
    r('data <- read.table("data_test.txt", header=TRUE)')
    print "Clustering data"
    r('result <- pvclust(data[130:187,2:31], method.dist="euclidian", method.hclust="average", nboot=10)')
    r('plot(result)')
    result_rpy = r("result")
    # parse its result as a python object pvclust_tree
    print "Parsing clustering"
    ct = pvclust_tree(result_rpy)
    a = ct.tree.nodes_alphabetically_ordered[1][:]
    b = ct.tree.nodes_alphabetically_ordered[8][:]
    c = a[:]
    c.append(b[0])
    print "a", a 
    print ct.tree.check_node_existance(a)
    print "b", b
    print ct.tree.check_node_existance(b)
    print "c", c
    print ct.tree.check_node_existance(c)
"""


        


