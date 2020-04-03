import networkx as nx
import arnie.RNAGraph.RG_utils as utils
from copy import copy

class RNAGraph(object):

    def __init__(self, secstruct):
        '''Create a NetworkX graph representing an RNA secondary structure, where nodes represent
        hairpin loops, internal loops, bulges, anything not a helix, and stems represent helices.

        Edge weights / lengths are equal to the number of base pairs present in the helix.
        Nodes are sized according to the number of unpaired bases in the loop.

        Input: secondary structure in dot-parentheses notation. Currently cannot handle pseudoknots.

        Attributes:

        G: NetworkX Directed Graph object. Amenable to much more analysis.

        n_hairpins, n_internal_loops, n_multiloops : number of loops
        n_helices : number of helices
        max_ladder_distance : maximum end_to_end distance of helices
         present in structure, not counting lengths of loops.

        loop_sizes (dictionary) : keys correspond to loop numbering (node numbering),
         values hold the number of unpaired bases present in each loop.

        NB: node numberings don't per se correspond to their order in the tree!
         Will be fixed in parse_stems util eventually.

        '''

        self.secstruct=secstruct
        self.stems = utils.parse_stems_from_bps(utils.convert_structure_to_bps(secstruct))
        self.stem_assignment = utils.get_stem_assignment(secstruct)
        self.pairmap = utils.get_pairmap(secstruct)
        self.G = nx.DiGraph()
        self.loop_sizes = { 0:0 }

        self.setup_graph()

        self.n_helices = len(list(self.G.edges))
        self.MLD = nx.algorithms.dag.dag_longest_path_length(self.G)

        self.n_hairpins, self.n_internal_loops, self.n_3WJs, self.n_4WJs, self.n_5WJs_up = self.count_loops()

        
    def setup_graph(self):
        '''Create graph by reading pairmap array and recursively creating edges.'''
        
        self.G.add_node(0)
        jj = 0
        while (jj < len(self.pairmap)):

            if self.pairmap[jj] > jj:

                # we're in the structure here
                self.add_edges_r_(jj, self.pairmap[jj],0,0)
                jj = self.pairmap[jj]+1

            else:
                #we're in external loop here
                jj += 1
                self.loop_sizes[0] += 1


    def add_edges_r_(self, start_index, end_index, last_helix, last_loop):
        '''Recursive method to add edges to graph.'''

        if (start_index > end_index):
            print('Error, start_index > end_index')
            sys.exit(0)
            
        if(self.pairmap[start_index] == end_index):
            #print("we're helix %d! keep going" % stem_assignment[start_index])
            self.add_edges_r_(start_index+1, end_index-1, self.stem_assignment[start_index], last_loop)
        else:
            jj = start_index
            while jj <=end_index:
                if self.pairmap[jj] > jj:
                    self.add_edges_r_(jj, self.pairmap[jj],self.stem_assignment[jj], last_loop)
                    jj = self.pairmap[jj] + 1 #leaving helix
                else:
                    #print(jj,'in a loop!')
                    #print('Last helix was %d, last loop was %d' % (last_helix, last_loop))
                    if last_helix not in list(self.G.nodes):
                        #print(last_helix, len(stems[int(last_helix-1)]))
                        self.G.add_edge(int(last_loop), int(last_helix),
                                   length = len(self.stems[int(last_helix-1)]),
                                  weight = len(self.stems[int(last_helix-1)]))

                        last_loop = copy(last_helix)

                    if int(last_helix) not in self.loop_sizes.keys():
                        self.loop_sizes[int(last_helix)] = 0

                    self.loop_sizes[int(last_helix)] += 1
                        
                    jj+=1

    def nx_draw(self, node_scaling = 1):
        nx.draw_kamada_kawai(self.G, node_size = [node_scaling*100*x for x in self.loop_sizes.values()])

    def count_loops(self):
        n_1, n_2, n_3, n_4, n_5 = 0,0,0,0,0
        for x in list(self.G.degree):
            if x[1]==1:
                n_1 +=1
            elif x[1] == 2:
                n_2 += 1
            elif x[1] == 3:
                n_3 += 1
            elif x[1] == 4:
                n_4 += 1
            elif x[1] > 4:
                n_5 += 1
        return n_1, n_2, n_3, n_4, n_5

    def get_info(self):
        print("Max ladder distance: %d" % self.MLD)
        print("n_helices: %d" % self.n_helices)
        print("n_hairpins: %d" % self.n_hairpins)
        print("n_3WJs: %d" % self.n_3WJs)
        print("n_4WJs: %d" % self.n_4WJs)
        print("n_5WJs_up: %d" % self.n_5WJs_up)


                
                
