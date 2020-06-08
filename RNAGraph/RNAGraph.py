import networkx as nx
import arnie.RNAGraph.RG_utils as utils
from copy import copy
from networkx.drawing.nx_agraph import graphviz_layout

import numpy as np
class RNAGraph(object):

    def __init__(self, secstruct, colors=None):
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

        self.secstruct = secstruct
        self.stems = utils.parse_stems_from_bps(utils.convert_structure_to_bps(secstruct))
        self.stem_assignment = utils.get_stem_assignment(secstruct)
        self.pairmap = utils.get_pairmap(secstruct)
        self.G = nx.DiGraph()
        self.loop_sizes = { 0:0 }
        self.N = len(self.secstruct)

        if colors is not None:
            #color_dct = ['#3288bd','k','#5e4fa2'] # from Brewer Spectral10
            #color_dct = ['#003EF5','#ffff33']
            color_dct =  ['#000000', '#778899']
            #linewidth_dct = [2.5,1.0]
            self.colors=[color_dct[i] for i in colors]
            #self.linewidths = [linewidth_dct[i] for i in colors]
        else:
            self.colors=['black']*len(self.secstruct)
            #self.linewidths = [1.0]*len(self.secstruct)

        self.setup_graph()

        nodes = [n for n in list(self.G.nodes) if not isinstance(n,str)]
        subgraph = self.G.subgraph(nodes)

        self.MLD = nx.algorithms.dag.dag_longest_path_length(subgraph)

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

        stem_assignment_left = np.concatenate([np.array([-1]), self.stem_assignment[:-1]])
        stem_assignment_right = np.concatenate([self.stem_assignment[1:], np.array([-1])])

        for i,pair in enumerate(self.pairmap):
            if pair==-1:
                self.G.add_node('n%d' % i)

                if stem_assignment_left[i] > 0:
                    if self.secstruct[i-1]=='(':
                        letter='a'
                    elif self.secstruct[i-1]==')':
                        letter='b'
                    self.G.add_edge('n%d' % i,'h%d%s' % (self.stem_assignment[i-1], letter), len=1, weight=1, color=self.colors[i])
                    # add helix_a node here

                if stem_assignment_right[i] > 0:

                    if self.secstruct[i+1]==')':
                        letter='a'
                    elif self.secstruct[i+1]=='(':
                        letter='b'
                    self.G.add_edge('n%d' % i,'h%d%s' % (self.stem_assignment[i+1], letter), len=1, weight=1, color=self.colors[i])
                    # add helix_b node here

        nuc_nodes = [n for n in list(self.G.nodes) if isinstance(n,str) and n.startswith('n')] # hacky
        for nuc in nuc_nodes:
            ind=int(nuc.replace('n',''))
            if 'n%d' % (ind-1) in nuc_nodes:
                self.G.add_edge('n%d' % (ind-1), 'n%d' % ind, len=1, weight=1, color=self.colors[ind-1])
                #print('n%d' % (ind-1), 'n%d' % ind,)


        for stem_ind in range(1,len(self.stems)+1):
            stem_length = len(self.stems[stem_ind-1])
            color_ind = np.where(self.stem_assignment==stem_ind)[0][0]
            self.G.add_edge('h%da' % (stem_ind),'h%db' % (stem_ind), len=stem_length, weight=stem_length,\
             color=self.colors[color_ind])

        # fix for stems that have no unpaired bases between them
        # dbs = [i for i in range(len(self.secstruct)) if self.secstruct[i:i+2] == ')(']
        # for i in dbs:
        #     stem_ind_1 = self.stem_assignment[i]
        #     stem_ind_2 = self.stem_assignment[i+1]
        #     self.G.add_edge('h%db' % stem_ind_1,'h%db' % stem_ind_2, length=1, weight=1,color='black')

        for i in range(self.N-1):

            if self.stem_assignment[i+1] != self.stem_assignment[i]:

                if self.stem_assignment[i+1] != 0.0 and self.stem_assignment[i] != 0.0:
                    stem_ind_1 = self.stem_assignment[i]
                    stem_ind_2 = self.stem_assignment[i+1]
                    color_ind = np.where(self.stem_assignment==stem_ind_1)[0][0]

                    # print(stem_ind_1, stem_ind_2)
                    # print(self.secstruct[i:i+2])

                    if self.secstruct[i:i+2] == '((':
                        self.G.add_edge('h%da' % stem_ind_1,'h%db' % stem_ind_2, len=1, weight=1,color=self.colors[color_ind])
                    elif self.secstruct[i:i+2] == ')(':
                        self.G.add_edge('h%db' % stem_ind_1,'h%db' % stem_ind_2, len=1, weight=1,color=self.colors[color_ind])
                    elif self.secstruct[i:i+2] == '))':
                        self.G.add_edge('h%db' % stem_ind_1,'h%da' % stem_ind_2, length=1, weight=1,color=self.colors[color_ind])


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
                                   len = len(self.stems[int(last_helix-1)]),
                                  weight = len(self.stems[int(last_helix-1)]))

                        last_loop = copy(last_helix)

                    if int(last_helix) not in self.loop_sizes.keys():
                        self.loop_sizes[int(last_helix)] = 0
                        #self.G.add_edge('n%d'%jj, int(last_helix))

                    self.loop_sizes[int(last_helix)] += 1

                    jj+=1

    def nx_draw(self, node_scaling = 1, loops=True, tree=False, colors=None):

        if loops:
        #node_color=['#ffa500']+['#7AC5CD']*(len(self.loop_sizes.keys())-1)
            plot_nodes = [n for n in list(self.G.nodes) if isinstance(n,str)]
        else:
            plot_nodes = [n for n in list(self.G.nodes) if not isinstance(n,str)]


        subgraph = self.G.subgraph(plot_nodes)
        if loops:
            subgraph = subgraph.to_undirected()
        if tree:
            prog = 'dot'
        else:
            prog='neato'
        pos =graphviz_layout(subgraph,prog=prog)

        colors = [subgraph[u][v]['color'] for u,v in subgraph.edges()]

        nx.draw(subgraph, pos, node_size=0, edge_color=colors, arrows=False)

    def count_loops(self):
        n_1, n_2, n_3, n_4, n_5 = 0,0,0,0,0
        nodes = [n for n in list(self.G.nodes) if not isinstance(n,str)]
        subgraph = self.G.subgraph(nodes)

        for x in list(subgraph.degree):
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
        print("n_hairpins: %d" % self.n_hairpins)
        print("n_internal_loops: %d" % self.n_internal_loops)
        print("n_3WJs: %d" % self.n_3WJs)
        print("n_4WJs: %d" % self.n_4WJs)
        print("n_5WJs_up: %d" % self.n_5WJs_up)


                
                
