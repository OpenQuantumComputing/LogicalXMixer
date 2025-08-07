import itertools
import networkx as nx
import numpy as np
import json
from scipy.special import comb
import sys
from PauliString import *
from PauliOperations import *
from Stabilizers import *
from GroupGraph import *
import math
#import openquantumcomputing._rust

from tqdm import tqdm

from sympy import *
from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
#from sympy.physics.quantum import TensorProduct

#from networkx.algorithms.approximation import clique

#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings
#TODO: normalize the Pauli-strings

class Graph:
    def __init__(self, G: nx.Graph, Xl, PS, cost, PS_reduced, cost_reduced, labels=None, positions=None):
        self.G = G
        self.Xl = Xl
        self.PS = PS
        self.PS_reduced = PS_reduced
        self.cost_reduced = cost_reduced
        self.cost = cost
        self.labels = labels
        self.positions = positions
    
    #def to_rust(self) -> openquantumcomputing._rust.Graph:
    #    return openquantumcomputing._rust.Graph(
    #        [list(e) for e in list(self.G.edges)],
    #        self.cost,
    #        self.cost_reduced
    #    )

class Mixer:
    def __init__(self, B, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, sort)
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0
        self.compute_commuting_pairs()

        self.compute_family_of_graphs(blacklist=blacklist)

        if self.reduced:
            if self.digraph:
                self.base_G_reduced=nx.DiGraph()
            else:
                self.base_G_reduced=nx.empty_graph()
            self.base_G_list_reduced=[]
        else:
            if self.digraph:
                self.base_G=nx.DiGraph()
            else:
                self.base_G=nx.empty_graph()
            self.base_G_list=[]

        if whitelist:
            self.select_graphs(whitelist)

    def setB(self, B, sort):
        if isinstance(B, set):
            B = list(B)
        elif isinstance(B, list):
            B = list(set(B))### to make it unique
        else:
            raise TypeError("B must be a list or a set.")
        self.nB = len(B)  ### number of B's
        if sort:
            B=sorted(B, key=lambda x: int(x, 2))

        if len(B) < 2:
            raise Exception("B must contain at least two elements.")

        self.nL = len(B[0])### number of literals
        for b in B:
            if len(b) != self.nL:
                raise Exception("All entries of B must have the same length.")

        self.B = []
        for b in B:
            self.B.append(BitString(1,b))

    @staticmethod
    def HtoString(H):
        ret=[]
        if isinstance(H, Add):
            for sumitem in H.args:### go through all items of the sum (weighted Pauli strings)
                sa=sumitem.args
                f=complex(sa[0])
                if not math.isclose(f.imag,0,abs_tol=1e-7):
                    raise Exception("The imaginary part should be zero", sa)
                ret.append(PauliString(f.real,sa[1]))
        else:### there is only one Pauli string
            sa=H.args
            f=complex(sa[0])
            if not math.isclose(f.imag,0,abs_tol=1e-7):
                raise Exception("The imaginary part should be zero", sa)
            ret.append(PauliString(f.real,sa[1]))
        return ret

    @staticmethod
    def convertLineProjectionToPauliString(state):
        for i in range(len(state.state)):
            if state.state[i]=="0":
                tmp=1/2*(1+Pauli(3))
            else:
                tmp=1/2*(1-Pauli(3))
            if i == 0:
                pauli_str=tmp
            else:
                pauli_str=TensorProduct(pauli_str,tmp)
        return pauli_str

    @staticmethod
    def PauliDecomposition(basisstates):
        # efficient algorithm to express
        # H = \sum_{s in basisstates} |s><s|
        # in Pauli basis
        PS=0
        for i in range(len(basisstates)):
            PS+=Mixer.convertLineProjectionToPauliString(basisstates[i])
        for i in range(len(basisstates[0].state)):
            PS = PS.expand(tensorproduct=True)
        PS=evaluate_pauli_product(PS)
        return Mixer.HtoString(PS)

    def compute_commuting_pairs(self):
        self.commuting_pairs = {}
        for i in range(self.nB):
            for j in range(i + 1, self.nB):
                Xij = Xoperator(self.B[i], self.B[j]).P
                self.commuting_pairs[Xij] = self.commuting_pairs.get(Xij, [])
                self.commuting_pairs[Xij].append([i, j])

    def __create_graph(self, pairs):
        if self.digraph:
            G = nx.DiGraph()
            positions={}
            labels={}
        else:
            G = nx.Graph()
        for b in self.B:
            G.add_node(b.state)
            if self.digraph:
                positions[b.state] = np.array(([int(b.state,2), 0]))
                labels[b.state] = b.state
        for e in pairs:
            G.add_edge(self.B[e[0]].state, self.B[e[1]].state)
            if self.digraph:
                G.add_edge(self.B[e[1]].state, self.B[e[0]].state)
        if self.digraph:
            return G, labels, positions
        else:
            return G

    def process_groups(self, group_graph, X):

    ### 1) minimal set of stabilizer generators
        for i in range(len(group_graph.Xgroup)):
            Xgroup=group_graph.Xgroup[i]
            subV=group_graph.subV[i]
            edges=group_graph.edges[i]
            stabilizer_generators=get_stabilizer_generators(Xgroup, subV[0])

        ### 2) create set of all stabilizers
            if not stabilizer_generators:
                stabilizer_elements=[PauliString(1,"I"*len(subV[0].state))]
            else:
                stabilizer_elements=get_group_elements(stabilizer_generators)
            cost=costPS(stabilizer_elements, PauliString(1,X))
            #print("-----------")
            #print("X=", X)
            #print(Xgroup)
            #print(subV)
            #print(edges)
            #print(stabilizer_generators)
            #print(stabilizer_elements)
            #print("-----------")

        ### 3) perform check
            #subVX=[]
            #for pair in self.commuting_pairs[X]:
            #    for p in pair:
            #        subVX.append(self.B[p])
            try:
                check_projector(stabilizer_elements, subV, self.B)#, subVX=subVX)
            except Exception as inst:
                #TODO: Why do we have to do this?
                continue


        ### 4) optimal projectors restricted to subspace
            stabilizergroup_elements_reduced, cost_reduced = restrict_projector(stabilizer_elements, subV, self.B, PauliString(1,X))
            for sg in stabilizergroup_elements_reduced:
                check_projector(sg, subV, self.B)#, subVX=subVX)

        ### 5) append to family of graphs
            if self.reduced:
                cost_tmp=cost_reduced
                cost_full=self.full_subspace_projector.cost_reduced
            else:
                cost_tmp=cost
                cost_full=self.full_subspace_projector.cost
            if cost_tmp > cost_full:
                if not self.full_subspace_projector_added:
                    self.graph_family.append(self.full_subspace_projector)
                    self.full_subspace_projector_added=True
            else:
                if self.digraph:
                    g, labels, positions =self.__create_graph(edges)
                    self.graph_family.append(Graph(g, X, stabilizer_elements, cost, stabilizergroup_elements_reduced, cost_reduced, labels, positions))
                else:
                    g=self.__create_graph(edges)
                    self.graph_family.append(Graph(g, X, stabilizer_elements, cost, stabilizergroup_elements_reduced, cost_reduced))

    def compute_full_subspace_projector(self, comm_pairs, X):
        subVX=[]
        for pair in comm_pairs:
            for p in pair:
                subVX.append(self.B[p])
        projector=Mixer.PauliDecomposition(subVX)
        check_projector(projector, subVX, self.B)
        cost=costPS(projector, PauliString(1,X))
        projector_reduced, cost_reduced = restrict_projector_nongroup(projector, subVX, self.B, PauliString(1,X))
        for sg in projector_reduced:
            check_projector(sg, subVX, self.B)

        if self.digraph:
            g, labels, positions =self.__create_graph(comm_pairs)
            self.full_subspace_projector = Graph(g, X, projector, cost, projector_reduced, cost_reduced, labels, positions)
        else:
            g=self.__create_graph(comm_pairs)
            self.full_subspace_projector = Graph(g, X, projector, cost, projector_reduced, cost_reduced)
        self.full_subspace_projector_added=False

    def compute_family_of_graphs(self, blacklist=[]):
        print("computing family of graphs")
        self.graph_family=[]
        for X in tqdm(self.commuting_pairs):
            if X in blacklist:
                continue
            comm_pairs=self.commuting_pairs[X]

            self.compute_full_subspace_projector(comm_pairs, X)

            group_graph = GroupGraph(comm_pairs, self.B)
            self.process_groups(group_graph, X)


    def addChain(self, i):

    ### 1) minimal set of stabilizer generators
        z0=self.B[i]
        z1=self.B[i+1]
        EX=[z0, z1]
        X=Xoperator(z0, z1)
        stabilizer_generators=get_stabilizer_generators([X], z0)

    ### 2) create set of all stabilizers
        if not stabilizer_generators:
            stabilizer_elements=[PauliString(1,"I"*len(z0.state))]
        else:
            stabilizer_elements=get_group_elements(stabilizer_generators)
        cost=costPS(stabilizer_elements, X)

    ### 3) perform check
        check_projector(stabilizer_elements, EX, self.B)

    ### 4) optimal projectors restricted to subspace
        stabilizergroup_elements_reduced, cost_reduced = restrict_projector(stabilizer_elements, EX, self.B, X)
        for sg in stabilizergroup_elements_reduced:
            check_projector(sg, EX, self.B)

    ### 5) append to family of graphs
        if self.digraph:
            g, labels, positions =self.__create_graph([[i, i+1]])
            self.graph_chain_family.append(Graph(g, X.P, stabilizer_elements, cost, stabilizergroup_elements_reduced, cost_reduced, labels, positions))
        else:
            g=self.__create_graph([[i, i+1]])
            self.graph_chain_family.append(Graph(g, X.P, stabilizer_elements, cost, stabilizergroup_elements_reduced, cost_reduced))


    def combine_graphs(self, subset):
        if self.reduced:
            G=self.base_G_reduced.copy()
            G_list=self.base_G_list_reduced.copy()
        else:
            G=self.base_G.copy()
            G_list=self.base_G_list.copy()
        for s in subset:
            G = nx.compose(G, s.G)
            G_list.append(s)
        return G, G_list


    def get_chain_mixer(self):

        self.graph_chain_family=[]
        for i in range(self.nB-1):
            X = Xoperator(self.B[i], self.B[i+1]).P
            self.addChain(i)

        G = nx.empty_graph()
        G_list = []
        cost = 0
        for s in self.graph_chain_family:
            G = nx.compose(G, s.G)
            G_list.append(s)
            if self.reduced:
                cost += s.cost_reduced
            else:
                cost += s.cost

        if self.reduced:
            self.solution_chain_reduced = G_list
            self.solution_chain_reduced_cost= cost
        else:
            self.solution_chain = G_list
            self.solution_chain_cost= cost

    def select_graphs(self, Xs):
        graph_family_copy=self.graph_family.copy()

        self.graph_family=[]
        cost={}
        for x in Xs:
            cost[x]=10**10

        subset={}
        for g in graph_family_copy:
            if self.reduced:
                self.base_cost += g.cost_reduced
            else:
                self.base_cost += g.cost
            self.base_nedges += g.G.number_of_edges()
            if g.Xl in Xs:
                if self.reduced:
                    newcost=g.cost_reduced
                else:
                    newcost=g.cost
                if newcost<cost[g.Xl]:
                    subset[g.Xl]=g
                    cost[g.Xl]=newcost
            else:
                self.graph_family.append(g)
        #print(subset)
        #print(list(subset.values()))
        G, G_list = self.combine_graphs(list(subset.values()))

        if self.reduced:
            self.base_G_reduced=G
            self.base_G_list_reduced=G_list
        else:
            self.base_G=G
            self.base_G_list=G_list

    def get_best_mixer_commuting_graphs(self):
        L = len(self.graph_family)
        first = True
        found = False
        for i in range(0, L + 1):
            print(
                i,
                "/",
                L,
                "Number of combinations ",
                L,
                " choose ",
                i,
                "=",
                comb(L, i),
            )
            for subset in tqdm(itertools.combinations(self.graph_family, i)):
                cost = self.base_cost
                for s in subset:
                    if self.reduced:
                        cost += s.cost_reduced
                    else:
                        cost += s.cost
                nedges = self.base_nedges
                for s in subset:
                    nedges += s.G.number_of_edges()
                if self.digraph:
                    nedges = int(nedges / 2)
                # a graph can not be connected if the number of edges is less then the number of nodes-1
                if nedges >= self.nB - 1:
                    if first:
                        G, G_list = self.combine_graphs(subset)
                        if nx.is_connected(G.to_undirected()):
                            if self.reduced:
                                self.solution_reduced=[]
                                self.solution_reduced_Xl=[]
                                self.solution_reduced.append(G_list)
                                self.solution_reduced_cost=cost
                            else:
                                self.solution=[]
                                self.solution_Xl=[]
                                self.solution.append(G_list)
                                self.solution_cost=cost
                            first = False
                            found = True
                    else:
                        if self.reduced:
                            sol_cost=self.solution_reduced_cost
                        else:
                            sol_cost=self.solution_cost
                        if cost <= sol_cost:
                            G, G_list = self.combine_graphs(subset)
                            if nx.is_connected(G.to_undirected()):
                                if cost<sol_cost:
                                    if self.reduced:
                                        self.solution_reduced=[]
                                        self.solution_reduced_Xl=[]
                                    else:
                                        self.solution=[]
                                        self.solution_Xl=[]
                                if self.reduced:
                                    self.solution_reduced.append(G_list)
                                    self.solution_reduced_cost=cost
                                else:
                                    self.solution.append(G_list)
                                    self.solution_cost=cost
                                found = True
                    if cost == 0 and found:
                        break
            if found:
                break


    #def get_best_mixer_rust(self) -> None:
    #    graphs = [g.to_rust() for g in self.graph_family]
    #    problem = openquantumcomputing._rust.Problem(graphs)
    #    solutions = openquantumcomputing._rust.find_minimal_mixer_enumeration(problem, self.reduced)

    #    if self.reduced:
    #        self.solution_reduced: List[Graph] = []
    #        self.solution_reduced_cost = solutions[0].cost
    #    else:
    #        self.solution: List[Graph] = []
    #        self.solution_cost = solutions[0].cost

    #    for solution in solutions:

    #        # TODO: Is there a way to use the same objects throughout?
    #        selected_graphs: List[Graph] = []
    #        for sg in solution.graphs:

    #            for idx, g in enumerate(graphs):
    #                if sg is g:
    #                    selected_graphs.append(self.graph_family[idx])
    #                    break

    #        if self.reduced:
    #            self.solution_reduced.append(selected_graphs)
    #        else:
    #            self.solution.append(selected_graphs)
