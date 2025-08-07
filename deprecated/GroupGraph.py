import networkx as nx
from PauliOperations import *
import numpy as np

class GroupGraph():

    def __init__(self, indices_pairs, B):
        self.indices_pairs = indices_pairs
        self.B=B
        self.n=len(self.indices_pairs)

        pair0=self.indices_pairs[0]
        self.X=Xoperator(self.B[pair0[0]], self.B[pair0[1]])

        self.Xgroup=[]
        self.edges=[]
        self.subV=[]

        if self.n==1:
            Xgroup=[PauliString(1,"I"*len(self.B[0].state))]
            Xgroup.append(self.X)

            V=[]
            pair=self.indices_pairs[0]
            V.append(self.B[pair[0]])
            V.append(self.B[pair[1]])

            self.Xgroup.append(Xgroup)
            self.edges.append([pair])
            self.subV.append(V)
        else:

            self.E_list={}
            self.XE_list={}
            for i in range(self.n):
                # create a list of X-operators from index i to all other index pairs
                self.compute_E_list(i)
            #print(self.E_list)
            #print(self.XE_list)
            self.Xgroup=[]
            self.subV=[]
            self.compute_maxgroups()

    def compute_E_list(self, from_index):
        # creates a list of X-operators from index from_index to all other index pairs
        pair1=self.indices_pairs[from_index]
        self.E_list[from_index]=[]
        self.XE_list[from_index]=[]
        for i in range(self.n):
            if i==from_index:
                continue
            pair2=self.indices_pairs[i]
            E=Xoperator(self.B[pair1[0]], self.B[pair2[0]])
            self.E_list[from_index].append(E)
            self.XE_list[from_index].append(E*self.X)

    def compute_adjacencymatrix(self, from_index):
        A=np.zeros((self.n-1,self.n-1))
        for i in range(self.n-1):
            for j in range(i+1,self.n-1):
                Ei_Ej=self.E_list[from_index][i]*self.E_list[from_index][j]
                if (Ei_Ej in self.E_list[from_index]) or (Ei_Ej in self.XE_list[from_index]):
                    A[i,j]=1
                    A[j,i]=1
        #print(A)
        return A

    def compute_maxgroups(self):
        for from_index in range(self.n):
            A = self.compute_adjacencymatrix(from_index)
            self.G=nx.Graph(A)
            for cl in list(nx.find_cliques(self.G)):
                #print("clique=", cl)
                Xgroup=[PauliString(1,"I"*len(self.B[0].state))]
                Xgroup.append(self.X)
                V=[]
                pair=self.indices_pairs[from_index]
                V.append(self.B[pair[0]])
                V.append(self.B[pair[1]])
                edges=[]
                edges.append(pair)
                for i in cl:
                    ind=i
                    if i >= from_index:
                        ind+=1
                    Xgroup.append(self.E_list[from_index][i])
                    Xgroup.append(self.XE_list[from_index][i])
                    pair=self.indices_pairs[ind]
                    V.append(self.B[pair[0]])
                    V.append(self.B[pair[1]])
                    edges.append(pair)
                Xgroup = [a for a in sorted(Xgroup, key=lambda p: p.P)]
                if Xgroup not in self.Xgroup:
                    self.Xgroup.append(Xgroup)
                    self.subV.append(V)
                    self.edges.append(edges)


    #def compute_Xs(self,i,j):
    #    ret=set()
    #    for a in self.EX[i]:
    #        for b in self.EX[j]:
    #            Xop=Xoperator(a,b)
    #            if Xop*a in self.VX:
    #                ret.add(Xop)
    #    return ret
    #
    #def compute_vertices(self):
    #    self.vertices={}
    #    self.vertices_set={}
    #    n=len(self.indices_pairs)
    #    for i in range(n):
    #        self.vertices[i]={}
    #        self.vertices_set[i]=[]
    #        for j in range(n):
    #            if i!=j:
    #                tmp=self.compute_Xs(i,j)
    #                self.vertices[i][j]=tmp
    #                for x in tmp:
    #                    self.vertices_set[i].append(x)
    #
    #def createGraphs(self):
    #    self.graphs={}
    #    if self.compute_labels:
    #        self.graphs_labels={}
    #    n=len(self.indices_pairs)
    #    for i in range(n):
    #        self.graphs[i]=nx.Graph()
    #        if self.compute_labels:
    #            self.graphs_labels[i]={}
    #        for vn, Xls in self.vertices[i].items():
    #            self.graphs[i].add_node(vn)#num)
    #            if self.compute_labels:
    #                self.graphs_labels[i][vn]=r"$\{"+list(Xls)[0].P+","+list(Xls)[1].P+"\}e_"+str(vn+1)+"$"
    #        for vn_a, Xls_a in self.vertices[i].items():
    #            for vn_b, Xls_b in self.vertices[i].items():
    #                if list(Xls_a)[0]*list(Xls_b)[0] in self.vertices_set[i]:
    #                    self.graphs[i].add_edge(vn_a, vn_b)

