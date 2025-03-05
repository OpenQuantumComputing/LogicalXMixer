from BruteForcePermutation import *
import unittest
import numpy as np
from random import *

# chose these values for better testing:
cost_of_x = 1
cost_of_cx = 10
cost_of_ccx = 40 

ket0 = np.array([[1],[0]])

ket1 = np.array([[0],[1]])

I = np.array([[1,0],[0,1]])

X = np.array([[0,1],[1,0]])

ketbra0 = np.array([[1,0],[0,0]])

ketbra1 = np.array([[0,0],[0,1]])


class Test(unittest.TestCase):  


    def matplus(self,U,V):
        M = U.copy()
        for i in range(len(U)):
            for j in range(len(U[i])):
                M[i][j] += V[i][j]
        return M

    def x(self,q,n):
        M = X
        for i in range(q):
            M = np.kron(I,M)
        for i in range(n-q-1):
            M = np.kron(M,I)
        return M

    def cx(self,ctrl,targ,n):
        M0 = np.array([[1]])
        M1 = np.array([[1]])
        
        for i in range(n):
            if i == ctrl:
                M0 = np.kron(M0,ketbra0)
                M1 = np.kron(M1,ketbra1)
            elif i == targ:
                M0 = np.kron(M0,I)
                M1 = np.kron(M1,X)
            else:
                M0 = np.kron(M0,I)
                M1 = np.kron(M1,I)
        M = self.matplus(M0,M1)
        return M

    def ccx(self,ctrl1,ctrl2,targ,n):
        
        M00 = np.array([[1]])
        M01 = np.array([[1]])
        M10 = np.array([[1]])
        M11 = np.array([[1]])
        
        for i in range(n):
            if i == ctrl1:
                M00 = np.kron(M00,ketbra0)
                M01 = np.kron(M01,ketbra0)
                M10 = np.kron(M10,ketbra1)
                M11 = np.kron(M11,ketbra1)
            elif i == ctrl2:
                M00 = np.kron(M00,ketbra0)
                M01 = np.kron(M01,ketbra1)
                M10 = np.kron(M10,ketbra0)
                M11 = np.kron(M11,ketbra1)
            elif i == targ:
                M00 = np.kron(M00,I)
                M01 = np.kron(M01,I)
                M10 = np.kron(M10,I)
                M11 = np.kron(M11,X)
            else:
                M00 = np.kron(M00,I)
                M01 = np.kron(M01,I)
                M10 = np.kron(M10,I)
                M11 = np.kron(M11,I)
        M = self.matplus(M00,M01)
        M = self.matplus(M,M10)
        M = self.matplus(M,M11)
        return M

    def i(self,n):
        M = I
        for i in range(n-1):
            M = np.kron(M,I)
        return M

    def list_to_unitary(self,gatelist,n):
        U = self.i(n)
        for gate in gatelist:
            if gate[0] == 'x':
                U = np.matmul(self.x(gate[1],n),U)
            elif gate[0] == 'cx':
                U = np.matmul(self.cx(gate[1],gate[2],n),U)
            elif gate[0] == 'ccx':
                U = np.matmul(self.ccx(gate[1],gate[2],gate[3],n),U)
        return U

    def bitstring_to_sv(self,bitstring):
        sv = np.array([[1]])
        for b in bitstring:
            if b == '0':
                sv = np.kron(sv,ket0)
            elif b == '1':
                sv = np.kron(sv,ket1)
        return sv

    def sv_to_bitstring(self,sv):
        #print(sv)
        for i in range(len(sv)):
            #print(i, sv[i,0])
            if sv[i,0] == 1:
                j = i
        n = np.log2(len(sv))
        s = '0' + str(int(n)) + 'b'
        #print(j,s)
        return format(j, s)

    def generate_random_list_of_bitstrings(self,k,n):
        l = sample(list(np.arange(2**n)),2**k)
        l2 = []
        s = '0' + str(int(n)) + 'b'
        for i in l:
            l2.append(format(i,s))
        return l2

    def test_all(self):
        for n in np.arange(2,4):
            for k in range(n):
                bitstrings = self.generate_random_list_of_bitstrings(k,n)
                #print(bitstrings)
                best_solution, best_solution_cost, best_solution_result, best_sorted_bits = find_solution(bitstrings,maxcost=20000,maxdepth=6)
                #print(best_solution, best_solution_result, best_sorted_bits)
                
                #test:
                U = self.list_to_unitary(best_solution,n)
                #print(U)
                result_list = [self.sv_to_bitstring(np.matmul(U,self.bitstring_to_sv(b))) for b in bitstrings]
                #print(result_list)
                self.assertEqual(result_list,best_solution_result)
                
                testlist = [1 for i in range(n)]
                for b in best_solution_result:
                    for i in range(n):
                        testlist[i] *= int(b[i])
                        
                testnumber = 1
                for i in range(n):
                    if best_sorted_bits[i] == '1':
                        testnumber *= testlist[i]
                
                self.assertEqual(testnumber,1)
        
        
if __name__ == '__main__':
    unittest.main()   

