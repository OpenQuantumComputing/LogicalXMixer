import itertools
import math
import numpy as np
from sympy import Matrix
from numpy.linalg import matrix_rank
from scipy.linalg import null_space

from BitString import *
from PauliString import *
from PauliOperations import *

from qrs2 import *


def stabilizer(bitstring, i):
    """ Take a bitstring and calculate stabilizer $P\in\pm Z$ of the i-th entry of bitstring, i.e., P*bitstring[i] = bitstring[i]
    @param bitstring BitString
    @return PauliString
    """
    if not isinstance(bitstring, BitString):
        raise Exception("Input must be of type BitString")
    ret = list("I" * len(bitstring.state))
    ret[i] = "Z"
    if bitstring.state[i] == "1":
        return PauliString(-1, "".join(ret))
    else:
        return PauliString(+1, "".join(ret))

def merge_stabilizer_generators(X, S):
    """ given a stabilizer S=<s1,...,s2> with code space C(S) orthogonal code space C(S^X) return C(S^<X>)
    @param X PauliString
    @param S list of PauliStrings
    @return generator of the stabiler of codespace $C(S) + C(S^X)$
    """
    if not isinstance(X, PauliString):
        raise Exception("Input must be of type PauliString")
    if not isinstance(S, list):
        raise Exception("Input must be a list of PauliString")
    ret=[]
    cand = []
    for ind in range(len(S)):
        if X*S[ind]==S[ind]*X:
            ret.append(S[ind])
        else:
            cand.append(S[ind])
    for ind in range(len(cand) - 1):
        ret.append(cand[ind] * cand[ind + 1])
    return ret

def get_stabilizer_generators(Xlist, bitstring):
    """ given a set of bitstrings <Xlist> |bitstring> compute the according code space
    @param Xlist list of PauliStrings
    @param bitstring BitString
    @return generator of the stabiler of codespace $C(S^{<Xlist>})$
    """
    if not isinstance(Xlist, list):
        raise Exception("Input must a list of PauliStrings")
    if not isinstance(bitstring, BitString):
        raise Exception("Input must a BitString")
    S=[]
    for i in range(len(bitstring.state)):
        S.append(stabilizer(bitstring, i))
    for X in Xlist:
        S=merge_stabilizer_generators(X, S)
    return S

def get_group_elements(S):
    """ let S be a list of Pauli strings
    The function returns all elements that can be generated from S
    """
    H = set()
    ### add identity
    Id=PauliString(1,"I"*len(S[0].P))
    H.add(Id)
    for i in range(1, len(S) + 1):
        for subset in itertools.combinations(S, i):
            tmp = Id
            for s in subset:
                tmp = tmp * s
            if not tmp.checkreal():
                raise Exception(
                    "Pauli strings for mixers can not have an imaginary part."
                )
            H.add(tmp)
    return list(H)

def check_V_PS(bitstrings, stabilizers, B, tozero):
    for s in bitstrings:
        res={}
        for ps in stabilizers:
            tmp=ps*s#BitString(1,s)
            res[tmp.state]=res.get(tmp.state,0)
            res[tmp.state]+=tmp.scalar
        for key in res:
            if BitString(1,key) not in B:
                raise Exception("outside B", res, key, B)
            if tozero:
                if not math.isclose(res[key],0,abs_tol=1e-7):
                    raise Exception("must project state",key,"to zero, but is ", res[key], stabilizers)
            else:
                if math.isclose(res[key],0,abs_tol=1e-7):
                    raise Exception("must project state",key,"to non-zero, but is ", res[key], stabilizers)

def check_projector(stabilizers, subV, B, subVX=None):
    other=[]
    if subVX:
        V=subVX
    else:
        V=subV
    for b in B:
        if b not in V:
            other.append(b)
    check_V_PS(other, stabilizers, B, tozero=True)
    check_V_PS(subV, stabilizers, B, tozero=False)

def getNullSpace(A):
    ns= Matrix(A).nullspace()
    ns=np.array(ns).astype(np.float64)
    ret=np.zeros((len(ns),len(ns[0])))
    for i in range(len(ns)):
        ret[i,:]=ns[i].flatten()
    return ret

def getPSs_from_nullspace(PSs,weights):
    if not len(PSs)==len(weights):
        raise Exception("lenghts must be equal")
    ret=[]
    for i in range(len(PSs)):
        tmp=weights[i]*PSs[i]
        if not math.isclose(tmp.scalar,0,abs_tol=1e-7):
            ret.append(tmp)
    return ret
    

def KerA(stabilizers, bitstrings, latex=False):
    n=len(bitstrings)
    m=len(stabilizers)
    A=np.zeros((n,m))
    if latex:
        p_line=" & "
        for j in range(m):
            p_line+=" \\rotatebox[]{90}{$"
            if stabilizers[j].scalar==1:
                p_line+="+"
            else:
                p_line+="-"
            p_line+=stabilizers[j].P+"$}"
            if j==m-1:
                p_line+="\\\\"
            else:
                p_line+=" & "
        print(p_line)
    for i in range(len(bitstrings)):
        if latex:
            p_line="$\\ket{"+str(bitstrings[i].state)+"}$ & "
        for j in range(m):
            A[i,j]=(stabilizers[j]*bitstrings[i]).scalar
            if latex:
                if A[i,j]==1:
                    p_line+="+1"
                else:
                    p_line+="-1"
                if j==m-1:
                    p_line+=" \\\\"
                else:
                    p_line+=" & "
        if latex:
            print(p_line)
    ret=getNullSpace(A)
    return ret

def restrict_projector(stabilizers, bitstrings, B, Xop):
    ret=[]
    bestcnot=10**100

    other=[]
    for b in B:
        if b not in bitstrings:
            other.append(b)
    if not other:
        ret.append([PauliString(1,"I"*len(b.state))])
        bestcnot=Xop.ncnot()
    else:
        # all "other" bitstrings should be projected to zero
        other_null = KerA(stabilizers, other, latex=False)

        for l in range(1,other_null.shape[0]+1):
            for null_subset in itertools.combinations(other_null, l):
                ns=np.zeros_like(null_subset[0])
                for n in null_subset:
                    ns+=n
                # "bitstrings" should not be projected to zero
                if sum(ns)!=0:
                    candidate=getPSs_from_nullspace(stabilizers, ns)
                    ncnot=costPS(candidate, Xop)
                    if ncnot<bestcnot:
                        bestcnot=ncnot
                        ret=[]
                    if ncnot==bestcnot:
                        ret.append(candidate)
    return ret, bestcnot

def isLinearDependent(A, v):
    rankA=matrix_rank(A)
    rankAv=matrix_rank(np.vstack((A,v)))
    if rankAv>rankA:
        return False
    else:
        return True


#def KerUnionB(stabilizers, bitstrings):
#    m=len(stabilizers)
#    for i in range(len(bitstrings)):
#        A=np.zeros((1,m))
#        for j in range(m):
#            A[0,j]=(stabilizers[j]*bitstrings[i]).scalar
#        tmp=getNullSpace(A)
#        if i==0:
#            ret=tmp
#        else:
#            for l in range(tmp.shape[0]):
#                if not isLinearDependent(ret, tmp[l,:]):
#                    ret=np.vstack((ret, tmp[l,:]))
#        if matrix_rank(ret)==m:
#            return ret, True
#    return ret, False


def restrict_projector_nongroup(stabilizers, bitstrings, B, Xop):
    ret=[]

    other=[]
    for b in B:
        if b not in bitstrings:
            other.append(b)

    if not other:
        ret.append([PauliString(1,"I"*len(b.state))])
        bestcnot=Xop.ncnot()
    else:
        bestcnot=10**100

        # all "other" bitstrings should be projected to zero
        other_null = KerA(stabilizers, other, latex=False)

        for l in range(1,other_null.shape[0]+1):
            for null_subset in itertools.combinations(other_null, l):
                ns=np.zeros_like(null_subset[0])
                for n in null_subset:
                    ns+=n
                # "bitstrings" should not be projected to zero
                candidate=getPSs_from_nullspace(stabilizers, ns)
                try:
                    check_V_PS(bitstrings, candidate, B, tozero=False)
                    ncnot=costPS(candidate, Xop)
                    if ncnot<bestcnot:
                        bestcnot=ncnot
                        ret=[]
                    if ncnot==bestcnot:
                        ret.append(candidate)
                except Exception as e:
                    continue
    return ret, bestcnot
