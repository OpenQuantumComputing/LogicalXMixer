#imports

import numpy as np
import itertools

from random import *


# Determine the relative costs of X, CX, and CCX. These may vary depending on the hardware and error mitigation/correction methods
cost_of_x = 1
cost_of_cx = 10
cost_of_ccx = 40



def check_if_redundant(new_gate, gatelist):
    redundant = False
    
    # This part checks if the new gate inverts one of the previous gates:
    for g in gatelist:
        if g == new_gate:
            redundant = True
        elif g[0] == 'cx' and (g[1] in new_gate or g[2] in new_gate):
            redundant = False
        elif g[0] == 'ccx' and (g[1] in new_gate or g[2] in new_gate or g[3] in new_gate):
            redundant = False
    
    # This part makes sure all gates are listed in the right order to avoid multiple equivalent circuits:
    if len(gatelist) > 0:
        g2 = new_gate.copy()[1:]
        last_layer = []
        last_layer_incomplete = True
        i = -1
        while last_layer_incomplete and len(gatelist) + i >= 0:
            g1 = gatelist[i][1:]
            if [value for value in g1 if value in g2] == []:
                i -= 1
                if min(g2) < min(g1):
                    last_layer_incomplete = False
                    redundant = True
            else:
                last_layer_incomplete = False
            
    return redundant

def flip(bitstr):
    new_bitstr = ''
    for b in bitstr:
        if b == '0':
            new_bitstr += '1'
        elif b == '1':
            new_bitstr += '0'
    return new_bitstr


def apply_gate(basis_states, gate):
    new_basis_states = []
    for b in basis_states:
        if gate[0] == 'x':
            q = gate[1]
            new_b = b[:q] + flip(b[q]) + b[q+1:]
        elif gate[0] == 'cx':
            q1, q2 = gate[1], gate[2]
            if b[q1] == '1':
                new_b = b[:q2] + flip(b[q2]) + b[q2+1:]
            else:
                new_b = b
        elif gate[0] == 'ccx':
            q1, q2, q3 = gate[1], gate[2], gate[3]
            if b[q1] == '1' and b[q2] == '1':
                new_b = b[:q3] + flip(b[q3]) + b[q3+1:]
            else:
                new_b = b
        new_basis_states.append(new_b)
    return new_basis_states

def bitwise_mult(bitstr1, bitstr2):
    new_bitstr = ''
    for b in range(len(bitstr1)):
        if bitstr1[b] == '1' and bitstr2[b] == '1':
            new_bitstr += '1'
        else:
            new_bitstr += '0'
    return new_bitstr

def bitsum(bitstr):
    s = 0
    for b in bitstr:
        if b == '1':
            s += 1
    return s
    
def check_if_solution(basis_states, gatelist, num_qubits_to_sort):
    new_basis_states = basis_states.copy()
    overlap = ''
    for i in range(len(basis_states[0])):
        overlap += '1'
    solution = True
    for b in new_basis_states:
        overlap = bitwise_mult(overlap, b)
        if bitsum(overlap) < num_qubits_to_sort:
            solution = False
            return False, new_basis_states, overlap
    return True, new_basis_states, overlap
    
def iterate_through_node(num_qubits, maxdepth, basis_states, num_qubits_to_sort, gatelist, cost, best_solution, best_solution_cost, best_solution_result, best_sorted_bits):
    #print(gatelist, basis_states)
    
    solution, basis_states, sorted_bits = check_if_solution(basis_states, gatelist, num_qubits_to_sort)
    if solution:
        #print('Found a solution! Cost:', cost, gatelist)
        #print('Checking if I can find a better solution...')
        if cost < best_solution_cost:
            best_solution = gatelist
            best_solution_cost = cost
            best_solution_result = basis_states
            best_sorted_bits = sorted_bits
    
    order = ['x','cx','ccx']
    shuffle(order)
    for o in order:
    
        # try all X-gates:
        if o == 'x' and cost + cost_of_x < best_solution_cost and len(gatelist) < maxdepth:
            for q in range(num_qubits):
                new_gate = ['x', q]
                redundant = check_if_redundant(new_gate, gatelist)
                if not redundant:
                    new_basis_states = basis_states.copy()
                    
                    new_gatelist = gatelist.copy()
                    new_gatelist.append(new_gate)
                    
                    #for g in new_gatelist:
                    #    new_basis_states = apply_gate(new_basis_states, g)
                    new_basis_states = apply_gate(new_basis_states, new_gate)
                    
                    new_cost = cost + cost_of_x
                    
                    best_solution, best_solution_cost, best_solution_result, best_sorted_bits = iterate_through_node(num_qubits, maxdepth, new_basis_states, num_qubits_to_sort, new_gatelist, new_cost, best_solution, best_solution_cost, best_solution_result, best_sorted_bits)
                    
        
        # try all CX-gates:
        if o == 'cx' and cost + cost_of_cx < best_solution_cost and len(gatelist) < maxdepth:
            for q1 in range(num_qubits):
                for q2 in range(num_qubits):
                    if q1 != q2:
                        new_gate = ['cx', q1, q2]
                        redundant = check_if_redundant(new_gate, gatelist)
                        if not redundant:
                            new_basis_states = basis_states.copy()
                            
                            new_gatelist = gatelist.copy()
                            new_gatelist.append(new_gate)
                            
                            #for g in new_gatelist:
                            #    new_basis_states = apply_gate(new_basis_states, g)
                            new_basis_states = apply_gate(new_basis_states, new_gate)
                            
                            new_cost = cost + cost_of_cx
                            
                            best_solution, best_solution_cost, best_solution_result, best_sorted_bits = iterate_through_node(num_qubits, maxdepth, new_basis_states, num_qubits_to_sort, new_gatelist, new_cost, best_solution, best_solution_cost, best_solution_result, best_sorted_bits)           
        
        
        # try all CCX-gates:
        if o == 'ccx' and cost + cost_of_ccx < best_solution_cost and len(gatelist) < maxdepth:
            for (q1, q2) in itertools.combinations(range(num_qubits),2):
                for q3 in range(num_qubits):
                    if q3 not in [q1, q2]:
                        new_gate = ['ccx', q1, q2, q3]
                        redundant = check_if_redundant(new_gate, gatelist)
                        if not redundant:
                            new_basis_states = basis_states.copy()
                            
                            new_gatelist = gatelist.copy()
                            new_gatelist.append(new_gate)
                            
                            #for g in new_gatelist:
                            #    new_basis_states = apply_gate(new_basis_states, g)
                            new_basis_states = apply_gate(new_basis_states, new_gate)  
                                                    
                            new_cost = cost + cost_of_ccx
                            
                            best_solution, best_solution_cost, best_solution_result, best_sorted_bits = iterate_through_node(num_qubits, maxdepth, new_basis_states, num_qubits_to_sort, new_gatelist, new_cost, best_solution, best_solution_cost, best_solution_result, best_sorted_bits)
        
                   
            
    return best_solution, best_solution_cost, best_solution_result, best_sorted_bits

def find_solution(basis_states,maxcost=100,maxdepth=5):
    num_qubits = len(basis_states[0])
    k = len(basis_states)
    num_qubits_to_sort = num_qubits - int(np.log2(k))
    
    best_solution = []
    
    #print('Checking for solutions with cost up to', maxcost, 'and depth up to', maxdepth)
    best_solution, best_solution_cost, best_solution_result, best_sorted_bits = iterate_through_node(num_qubits, maxdepth, basis_states, num_qubits_to_sort, [], 0, [], maxcost, [], '')
    
    return best_solution, best_solution_cost, best_solution_result, best_sorted_bits

def find_optimal_permutation(basis_states,selection,maxcost=100,maxdepth=5):
    bs = [basis_states[s] for s in selection]
    best_solution, best_solution_cost, best_solution_result, best_sorted_bits = find_solution(bs,maxcost,maxdepth)
    return best_solution, best_solution_cost

# probably don't need this anymore:
"""
def transform_pauli_string(pauli_string, gatelist):
    pauli_list = []
    for p in pauli_string:
        pauli_list.append(p)
    for g in gatelist:
        if g[0] == 'x':
            pass
        elif g[0] == 'cx':
            if pauli_list[g[1]] == 'I' and pauli_list[g[2]] == 'I':
                pass
            elif pauli_list[g[1]] == 'I' and pauli_list[g[2]] == 'X':
                pass
            elif pauli_list[g[1]] == 'X' and pauli_list[g[2]] == 'I':
                pauli_list[g[2]] = 'X'
            elif pauli_list[g[1]] == 'X' and pauli_list[g[2]] == 'X':
                pauli_list[g[2]] = 'I'
        elif g[0] == 'ccx':      # not sure what to do here
            if pauli_list[g[1]] == 'I' and pauli_list[g[2]] == 'I' and pauli_list[g[3]] == 'I':
                pass
            elif pauli_list[g[1]] == 'I' and pauli_list[g[2]] == 'I' and pauli_list[g[3]] == 'X':
                pass
            elif pauli_list[g[1]] == 'I' and pauli_list[g[2]] == 'X' and pauli_list[g[3]] == 'I':
                pass
            elif pauli_list[g[1]] == 'I' and pauli_list[g[2]] == 'X' and pauli_list[g[3]] == 'X':
                pass
            elif pauli_list[g[1]] == 'X' and pauli_list[g[2]] == 'I' and pauli_list[g[3]] == 'I':
                pass
            elif pauli_list[g[1]] == 'X' and pauli_list[g[2]] == 'I' and pauli_list[g[3]] == 'X':
                pass
            elif pauli_list[g[1]] == 'X' and pauli_list[g[2]] == 'X' and pauli_list[g[3]] == 'I':
                pass
            elif pauli_list[g[1]] == 'X' and pauli_list[g[2]] == 'X' and pauli_list[g[3]] == 'X':
                pass
    new_pauli_string = ''
    for p in pauli_list:
        new_pauli_string += p
    return new_pauli_string

def logical_X(b1, b2):
    gatelist = []
    pauli_string = ''
    n = len(b1)
    for q in range(n):
        if b1[q] == b2[q]:
            pauli_string += 'I'
        else:
            gatelist.append(['X', q])
            pauli_string += 'X'
    return gatelist, pauli_string
    
"""