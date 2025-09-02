import math
from itertools import combinations
from functools import reduce
import operator

is_power_of_two = lambda x: (x > 0) and (x & (x - 1)) == 0

def ncnot(P) :
    """
    Calculate the number of CNOT gates required to implement a Pauli string.

    Args:
        P (int): Pauli string represented as a binary integer, where each 1 represents a qubit that is acted upon by a Pauli operator (X or Z).

    Returns:
        int: Number of CNOT gates required to implement the Pauli string.
    """
    ncnot = P.bit_count()
    return (ncnot > 1)*2*(ncnot - 1)

def pauli_int_to_str(P, nL, operator="X"):
    """
    Converts a Pauli string represented as an integer to its string representation.

    Args:
        P (int): Pauli string.
        nL (int): Number of qubits.
        operator (str, optional): Whether the 1s represent "X" or "Z" operators. Defaults to "X".

    Raises:
        ValueError: If the operator is not "X" or "Z".

    Returns:
        str: Formatted Pauli string.
    """
    P = f"{P:0{nL}b}"
    P = P.replace("0", "I")
    if operator == "X":
        P = P.replace("1", "X")
    elif operator == "Z":
        P = P.replace("1", "Z")
    else:
        raise ValueError("Operator must be 'X' or 'Z'.")
    return P

def parity(n):
    #using Brian Kernighan's algorithm to check parity (commutation/anti-commutation), parity = 1 if even (commutes), and parity = -1 if odd (anti-commutes)
    parity = 0
    while n:
        parity ^= 1
        n &= n - 1 
    if parity == 0:
        parity = 1
    else:
        parity = -1
    return parity

def convert_to_binary_string(int_values, n):
    """
    Convert integers to binary strings with n bits.
    
    - If int_values is a single int: return binary string
    - If int_values is a list: return list of binary strings or (bin, bin) tuples
    - If int_values is a dict: return dict with binary string keys and values
    """
    
    if isinstance(int_values, int):
        return format(int_values, f'0{n}b')

    elif isinstance(int_values, list):
        result = []
        for i in int_values:
            if isinstance(i, int):
                result.append(format(i, f'0{n}b'))
            elif isinstance(i, tuple) and len(i) == 2:
                if all(isinstance(x, int) for x in i):
                    result.append((format(i[0], f'0{n}b'), format(i[1], f'0{n}b')))
        return result

    elif isinstance(int_values, dict):
        result = {}
        for key, val in int_values.items():
            if isinstance(key, int) and isinstance(val, int):
                result[format(key, f'0{n}b')] = format(val, f'0{n}b')
        return result

    else:
        raise TypeError("Input must be int, list, or dict of ints")

# def is_connected(orbits):
#     """
#     Check if the list of orbits forms a connected structure.
#     Connectivity means you can move from any orbit to any other via overlapping nodes.  

#     Args:
#         orbits (List[Tuple[int, ...]]): Node indices representing the group generated sets.

#     Returns:
#         bool: True if the orbits are connected, False otherwise.
#     """
#     if len(orbits) <= 1:
#         return True

#     # Convert orbits to sets for quick intersection
#     orbit_sets = [set(o) for o in orbits]

#     # Build adjacency list: two orbits connected if their sets overlap
#     adjacency = {i: set() for i in range(len(orbit_sets))}
#     for i in range(len(orbit_sets)):
#         for j in range(i + 1, len(orbit_sets)):
#             if orbit_sets[i].intersection(orbit_sets[j]):
#                 adjacency[i].add(j)
#                 adjacency[j].add(i)

#     visited = set()
#     stack = [0]  # start from first orbit

#     while stack:
#         node = stack.pop()
#         if node not in visited:
#             visited.add(node)
#             stack.extend(adjacency[node] - visited)

#     return len(visited) == len(orbit_sets)

def is_connected(orbits):
    """
    Check if the list of orbits forms a connected structure.
    Connectivity means you can move from any orbit to any other via overlapping nodes.
    
    Args:
        orbits (List[Tuple[int, ...]]): Node indices representing the group generated sets.

    Returns:
        bool: True if the orbits are connected, False otherwise.
    """
    n = len(orbits)
    if n <= 1:
        return True

    orbit_sets = [set(o) for o in orbits]  # Precompute sets for O(1) lookup.
    visited = [False] * n
    stack = [0]
    visited[0] = True
    seen_count = 1

    while stack:
        i = stack.pop()
        si = orbit_sets[i]
        for j in range(n):
            if not visited[j] and si & orbit_sets[j]:  # Direct intersection check.
                visited[j] = True
                seen_count += 1
                if seen_count == n:  # Early exit when fully connected.
                    return True
                stack.append(j)

    return False


def find_best_cost(Xs, Zs_operators):
    """
    Finds the best combination of logical X operators that generates an orbit, minimizing the cost function.
    
    Args:
        Xs (List[int]): A list of logical X operators (int representations) that generates an orbit.
        Zs_operators (List[Tuple[int, str]]): A list of tuples representing the Z-operators, where the first element is the sign(+1/-1) and the second is the int representation of the oeprator.
    
    Return:
        List[int]: A list of the best logical X operators (int representations).
    """
    all_x_operators = []
    n = len(Xs)
    Zs = [string[1] for string in Zs_operators] # Only extract the Z-strings without the sign, as the sign is not used in the cost function
    
    # Generate all combinations of Xs (n-1 Xs for an orbit that has 2^n states) and their corresponding hats
    for r in range(1, n + 1):  # Start from 1 to include single elements
        for combo in combinations(Xs, r):
            # Use XOR operation between all Xs in the combination to find a new X that could be used to generate an orbit
            hat = reduce(operator.xor, combo)
            
            # Append a list that contains first the combination of Xs used, and then the XOR of those Xs (hat)
            all_x_operators.append([combo, hat])
    
    all_costs = {}

    # for each inner list in all_x_operators we have the used Xs and the combination of these Xs. 
    # We take each combination (X_combos) and calculate the cost of the ncnot for each Z in Zs to get the total cost of using that combination of Xs with the projector (Zs)
    for used_Xs, X_combos in all_x_operators:
        total_cost = 0
        for Z in Zs:
            cost = ncnot(X_combos | Z) # We use the bitwise OR to find which qubits are acted upon by the Xs and Zs together and pass it to calculate the cost
            total_cost += cost #add up the cost of all Zs for that combination of Xs
        
        #NB I CHANGED COMBOS TO THEIR XOR
        all_costs[X_combos] = total_cost #Here we store which Xs were usied and the total cost of using them with the Zs
    
    # Brute force way
    # Sorting the keys from lowest to highest cost
    brute_all_keys = sorted(all_costs.keys(), key=lambda k: all_costs[k])
    brute_cost = float('inf')
    brute_Xs = []

    # Finding all combinations of n Xs (to generate an orbit of size 2^n) and checking which combination has the lowest cost and is linearly independent
    for combo in combinations(brute_all_keys, n):
        total_cost = sum(all_costs[key] for key in combo)
        
        # Disregard the combination if its cost is already higher than the best found cost
        if total_cost >= brute_cost:
            continue
        
        # Check if the combination is linearly independent and therefore valid
        if not is_independent(combo):
            continue
        
        else: 
            brute_cost = total_cost
            brute_Xs = combo
            
    """
    if n == 1:
        best_Xs_reduced = all_x_operators[0][1]  
        best_cost = all_costs[all_x_operators[0][1]]
        return best_Xs_reduced, best_cost

    best_Xs_reduced = []  
    best_cost = float('inf')  
    
    #the sorted keys, from highest to lowest cost
    all_keys = sorted(all_costs.keys(), key=lambda k: all_costs[k], reverse=True)
    processed_keys = set()
    semi_processed_keys = set()

    for i in range(len(all_keys)):
        combos = []
        key_copy = [key for key in all_keys if key not in processed_keys]
        combos.append(([], key_copy, [], 0))
        print("We go into this loop with all_keys:", key_copy, "\nThis is iteration:", i)
        
        while combos:
            #pops the Xs with the lowest cost first
            best_Xs, available_keys, unavailable_combos, total_cost = combos.pop()
            print("this was just popped from combos:", best_Xs, available_keys, unavailable_combos, total_cost)

            if len(available_keys) < n:
                print("processed_keys:", processed_keys, "semi_processed_keys:", semi_processed_keys, "available_keys:", available_keys)
                semi_processed_keys = set()
                available_keys = [key for key in key_copy if key not in processed_keys]
                unavailable_combos = [key for key in processed_keys]
                print("Available keys after filtering:", available_keys)
                combos.append(([], available_keys, unavailable_combos, 0))
                break

            if len(best_Xs) == 0:
                print("This loop is triggered, processed_keys:", processed_keys, "semi_processed_keys:", semi_processed_keys)
                popped_key = available_keys.pop()
                processed_keys.add(popped_key)
                best_Xs = [popped_key]
                unavailable_combos = [popped_key]
                total_cost += all_costs[popped_key]
                print("the new combo is: best_Xs:", best_Xs, "available_keys:", available_keys, "unavailable_combos:", unavailable_combos, "total_cost:", total_cost)
                combos.append((best_Xs, available_keys, unavailable_combos, total_cost))
                continue
            
            
            # best_key = available_keys.pop()
            # best_Xs.append(best_key)
            # processed_keys.add(best_key) 
            # unavailable_keys.append(best_key) 
            # total_cost = 0 
            # total_cost += all_costs[best_key] 

            # print("Best key:", best_key, "available_keys:", available_keys, "best_Xs:", best_Xs, "unavailable_keys:", unavailable_keys)
            
            for key in reversed(available_keys):  
                print("we are now in this loop")
                # if set(key) in unavailable_keys or set(best_Xs) in unavailable_keys:
                #     print("THIS IF-STATEMENT TRIGGERED")
                #     continue
                if key in (processed_keys or semi_processed_keys):
                    print("This is triggered, key is: ", key, "and processed_keys:", processed_keys, "semi_processed_keys:", semi_processed_keys)
                    continue
                
                best_Xs.append(key)
                
                if set(best_Xs) in unavailable_combos:
                    print("THIS IF-STATEMENT TRIGGERED")
                    best_Xs.pop()  # Remove the last added key if it is already in unavailable_combos
                    continue
                
                # processed_keys.add(key)
                print("Best Xs after appending best_key and key:", best_Xs)
                # adding_key = set(key)
                # new_unavailable_combos = [tuple(sorted(set().union(*combo, adding_key))) for r in range(1, len(unavailable_keys) + 1) for combo in combinations(unavailable_keys, r)]
                # print("New unavailable combos:", new_unavailable_combos)
                # unavailable_combos.extend(new_unavailable_combos)
                # unavailable_keys.append(key)

                total_cost += all_costs[key]
                new_unavailable_combos = [key ^ x for x in unavailable_combos if (key^x) not in unavailable_combos]  # XOR the new key with all unavailable combos
                unavailable_combos.extend(new_unavailable_combos + ([key] if (key not in unavailable_combos) else []))

                if total_cost > best_cost:
                    print("This resulted in a cost that is higher than the best cost, so we break", "total_cost:", total_cost, "best_cost:", best_cost)
                    best_Xs.pop()
                    total_cost -= all_costs[key]
                    combos.append((best_Xs, available_keys, unavailable_combos, total_cost))
                    continue
                
                if len(best_Xs) == n:  
                    best_Xs_reduced = best_Xs.copy()  # Create a copy instead of sharing reference
                    best_cost = total_cost
                    # unavailable_combos.append(best_Xs)
                    print("Best Xs reduced:", best_Xs_reduced, "Total cost:", best_cost)
                    # combos.append((best_Xs, available_keys, unavailable_keys, total_cost))
                    # best_Xs = []
                    # available_keys = [key for key in key_copy if key not in processed_keys]
                    # total_cost = 0
                    # break
                    best_Xs.pop()
                    total_cost -= all_costs[key]
                    combos.append((best_Xs, available_keys, unavailable_combos, total_cost))
            
            backtracking = best_Xs.pop()
            total_cost -= all_costs[backtracking]
            semi_processed_keys.add(backtracking)
            combos.append((best_Xs, available_keys, unavailable_combos, total_cost))
        
        

    # while len(best_Xs) < n:
    #     # We start by selecting the lowest cost from all_costs as we want to minimize the cost
    #     lowest_cost = min(all_costs.values())
    #     keys = [k for k, v in all_costs.items() if v == lowest_cost]
        
    #     # iterate through the keys with the lowest cost
    #     for key in keys:
    #         # Checks that either the key adds to the subset or that it is already covered (i.e. that we are actually creating an orbit)
    #         if (not set(key).issubset(covered)) or (required == covered):  # If key has *any* uncovered elements
    #             # Checks that if it is already covered, we use the lowest cost from maybe_later
    #             if required == covered:
    #                 new_key_and_cost = maybe_later.pop(0) if maybe_later else [key, lowest_cost]
    #                 best_Xs.append(new_key_and_cost[0])
    #                 best_cost += new_key_and_cost[1]
                
    #             # if the required set is not covered, we add the key to the best_Xs and update the covered set
    #             else:
    #                 covered.update(key)
    #                 best_Xs.append(key)
    #                 best_cost += lowest_cost
                
    #             # If we have enough Xs to cover the orbit, we break the loop
    #             if len(best_Xs) == n:
    #                 break
    #             # We delete the key from all_costs as we have used it
    #             del all_costs[key]
            
    #         # If the key does not add to the subset, we store it in maybe_later for later use (if we get a covered set and need to add more Xs)
    #         else:
    #             # We delete the key from all_costs as it does not add to the subset
    #             del all_costs[key]  
    #             maybe_later.append([key, lowest_cost])
                
    # # The best_Xs are reduced by applying XOR to the tuples to get the string for the combination of the original Xs
    # best_Xs_reduced = [reduce(operator.xor, x) for x in best_Xs]
    """
    return brute_Xs, brute_cost #best_Xs_reduced, best_cost

def is_independent(combo):
    """Takes in a list of x-operators and checks if they are linearly independent"""
    n = len(combo)
    # Check all non-empty subsets of the combo to see that they are linearly independent
    for r in range(1, n + 1):
        for subset in combinations(combo, r):
            hat = reduce(operator.xor, subset)
            if hat == 0:
                # Found a dependent subset
                return False
    return True



# #checks all combinations of the n Xs to see if any of them are not linearly independent
        # for r in range(1, n + 1):
        #     for hats in combinations(combo, r):
        #         hat = reduce(operator.xor, hats)
        #         print("The hat for this combo is:", hat)
        #         if hat == 0:
        #             continue
    
if __name__ == '__main__':
    # results = find_best_cost([0b0010, 0b0110, 0b1000], [(1, 0b0010), (1, 0b0110), (1, 0b1000), (1, 0b1010), (1, 0b1100), (1, 0b1110)])
    results = find_best_cost([14, 12, 11], [(1, 0), (1, 13)])
    # results = find_best_cost([9, 21], [(1, 0), (1, 25), (-1, 20), (-1, 13), (1, 2), (1, 27), (-1, 22), (-1, 15)])
    # results = find_best_cost([0b1101, 0b1010], [(1, 0), (1, 11), (-1, 14), (-1, 5)])
    # results = find_best_cost([0b0010], [(1, 0), (1, 1), (-1, 4), (-1, 5), (-1, 8), (-1, 9), (1, 12), (1, 13)])
    print("Best combo of Xs (heuristic):", results[0],"\nBest cost (heuristic):", results[1])#, "\nBest combo of Xs (exact):", results[2], "\nBest cost (exact):", results[3])
    #print("Best combco of Xs reduced (heuristic):", results[1])
    
    # family_of_valid_graphs_1 = {0b0010 : [(0, 1), (2, 7), (3, 9), (4, 13), (5, 10), (6, 8), (11, 14)],
    # 0b0111 : [(0, 2), (1, 7), (3, 4), (5, 14), (6, 12), (9, 13), (10, 11)],
    # 0b1010 : [(0, 3), (1, 9), (2, 4), (6, 11), (7, 13), (8, 14), (10, 12)],
    # 0b1101 : [(0, 4), (1, 13), (2, 3), (5, 8), (6, 10), (7, 9), (11, 12)],
    # 0b1110 : [(0, 5), (1, 10), (2, 14), (4, 8), (6, 13), (7, 11), (9, 12)],
    # 0b0001 : [(0, 6), (1, 8), (2, 12), (3, 11), (4, 10), (5, 13), (9, 14)],
    # 0b0101 : [(0, 7), (1, 2), (3, 13), (4, 9), (5, 11), (8, 12), (10, 14)],
    # 0b0011 : [(0, 8), (1, 6), (3, 14), (4, 5), (7, 12), (9, 11), (10, 13)],
    # 0b1000 : [(0, 9), (1, 3), (2, 13), (4, 7), (5, 12), (6, 14), (8, 11)],
    # 0b1100 : [(0, 10), (1, 5), (2, 11), (3, 12), (4, 6), (7, 14), (8, 13)],
    # 0b1011 : [(0, 11), (1, 14), (2, 10), (3, 6), (4, 12), (5, 7), (8, 9)],
    # 0b0110 : [(0, 12), (2, 6), (3, 10), (4, 11), (5, 9), (7, 8), (13, 14)],
    # 0b1111 : [(0, 13), (1, 4), (2, 9), (3, 7), (5, 6), (8, 10), (12, 14)],
    # 0b1001 : [(0, 14), (1, 11), (2, 5), (3, 8), (6, 9), (7, 10), (12, 13)],
    # 0b0100 : [(1, 12), (2, 8), (3, 5), (4, 14), (6, 7), (9, 10), (11, 13)]}
    # operators_1 = [0b1001, 0b0110, 0b111]
    # nodes_1 = (0, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14)

    # operators_2 = [0b0010, 0b0110, 0b1000, 0b1010, 0b1100, 0b1110]
    # nodes_2 = (2, 6, 7, 8)
    
    # suborbits = split_into_suborbits(family_of_valid_graphs=family_of_valid_graphs_1, operators=operators_2, nodes=nodes_2)
    # #print("suborbits: ", suborbits)
