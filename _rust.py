from typing import List

class Graph:
    edges: List[List[str]]
    cost: int
    cost_reduced: int

class Problem:
    graphs: List[Graph]

class Solution:
    graphs: List[Graph]
    cost: int

def find_minimal_mixer(problem: Problem, reduced: bool) -> Solution: # type: ignore
    pass

def find_minimal_mixer_enumeration(problem: Problem, reduced: bool) -> Solution: # type: ignore
    pass
