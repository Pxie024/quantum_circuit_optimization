#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 15:06:30 2021

@author: paulxie
"""


import numpy as np
from gurobipy import *
import csv
import pandas as pd
import ast
from mns_nnc_solver import mns_solver
from mns_nnc_checker import mns_checker


# read .dat to a list of lists
# datContent = [i.split() for i in open("5qubit_instances/out-circle_rand_q5-2.dat").readlines()]


ARCS = [(1,2), 
        (1,3), 
        (2,3), 
        (4,3),
        (4,5),
        (5,3)]

# ARCS = [(1,2),  #reduced linear graph
#         (2,3),
#         (4,3),
#         (4,5)]


OPERATIONS = {
23:    [0,      4,       11,      1],
26:    [4,      0,       12,      2],
28:    [0,      4,       11,      3],
33:    [3,      2,       11,      1],
36:    [2,      3,       12,      2],
38:    [3,      2,       11,      3],
43:    [2,      4,       11,      4],
46:    [4,      2,       12,      5],
48:    [2,      4,       11,      6],
53:    [3,      0,       11,      4],
56:    [0,      3,       12,      5],
58:    [3,      0,       11,      6],
63:    [0,      4,       11,      7],
66:    [4,      0,       12,      8],
68:    [0,      4,       11,      9],
73:    [3,      2,       11,      7],
76:    [2,      3,       12,      8],
78:    [3,      2,       11,      9],
83:    [1,      4,       11,      10],
86:    [4,      1,       12,      11],
88:    [1,      4,       11,      12],
93:    [3,      2,       11,      10],
96:    [2,      3,       12,      11],
98:    [3,      2,       11,      12],
103:   [2,      4,       11,      13],
106:   [4,      2,       12,      14],
108:   [2,      4,       11,      15],
113:   [3,      1,       11,      13],
116:   [1,      3,       12,      14],
118:   [3,      1,       11,      15]


}


soln = mns_solver(ARCS, OPERATIONS, max_instance = 4, K_dir = 10, S = 34)
m,y_dir,y_rev,x,z = soln.solve_it()
soln.print_results(m,y_dir,y_rev,x,z)




# check = mns_checker(y_dir, y_rev, x, z, ARCS, OPERATIONS, 2)
# loc = check.get_qubits_locations()
# check.swap_gates_checker(loc)
# check.two_qubit_gates_checker(loc)
# check.qubit_location_checker(loc)
# check.first_inst_checker()
# check.adjacency_checker(loc)
