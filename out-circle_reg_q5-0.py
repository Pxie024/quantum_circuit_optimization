#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 21:13:58 2021

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
        (2,3), 
        (3,4), 
        (4,5), 
        (5,1)]


OPERATIONS = {
23:    [0,      1,       11,      1],
26:    [1,      0,       12,      2],
28:    [0,      1,       11,      3],
33:    [2,      3,       11,      1],
36:    [3,      2,       12,      2],
38:    [2,      3,       11,      3],
43:    [1,      2,       11,      4],
46:    [2,      1,       12,      5],
48:    [1,      2,       11,      6],
53:    [0,      4,       11,      4],
56:    [4,      0,       12,      5],
58:    [0,      4,       11,      6],
63:    [3,      0,       11,      7],
66:    [0,      3,       12,      8],
68:    [3,      0,       11,      9],
73:    [1,      2,       11,      7],
76:    [2,      1,       12,      8],
78:    [1,      2,       11,      9],
83:    [2,      4,       11,      10],
86:    [4,      2,       12,      11],
88:    [2,      4,       11,      12],
93:    [0,      1,       11,      10],
96:    [1,      0,       12,      11],
98:    [0,      1,       11,      12],
103:   [1,      2,       11,      13],
106:   [2,      1,       12,      14],
108:   [1,      2,       11,      15],
113:   [3,      0,       11,      13],
116:   [0,      3,       12,      14],
118:   [3,      0,       11,      15]

}


soln = mns_solver(ARCS, OPERATIONS, max_instance = 7, K_dir = 10, S = 34)
m,y_dir,y_rev,x,z = soln.solve_it()
soln.print_results(m,y_dir,y_rev,x,z)


# Minimum objective value is  354.0
# Solution time is 1197.9 sec


check = mns_checker(y_dir, y_rev, x, z, ARCS, OPERATIONS, 2)
loc = check.get_qubits_locations()
check.swap_gates_checker(loc)
check.two_qubit_gates_checker(loc)
check.qubit_location_checker(loc)
check.first_inst_checker()
check.adjacency_checker(loc)
