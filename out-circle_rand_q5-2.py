#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 14:54:00 2021

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


ARCS = [(1,2), (3,4), (3,1), (4,5), (2,5)]
# ARCS = [(1,2), (2,5), (4,5), (3,4)]

OPERATIONS = {
23:    [3,      2,       11,      1],
26:    [2,      3,       12,      2],
28:    [3,      2,       11,      3],
33:    [0,      1,       11,      1],
36:    [1,      0,       12,      2],
38:    [0,      1,       11,      3],
43:    [2,      4,       11,      4],
46:    [4,      2,       12,      5],
48:    [2,      4,       11,      6],
53:    [1,      0,       11,      4],
56:    [0,      1,       12,      5],
58:    [1,      0,       11,      6],
63:    [4,      3,       11,      7],
66:    [3,      4,       12,      8],
68:    [4,      3,       11,      9],
73:    [1,      0,       11,      7],
76:    [0,      1,       12,      8],
78:    [1,      0,       11,      9],
83:    [2,      0,       11,     10],
86:    [0,      2,       12,      11],
88:    [2,      0,       11,      12],
93:    [1,      3,       11,      10],
96:    [3,      1,       12,      11],
98:    [1,      3,       11,      12],
103:   [0,      2,       11,      13],
106:   [2,      0,       12,      14],
108:   [0,      2,       11,      15],
113:   [4,      3,       11,      13],
116:   [3,      4,      12,      14],
118:   [4,      3,       11,      15]
}


soln = mns_solver(ARCS, OPERATIONS, max_instance = 7, K_dir = 10, S = 34)
m,y_dir,y_rev,x,z = soln.solve_it()
soln.print_results(m,y_dir,y_rev,x,z)


# Minimum objective value is  353.0
# Solution time is 1593.5 sec


# check = mns_checker(y_dir, y_rev, x, z, ARCS, OPERATIONS, 2)
# loc = check.get_qubits_locations()
# check.swap_gates_checker(loc)
# check.two_qubit_gates_checker(loc)
# check.qubit_location_checker(loc)
# check.first_inst_checker()
# check.adjacency_checker(loc)


# 1 -- 2 -- 3  -- 4 -- 5


