#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 21:11:32 2021

@author: paulxie
"""


import numpy as np
from gurobipy import *
import csv
import pandas as pd
import ast
from mns_nnc_solver import mns_solver

# read .dat to a list of lists
# datContent = [i.split() for i in open("5qubit_instances/out-circle_rand_q5-2.dat").readlines()]


ARCS = [(1,2), 
        (3,2), 
        (3,1), 
        (3,4), 
        (5,3),
        (5,4)]

OPERATIONS = {
23:    [3,      4,       11,      1],
26:    [4,      3,       12,      2],
28:    [3,      4,       11,      3],
33:    [0,      1,       11,      1],
36:    [1,      0,       12,      2],
38:    [0,      1,       11,      3],
43:    [2,      0,       11,      4],
46:    [0,      2,       12,      5],
48:    [2,      0,       11,      6],
53:    [4,      1,       11,      4],
56:    [1,      4,       12,      5],
58:    [4,      1,       11,      6],
63:    [4,      3,       11,      7],
66:    [3,      4,       12,      8],
68:    [4,      3,       11,      9],
73:    [0,      2,       11,      7],
76:    [2,      0,       12,      8],
78:    [0,      2,       11,      9],
83:    [4,      0,       11,      10],
86:    [0,      4,       12,      11],
88:    [4,      0,       11,      12],
93:    [2,      1,       11,      10],
96:    [1,      2,       12,      11],
98:    [2,      1,       11,      12],
103:   [2,      1,       11,      13],
106:   [1,      2,       12,      14],
108:   [2,      1,       11,      15],
113:   [4,      0,       11,      13],
116:   [0,      4,       12,      14],
118:   [4,      0,       11,      15]

}


soln = mns_solver(ARCS, OPERATIONS, max_instance = 4, K_dir = 10, S = 34)
m,y_dir,y_rev,x,z = soln.solve_it(TimeLimit = 15800)
soln.print_results(m,y_dir,y_rev,x,z)


# Minimum objective value is  353.0
# Solution time is 531.5 sec
