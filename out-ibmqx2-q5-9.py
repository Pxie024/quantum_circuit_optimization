#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 21:57:44 2021

@author: paulxie
"""


import numpy as np
from gurobipy import *
import csv
import pandas as pd
import ast
from mns_nnc_solver import mns_solver

ARCS = [(1,2),
        (1,3),
        (2,3), 
        (4,3), 
        (4,5),
        (5,3)]

OPERATIONS = {
23:    [0,      3,       11,      1],
26:    [3,      0,       12,      2],
28:    [0,      3,       11,      3],
33:    [1,      2,       11,      1],
36:    [2,      1,       12,      2],
38:    [1,      2,       11,      3],
43:    [2,      1,       11,      4],
46:    [1,      2,       12,      5],
48:    [2,      1,       11,      6],
53:    [0,      4,       11,      4],
56:    [4,      0,       12,      5],
58:    [0,      4,       11,      6],
63:    [2,      4,       11,      7],
66:    [4,      2,       12,      8],
68:    [2,      4,       11,      9],
73:    [3,      1,       11,      7],
76:    [1,      3,       12,      8],
78:    [3,      1,       11,      9],
83:    [4,      2,       11,      10],
86:    [2,      4,       12,      11],
88:    [4,      2,       11,      12],
93:    [3,      0,       11,      10],
96:    [0,      3,       12,      11],
98:    [3,      0,       11,      12],
103:   [0,      4,       11,      13],
106:   [4,      0,       12,      14],
108:   [0,      4,       11,      15],
113:   [3,      2,       11,      13],
116:   [2,      3,       12,      14],
118:   [3,      2,       11,      15]

}


soln = mns_solver(ARCS, OPERATIONS, max_instance = 4, K_dir = 10, S = 34)
m,y_dir,y_rev,x,z = soln.solve_it(TimeLimit = 10800)
soln.print_results(m,y_dir,y_rev,x,z)


# Minimum objective value is  420.0
# Solution time is 933.2 sec