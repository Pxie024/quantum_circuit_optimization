#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 21:44:23 2021

@author: paulxie
"""


import numpy as np
from gurobipy import *
import csv
import pandas as pd
import ast
from mns_nnc_solver import mns_solver

ARCS = [(1,2), 
        (2,3), 
        (3,4), 
        (4,5)]

OPERATIONS = {
23:    [2,      1,       11,      1],
26:    [1,      2,       12,      2],
28:    [2,      1,       11,      3],
33:    [0,      3,       11,      1],
36:    [3,      0,       12,      2],
38:    [0,      3,       11,      3],
43:    [2,      1,       11,      4],
46:    [1,      2,       12,      5],
48:    [2,      1,       11,      6],
53:    [0,      4,       11,      4],
56:    [4,      0,       12,      5],
58:    [0,      4,       11,      6],
63:    [0,      3,       11,      7],
66:    [3,      0,       12,      8],
68:    [0,      3,       11,      9],
73:    [1,      2,       11,      7],
76:    [2,      1,       12,      8],
78:    [1,      2,       11,      9],
83:    [0,      3,       11,      10],
86:    [3,      0,       12,      11],
88:    [0,      3,       11,      12],
93:    [4,      1,       11,      10],
96:    [1,      4,       12,      11],
98:    [4,      1,       11,      12],
103:   [3,      1,       11,      13],
106:   [1,      3,       12,      14],
108:   [3,      1,       11,      15],
113:   [4,      2,       11,      13],
116:   [2,      4,       12,      14],
118:   [4,      2,       11,      15]
}


soln = mns_solver(ARCS, OPERATIONS, max_instance = 7, K_dir = 10, S = 34)
m,y_dir,y_rev,x,z = soln.solve_it(TimeLimit = 15800)
soln.print_results(m,y_dir,y_rev,x,z)


# Minimum objective value is  386.0
# Solution time is 144.9 sec