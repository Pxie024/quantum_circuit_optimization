#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 13:42:23 2021

@author: paulxie
"""


import numpy as np
from gurobipy import *
from networkx import DiGraph
import networkx as nx




class kt_solver:
    def __init__(self, arcs, operations, K_dir, S, K_rev = None):
        G = DiGraph()
        for a in arcs:
            G.add_edge(a[0], a[1])
        self.V = list(G.nodes())
        self.A = list(G.edges())
        self.A_rev = []
        for a in self.A:
            self.A_rev.append((a[1], a[0]))
        self.gates_layout = {}
        self.C = []
        self.gates_id = list(operations.keys())
        self.K_rev = {}
        self.Q = []
        Layers = []
        for o in operations:
            Layers.append(operations[o][-1])
            self.K_rev[o] = operations[o][2]
            self.C.append((operations[o][0], operations[o][1]))
            if operations[o][0] not in self.Q:
                self.Q.append(operations[o][0])
            if operations[o][1] not in self.Q:
                self.Q.append(operations[o][1])
        self.Q.sort()
        self.Layers = range(min(Layers), max(Layers)+1)
        for l in self.Layers:
            self.gates_layout[l] = []
            for o in operations:
                if operations[o][-1] == l:
                    self.gates_layout[l].append(o)
        if K_rev != None:
            for k in self.K_rev.keys():
                self.K_rev[k] = K_rev
        self.K_dir = K_dir
        self.S = S
        self.operations = operations
        
        
    def solve_it(self, TimeLimit = 10800):
        M = len(self.V) + 2
        
        m = Model('KT_Formulation')
        
        x = m.addVars(self.Q, self.Layers, lb = min(self.V), ub = max(self.V), vtype = GRB.INTEGER, name = 'x')
        z = m.addVars(self.Q, self.Q, self.Layers, vtype = GRB.BINARY, name = 'z')
        # for i,j,t in z.keys():
        #     z[i,j,t] = -z[j,i,t]
            
        y_dir = m.addVars(self.A, self.gates_id, self.Layers, vtype = GRB.BINARY, name ='y')
        y_rev = m.addVars(self.A_rev, self.gates_id, self.Layers, vtype = GRB.BINARY, name = 'y*')
        r = m.addVars(self.Q, self.V, self.Layers, vtype = GRB.BINARY, name = 'r')
        k = m.addVars(self.Q, self.Q, self.Layers[:-1], vtype = GRB.BINARY, name = 'k')
        # for i,j,t in k.keys():
        #     k[i,j,t] = k[j,i,t]
        
        m._cost_1 = 0
        for a in self.A:
            for c in self.gates_id:
                for t in self.Layers:
                    m._cost_1 += self.K_dir * y_dir[a[0],a[1],c,t] + self.K_rev[c] * y_rev[a[1],a[0],c,t]
                    
                    
        m._cost_2 = 0
        for i in self.Q:
            for j in self.Q:
                if i != j:
                    for t in self.Layers[:-1]:
                        m._cost_2 += 0.5 * self.S * k[i,j,t]
        
        m.setObjective(m._cost_1 + m._cost_2, GRB.MINIMIZE)
        
        
        for t in self.Layers:
            for c in self.gates_layout[t]:
                p,q = self.operations[c][0], self.operations[c][1]
                for (i,j) in self.A:
                    m.addConstr((r[p,i,t] + r[q,j,t] >= 2*y_dir[i,j,c,t]), 'constr_1')
                    m.addConstr((r[p,j,t] + r[q,i,t] >= 2*y_rev[j,i,c,t]), 'constr_2')
                    
                # m.addConstr((x[p,t] - x[j,t] >= -1), 'constr_6_L')
                # m.addConstr((x[p,t] - x[j,t] <= 1), 'constr_6_R')
                    
                m.addConstr((sum(y_dir[a[0],a[1],c,t]+ y_rev[a[1],a[0],c,t] for a in self.A) == 1), 'constr_3')
                
            for i in self.Q:
                for j in self.Q:
                    if i != j:
                        m.addConstr((x[i,t] - x[j,t] <= M*z[i,j,t] -1), 'constr_4')
                        # m.addConstr((x[j,t] - x[i,t] <= M*(1-z[i,j,t]) -1), 'constr_5')
                        
            for q in self.Q:
                m.addConstr((x[q,t] == quicksum(i* r[q,i,t] for i in self.V)), 'constr_7')
                m.addConstr((quicksum(r[q,i,t] for i in self.V) == 1), 'constr_10')
                
            m.addConstrs((quicksum(r[q,i,t] for q in self.Q) <= 1 for i in self.V), 'constr_9')
            
            for i in self.Q:
                for j in self.Q:
                    if i != j :
                        m.addConstr((z[i,j,t] + z[j,i,t] == 1), 'constr_11')
            
                
        
        for t in self.Layers[:-1]:
            for i in self.Q:
                for j in self.Q:
                    if i != j :
                        m.addConstr((z[i,j,t] - z[i,j,t+1] >= -k[i,j,t]), 'constr_8_L')
                        m.addConstr((z[i,j,t] - z[i,j,t+1] <= k[i,j,t]), 'constr_8_R')
        
        
        m.Params.TimeLimit = TimeLimit
        m.optimize()
        
        return m, x, z, y_dir, y_rev, r, k
        
        
    def print_results(self, m, x, z, y_dir, y_rev, r, k):
        if m.status != GRB.OPTIMAL:
            print('\nThe model is infeasible or was not solved to optimality.')
        else:
            print('\n')
            print('Optimal objective value is {}'.format(m.ObjVal, 2))
            swap_count = 0
            for q1, q2, t in k.keys():
                if (k[q1,q2,t].x == 1) and (q1<q2):
                    swap_count += 1
                    print('A swap gate was implemented on qubits {} between layer {} and {}'.format((q1,q2), t, t+1))
            
            dir_count = 0
            m_value = []  # record the m value needed for each layer 
            for t in self.Layers[:-1]:
                steps = {1:[]}  # track which qubits were swaped at each step 
                for q1 in self.Q:
                    for q2 in self.Q:
                        if q1 < q2:
                            max_inst = max(steps.keys())
                            if (k[q1,q2,t].x == 1):
                                if (q1 not in steps[max_inst]) and (q2 not in steps[max_inst]):
                                    steps[max_inst].extend([q1, q2])
                                else:
                                    steps[max_inst + 1] = [q1, q2]
                max_inst = max(steps.keys())
                if steps[max_inst] != []:
                    m_value.append(max_inst)
                    print('{} steps is needed between layers {} and {}'.format(max_inst, t, t+1))
                            
                    
                
            for t in self.Layers:
                for a in self.A:
                    for c in self.gates_id:
                        if y_dir[a[0], a[1], c, t].x == 1:
                            print('Gate {} is directly implemented on arc {} at layer {}'.format(c, a, t))
                            dir_count += 1
                        elif y_rev[a[1], a[0], c, t].x == 1:
                            print('Gate {} is reversely implemented on arc {} at layer {}'.format(c, a, t))
                            
            print('\nGates directly implemented: {}/{}'.format(dir_count, len(self.operations)))
            print('Number of swap gates used: {}'.format(swap_count))
            print('m value needed to for feasible soluiton: ', max(m_value))


    def get_m_value(self, k):
        m_value = []  # record the m value needed for each layer 
        for t in self.Layers[:-1]:
            steps = {1:[]}  # track which qubits were swaped at each step 
            for q1 in self.Q:
                for q2 in self.Q:
                    if q1 < q2:
                        max_inst = max(steps.keys())
                        if (k[q1,q2,t].x == 1):
                            if (q1 not in steps[max_inst]) and (q2 not in steps[max_inst]):
                                steps[max_inst].extend([q1, q2])
                            else:
                                steps[max_inst + 1] = [q1, q2]
            max_inst = max(steps.keys())
            if steps[max_inst] != []:
                m_value.append(max_inst)
            else:
                m_value.append(1)
        
        return max(m_value)
        





ARCS = [(1,3), (4,1), (2,4), (5,2)]

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






soln = kt_solver(ARCS, OPERATIONS, 10, 34)
m, x, z, y_dir, y_rev, r, k = soln.solve_it()
soln.print_results(m, x, z, y_dir, y_rev, r, k)
m_val = soln.get_m_value(k)
print('Runtime = ', m.runtime)