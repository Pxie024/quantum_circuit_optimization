#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 11:45:02 2021

@author: paulxie
"""


import numpy as np
from gurobipy import *
from networkx import DiGraph
import networkx as nx




class mns_solver:
    def __init__(self, arcs, operations, max_instance, K_dir, S, K_rev = None):
        G = DiGraph()
        for a in arcs:
            G.add_edge(a[0], a[1])
        self.V = list(G.nodes())
        self.A = list(G.edges())
        self.A_rev = []
        for a in self.A:
            self.A_rev.append((a[1], a[0]))
        self.gates_layout = []
        self.C = []
        self.C_id = list(operations.keys())
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
            self.gates_layout.append([])
            for o in operations:
                if operations[o][-1] == l:
                    self.gates_layout[-1].append((operations[o][0], operations[o][1]))
        self.max_inst = range(max_instance) 
        if K_rev != None:
            for k in self.K_rev.keys():
                self.K_rev[k] = K_rev
        self.K_dir = K_dir
        self.S = S
        self.time_range = np.arange(1, len(self.Layers)*len(self.max_inst)+1)
        self.first_inst = (np.array(self.Layers)-1) * len(self.max_inst) + 1
        self.operations = operations
        
        
    def detail_time(self, t):
        if t % len(self.max_inst) == 0:
            i = len(self.max_inst) # which instance the given time belongs to 
            h = t // len(self.max_inst) - 1 # which layer the time instant belongs to
        else:
            i = t % len(self.max_inst)
            h = t // len(self.max_inst) # which layer the time instant belongs to
       
        
        return h+1 , i
    
    
    
    def incident_arcs(self, node):
        inc_arcs = []
        for a in self.A:
            if node in a:
                inc_arcs.append(a)
        return inc_arcs
    
    
    
    def solve_it(self, TimeLimit = 10800, init_cond = False):  # max computing time set to 3 hours 
        m = Model('swap gates')
        y_dir = m.addVars(self.A, self.C_id, self.time_range, vtype = GRB.BINARY, name = "direct_imp")
        y_rev = m.addVars(self.A_rev, self.C_id, self.time_range, vtype = GRB.BINARY, name = "reverse_imp")
        x = m.addVars(self.Q, self.V, self.time_range, vtype = GRB.BINARY, name = "qubit_position")
        z = m.addVars(self.A, self.time_range, vtype = GRB.BINARY, name = "swap_gates")
        
        m._cost_1 = 0
        for a in self.A:
            for c in self.C_id:
                for t in self.time_range:
                    m._cost_1 += self.K_dir * y_dir[a[0],a[1],c,t] + self.K_rev[c] * y_rev[a[1],a[0],c,t]
        
        
        # cost_1 = quicksum(quicksum(quicksum([self.K_dir * y_dir[a[0],a[1],c,t], self.K_rev[c] * y_rev[a[1],a[0],c,t]] 
        #                                     for t in self.time_range) for c in self.C_id) for a in self.A)
        m._cost_2 = quicksum(quicksum(self.S * z[a[0],a[1],t] for t in self.time_range) for a in self.A)
        m.setObjective(m._cost_1 + m._cost_2, GRB.MINIMIZE)
        
        for t in self.time_range:    
            m.addConstrs((quicksum(x[q, i, t] for q in self.Q) <= 1 for i in self.V), 'constr_5')
            m.addConstrs((quicksum(x[q, i, t] for i in self.V) == 1 for q in self.Q), 'constr_6')
            
        for t in self.time_range[:-1]:
            for i in self.V:
                for j in self.V:
                    if (i, j) in self.A:
                        m.addConstrs((x[q,j,t+1] + x[q,i,t] <= 1 + z[i,j,t] for q in self.Q), 'constr_7')
                        m.addConstrs((x[q,i,t+1] + x[q,j,t] <= 1 + z[i,j,t] for q in self.Q), 'constr_8')
                    elif (j, i) not in self.A and i != j:    
                        m.addConstrs((x[q,i,t+1] + x[q,j,t] <= 1 for q in self.Q), 'constr_9')
                        
        for t_1 in self.first_inst:
            C_h_id = []
            for o in self.operations:
                if (self.operations[o][-1]) == (t_1-1) / len(self.max_inst) + 1:
                    C_h_id.append(o)
            #print(t_1, C_h_id)
            
        # for h in self.Layers:
        #     t_1 = h * len(self.max_inst)
        #     C_h = self.gates_layout[h]
            # if type(C_h) != list:
            #     C_h = [C_h]  ## will need some update if any layer is composed of multiple gates 
            for (i, j) in self.A:
                m.addConstrs((x[self.operations[c][0],i,t_1] + x[self.operations[c][1],j,t_1] >= 2*y_dir[i,j,c,t_1] for c in C_h_id), 'constr_2')
                m.addConstrs((x[self.operations[c][0],j,t_1] + x[self.operations[c][1],i,t_1] >= 2*y_rev[j,i,c,t_1] for c in C_h_id), 'constr_3')
            for c in C_h_id:
                # m.addConstrs((quicksum(y_dir[a[0],a[1],c,t_1]+ y_rev[a[1],a[0],c,t_1] for a in self.A) == 1), 'constr_4')
                
                val = 0
                for (i, j) in self.A:
                    val += (y_dir[i,j,c,t_1] + y_rev[j,i,c,t_1])
                #print(val)
                m.addConstr((val == 1), 'constr_4')
            
            t_h = range(t_1 + 1, t_1 + len(self.max_inst))
            for t in t_h:
                for (i, j) in self.A:
                    delta_i = self.incident_arcs(i)
                    delta_j = self.incident_arcs(j)
                    term_1 = quicksum(z[a,b,t-1] for (a,b) in delta_i)
                    term_2 = quicksum(z[c,d,t-1] for (c,d) in delta_j)
                    m.addConstr((z[i,j,t] <= quicksum([term_1, term_2])), 'constr_10')
            
        if init_cond == True:
            m.addConstrs((x[i,i,1] == 1 for i in self.Q), 'intial_condition')
        
        m.Params.TimeLimit = TimeLimit
        m.optimize()
        return m, y_dir, y_rev, x, z
                
                
        
        
        
    def print_results(self, m, y_dir, y_rev, x, z):
        print('\n*********************************** Summary *************************************')
        if m.status == GRB.OPTIMAL:
            print('\nMinimum objective value is ', round(m.ObjVal, 3))
            print('\nSolution time is {} sec'.format(round(m.RunTime, 1)))
            
            # quibit locations
            for t in self.time_range:
                h, ins = self.detail_time(t)
                if ins == 1:
                    print('\nAt first instance of layer {}:'.format(h))
                    for i in self.Q:
                        for j in self.V :
                            if x[i,j,t].x == 1:
                                print('\n    qubit {} is at location {}'.format(i, j))
                
            # swap gates implementation
            for t in self.time_range:
                for a in self.A:
                    if z[a[0],a[1],t].x == 1:
                        h, ins = self.detail_time(t)
                        print('\nA SWAP gate is implemented on arc {} at instance {} of layer {}'.format(a, ins, h))
                        
                        
            # 2-qubit quantum gates 
            for t in self.time_range:
                for a in self.A:
                    for c in self.C_id:
                        if y_dir[a[0],a[1],c,t].x == 1:
                            h, ins = self.detail_time(t)
                            print('\nGate {} is directly implemented on arc {} at instance {} of layer {}'.format(c, a, ins, h))
                        if y_rev[a[1],a[0],c,t].x == 1:
                            h, ins = self.detail_time(t)
                            print('\nGate {} is reversely implemented on arc {} at instance {} of layer {}'.format(c, a, ins, h))
                            
        else:
            print('\nThe model is unsolved or infeasible.')
                            


    def get_qubits_locations(self, x):  #return the qubit locations at each time instance 
        self.locations = {}
        for q in self.Q:
            self.locations[q] = []
            for t in self.time_range:
                for i in self.V:
                    if x[q,i,t].x == 1:
                        self.locations[q].append(i)
        return self.locations
    
    












