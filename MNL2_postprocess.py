import random
import math
import pandas as pd
from gurobipy import *
import networkx as nx
from itertools import combinations
import time
import copy
import datetime
import myDictionary as md

repNumber = 10

repArray = []
intValueArray = []
exactValueArray = []
percentArray = []
timeArray = []
solutionMatrix = {}
for rep in range(repNumber):
    print(datetime.datetime.now())
    print('Rep =',rep)
    tic = time.time()
    N = 8
    
    candidates = pd.read_csv('rep%s_CAN30.csv'%(rep))
    cbds = pd.read_csv('rep%s_CBD.csv'%(rep))
    nbhs = pd.read_csv('rep%s_NBH.csv'%(rep))
    ods = pd.read_csv('rep%s_OD40.csv'%(rep))
    
    odCoefficients = {}
    oidCoefficients = {}
    
    for odPairID in ods['OD_Pair_ID']:
        [x_origin] = ods.loc[ods['OD_Pair_ID'] == odPairID,'x_origin']
        [y_origin] = ods.loc[ods['OD_Pair_ID'] == odPairID,'y_origin']
        [x_destin] = ods.loc[ods['OD_Pair_ID'] == odPairID,'x_destin']
        [y_destin] = ods.loc[ods['OD_Pair_ID'] == odPairID,'y_destin']
        length_od = math.sqrt((x_origin-x_destin)**2 + (y_origin-y_destin)**2)
        coefficient_od = math.exp(-length_od)
        odCoefficients[odPairID] = coefficient_od
        
        for canID in candidates['ID']:
            [x_i] = candidates.loc[candidates['ID'] == canID, 'xCoordinate']
            [y_i] = candidates.loc[candidates['ID'] == canID, 'yCoordinate']
            length_oid = math.sqrt((x_origin-x_i)**2 + (y_origin-y_i)**2) + math.sqrt((x_i-x_destin)**2 + (y_i-y_destin)**2)
            coefficient_oid = math.exp(-length_oid)
            oidCoefficients[odPairID,canID] = coefficient_oid
    
    # ILP Model
    model = Model('MNL')
    model.setParam('OutputFlag', 0)
    #model.setParam('NumericFocus', 3)
    #model.setParam('IntFeasTol', 1e-9)
    
    ## Employ Variables
    x_vars = []
    x_names = []
    for canID in candidates['ID']:
        x_vars += [(canID)]
        x_names += ['X[%s]'%canID]    
    X = model.addVars(x_vars, vtype = GRB.BINARY, name = x_names)
    
    pp_vars = []
    pp_names = []
    for odPairID in ods['OD_Pair_ID']:
        for canID in candidates['ID']:
            pp_vars += [(odPairID,canID)]
            pp_names += ['PP[%s,%s]'%(odPairID,canID)]
    PP = model.addVars(pp_vars, vtype = GRB.CONTINUOUS, name = pp_names)
            
    pc_vars = []
    pc_names = []
    for odPairID in ods['OD_Pair_ID']:
        pc_vars += [(odPairID)]
        pc_names += ['PC[%s]'%(odPairID)]
    PC = model.addVars(pc_vars, vtype = GRB.CONTINUOUS, name = pc_names)
    
    
    ## Add Constraints
    ### Knapsack
    LHS = []
    for canID in candidates['ID']:
        LHS += [(1,X[canID])]
    model.addConstr(LinExpr(LHS)==N, name='Eq.knapsack')
                    
    ### PP does not exceed X
    for canID in candidates['ID']:
        for odPairID in ods['OD_Pair_ID']:
            LHS = [(1,PP[odPairID, canID]),(-1,X[canID])]
            model.addConstr(LinExpr(LHS)<=0, name='Eq.pp<x(%s,%s)'%(odPairID,canID))
    
    ### Entire Probability Equals 1 for each o-d pair
    for odPairID in ods['OD_Pair_ID']:
        LHS = [(1,PC[odPairID])]
        for canID in candidates['ID']:
            LHS += [(1,PP[odPairID,canID])]
        model.addConstr(LinExpr(LHS)==1, name='Eq.pp+pc=1(%s)'%(odPairID))
    
    ### 1: PP < PC and 2: PC < PP
    for odPairID in ods['OD_Pair_ID']:
        for canID in candidates['ID']:
            
            LHS1 = [(odCoefficients[odPairID],PP[odPairID,canID]),(-oidCoefficients[odPairID,canID],PC[odPairID])]
            model.addConstr(LinExpr(LHS1)<=0, name='Eq.pp<pc(%s,%s)'%(odPairID,canID))
    
            LHS2 = [(oidCoefficients[odPairID,canID],PC[odPairID]),(-odCoefficients[odPairID],PP[odPairID,canID]),(oidCoefficients[odPairID,canID],X[canID])]
            model.addConstr(LinExpr(LHS2)<=oidCoefficients[odPairID,canID], name='Eq.pc<pp(%s,%s)'%(odPairID,canID))
            
        
        
    # Objective Function
    objTerms = []
    for odPairID in ods['OD_Pair_ID']:
        for canID in candidates['ID']:
            objTerms += [(1,PP[odPairID,canID])]
    model.setObjective(LinExpr(objTerms), GRB.MAXIMIZE)
    
    # update and solve the model
    model.update()
    model.optimize()
    
    toc = time.time()
    print('CPU Time =',toc-tic)
    
    obj = model.getObjective()
    optValue = obj.getValue() 
    print('optimal value =',optValue)
      
    # read the optimal solution
    variableName = []
    variableValue = []
    for v in model.getVars():
        if v.x > 0:
            variableName += [v.varname]
            variableValue += [v.x]
    
    optSolution = pd.DataFrame(list(zip(variableName, variableValue)),columns =['varName', 'varVal'])
    #optSolution.to_csv(r'optSolution.csv', index = False)#Check
    
    theSelected = []
    for varName in optSolution['varName']:
        if varName[0] == 'X':
            theSelected += [int(varName[2:-1])]
    theSelected.sort()
    for i in range(len(theSelected)):
        solutionMatrix[rep,i] = theSelected[i]
    
    objValue = md.evaluator(theSelected, odCoefficients, oidCoefficients, ods)
    print('exact value =',objValue)
    percentError = (optValue - objValue)/objValue*100
    print('percent error =',percentError,'%')
    print('theSelected =',theSelected)
    print()
    print()
    
    
    repArray += [rep]
    intValueArray += [optValue]
    exactValueArray += [objValue]
    percentArray += [percentError]
    timeArray += [toc - tic]

    data = pd.DataFrame(list(zip(repArray,intValueArray,exactValueArray,percentArray,timeArray)),columns =['Rep','Int Opt','Exact','Percent Error (%)','CPU Time (s)'])
    data.to_csv(r'data.csv', index = False)#Check

for i in range(len(theSelected)):
    array_i = []
    for rep in repArray:
        array_i += [solutionMatrix[rep,i]]
    data['sol[%s]'%i] = array_i

data.to_csv(r'data.csv', index = False)#Check
    
    
    
    
    
    
    
    
    
    
