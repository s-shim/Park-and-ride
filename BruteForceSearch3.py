import random
import math
import pandas as pd
from gurobipy import *
import networkx as nx
from itertools import combinations
import time
import copy
import datetime

N = 8
repNumber = 1000

repArray = []
bestValArray = []
timeArray = []
bestMatrix = {}
for rep in range(repNumber):
    print(datetime.datetime.now())
    
    candidates = pd.read_csv('rep%s_CAN30.csv'%rep)
    cbds = pd.read_csv('rep%s_CBD.csv'%rep)
    nbhs = pd.read_csv('rep%s_NBH.csv'%rep)
    ods = pd.read_csv('rep%s_OD40.csv'%rep)
    
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
            
    tic = time.time()
    bestObjValue = 0.0     
    for comb in combinations(candidates['ID'],N):
        objValue = 0.0
        for odPairID in ods['OD_Pair_ID']:
            totalCoefficient = odCoefficients[odPairID]
            totalCoefficientPP = 0.0
            for canID in list(comb):
                totalCoefficient += oidCoefficients[odPairID,canID]
                totalCoefficientPP += oidCoefficients[odPairID,canID]
            objValue += totalCoefficientPP / totalCoefficient
        if bestObjValue < objValue:
            bestObjValue = objValue
            bestSolution = copy.deepcopy(list(comb))
    toc = time.time()
    bestSolution.sort()
    
    for i in range(len(bestSolution)):
        bestMatrix[rep,i] = bestSolution[i]
            
    print('bestObjValue =', bestObjValue)
    print('bestSolution =', bestSolution)
    print('computational time =',toc - tic)
    print('N =',N)
    print('# candidates =',len(candidates['ID']))
    print('# od pairs =',len(list(ods['OD_Pair_ID'])))
    
    repArray += [rep]
    bestValArray += [bestObjValue]
    timeArray += [toc - tic]
    
    data = pd.DataFrame(list(zip(repArray,bestValArray,timeArray)),columns =['Rep','optimal value','time (s)'])
        
    for i in range(len(bestSolution)):
        arraySolution = []
        for r in range(len(repArray)):
            arraySolution += [bestMatrix[r,i]]            
        data['sol[%s]'%i] = arraySolution

    data.to_csv(r'dataBFS_N%s.csv'%N, index = False)#Check
            
        
