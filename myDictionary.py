import random
import math
import pandas as pd
from gurobipy import *
import networkx as nx
from itertools import combinations
import time
import copy
import datetime


def evaluator(theSelected, odCoefficients, oidCoefficients,ods):
    objValue = 0.0
    for odPairID in ods['OD_Pair_ID']:
        totalCoefficient = odCoefficients[odPairID]
        totalCoefficientPP = 0.0
        for canID in theSelected:
            totalCoefficient += oidCoefficients[odPairID,canID]
            totalCoefficientPP += oidCoefficients[odPairID,canID]
        objValue += totalCoefficientPP / totalCoefficient
    return objValue


def constructMatrix(candidates,ods):
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

    return odCoefficients, oidCoefficients 

        
