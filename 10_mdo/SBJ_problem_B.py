import os
from DMDO import *
import math
import numpy as np
from numpy import cos, exp, pi, prod, sin, sqrt, subtract, inf
import copy
from typing import List, Dict, Any, Callable, Protocol, Optional

from SBJ_analysis import SBJ_A1, SBJ_A2, SBJ_A3, SBJ_A4, Wing_Mod, loads

user =USER

def get_inputs(V:List[variableData], names:List[str], sp_index:int):
  """This function returns a list of variable objects based on name and subproblem index"""
  variables = [vi for vi in V if vi.name in names and vi.sp_index==sp_index and vi.coupling_type not in [COUPLING_TYPE.FEEDFORWARD,COUPLING_TYPE.DUMMY]]
  ordered_variables = []
  for name in names:
    for vi in variables:
      if vi.name == name:
        ordered_variables += [vi]
  return ordered_variables

def get_outputs(V:List[variableData], names:List[str], sp_index:int):
  """This function returns a list of variable objects based on name and subproblem index"""
  variables = [vi for vi in V if vi.name in names and vi.sp_index==sp_index and vi.coupling_type in [COUPLING_TYPE.FEEDFORWARD,COUPLING_TYPE.DUMMY]]
  ordered_variables = []
  for name in names:
    for vi in variables:
      if vi.name == name:
        ordered_variables += [vi]
  return ordered_variables

def SBJ_opt1(x, y):

  inputs = {}
  for input_name,value in zip(["SFC","We","LD","Ws","Wf"],x): # these are your targets
    inputs[input_name] = value

  outputs = {}
  
  for outputs_name,value in zip(["Wt", "range"],y): # these are your responses
    outputs[outputs_name] = value

  # SBJ_constraint_range
  G = -outputs["range"]/2000 +1

  return [outputs["Wt"], [G]] # one objective, one constraint

def SBJ_opt2(x, y):
  # SBJ_constraint_power
  inputs = {}
  for input_name,value in zip(["D","T"],x): # these are your targets
    inputs[input_name] = value

  outputs = {}
  for outputs_name,value in zip(["SFC","We","ESF","Temp_E","Throttle_uA"],y): # these are your responses
    outputs[outputs_name] = value

  # -----THIS SECTION COMPUTES POLYNOMIAL CONSTRAINT FUNCTIONS-----
  Dim_Throttle = inputs["T"]*16168
  Temp_uA=1.02
  g1 = outputs["Temp_E"] /Temp_uA-1
  g2 = Dim_Throttle/outputs["Throttle_uA"]-1

  return[0, [g1, g2]] # one objective, two constraints

def SBJ_opt3(x, y):
  inputs = {}
  for input_name,value in zip(["Wt","ESF","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],x): # these are your targets
    inputs[input_name] = value

  outputs = {}
  for outputs_name,value in zip(["L", "D", "LD", "Pg", "CLo1", "CLo2"],y): # these are your responses
    outputs[outputs_name] = value

  # %------Constraints------%
  Pg_uA=1.1
  if outputs["CLo1"] > 0:
      g2=(2*(outputs["CLo2"]))-(outputs["CLo1"])
      g3=(2*(-outputs["CLo2"]))-(outputs["CLo1"])
  else:
      g2=(2*(-outputs["CLo2"]))-(outputs["CLo1"])
      g3=(2*(outputs["CLo2"]))-(outputs["CLo1"])

  g1=outputs["Pg"]/Pg_uA-1

  return [0, [g1, g2, g3]] # no objective, three constraints

def SBJ_opt4(x, y):
  inputs = {}
  for input_name,value in zip(["L","tc","ARw","LAMBDAw","Sref","Sht","ARht","lambda"],x): # these are your targets
    inputs[input_name] = value

  outputs = {}
  for outputs_name,value in zip(["Ws", "Wf", "theta"],y): # these are your responses
    outputs[outputs_name] = value

  index_t = len(inputs) # get second last 9 entries
  index_ts = len(inputs)+9 # get last 9 entries
  
  Z = [inputs["tc"],user.h,user.M,inputs["ARw"],inputs["LAMBDAw"],inputs["Sref"],inputs["Sht"],inputs["ARht"]]
  L = inputs["L"]
  # %----Inputs----%
  t= np.divide(x[index_t:index_t+9], 12.)#convert to feet
  ts= np.divide(x[index_ts:], 12.);#convert to feet
  LAMBDA = inputs["lambda"]

  # SBJ_constraint_weight
  # %----Inputs----%
  # if isinstance(t, np.ndarray):
  #   t=t/12 #convert to feet
  #   ts=ts/12 #convert to feet
  # else:
  #   t=np.array(t)/12 #convert to feet
  #   ts=np.array(ts)/12 #convert to feet

  t1=t[0:3] 
  t2=t[3:6]
  t3=t[6:9]
  ts1=ts[0:3]
  ts2=ts[3:6]
  ts3=ts[6:9]
  G=4000000*144
  E=10600000*144
  nu=0.3

  # %----Inputs----%

  beta=.9
  [c,c_box,Sweep_40,D_mx,b,a]=Wing_Mod(Z,LAMBDA)
  teq1=((t1**3)/4+(3*t1)*(ts1-t1/2)**2)**(1/3)
  teq2=((t2**3)/4+(3*t2)*(ts2-t2/2)**2)**(1/3)
  teq3=((t3**3)/4+(3*t3)*(ts3-t3/2)**2)**(1/3)
  l=0.6*c_box
  h=Z[0]*(beta*np.array(c[0:3]))-(0.5)*(ts1+ts3)
  A_top=(0.5)*(t1*l)+(1/6)*(t2*h)
  A_bottom=(0.5)*(t3*l)+(1/6)*(t2*h)
  Y_bar=h*(2*A_top)/(2*A_top+2*A_bottom)
  Izz=2*A_top*(h-Y_bar)**2+2*A_bottom*(-Y_bar)**2
  P,Mz,Mx,_,_ =loads(b,c,Sweep_40,D_mx,L,Izz,Z,E)

  sig_1=Mz*(0.95*h-Y_bar)/Izz
  sig_2=Mz*(h-Y_bar)/Izz
  sig_3=sig_1
  sig_4=Mz*(0.05*h-Y_bar)/Izz
  sig_5=Mz*(-Y_bar)/Izz
  sig_6=sig_4
  q=Mx/(2*l*h)

  # %-----THIS SECTION COMPUTES THE TOTAL WEIGHT OF A/C-----%

  k=6.09375


  # %----Point 1----%
  T1=P*(l-a)/l
  tau1_T=T1/(h*t2)
  tau1=q/t2+tau1_T
  sig_eq1=sqrt(sig_1**2+3*tau1**2)
  sig_cr1=((pi**2)*E*4/(12*(1-nu**2)))*(teq2/(0.95*h))**2
  tau_cr1=((pi**2)*E*5.5/(12*(1-nu**2)))*(teq2/(0.95*h))**2
  G = np.zeros((72))
  G[0:3]=k*sig_eq1
  G[3:6]=k*(((sig_1)/sig_cr1)+(tau1/tau_cr1)**2)
  G[6:9]=k*((-(sig_1)/sig_cr1)+(tau1/tau_cr1)**2)


  # %----Point 2----%

  tau2=q/t1
  sig_eq2=sqrt(sig_2**2+3*tau2**2)
  sig_cr2=((pi**2)*E*4/(12*(1-nu**2)))*(teq1/l)**2
  tau_cr2=((pi**2)*E*5.5/(12*(1-nu**2)))*(teq1/l)**2
  G[9:12]=k*sig_eq2
  G[12:15]=k*(((sig_2)/sig_cr2)+(tau2/tau_cr2)**2)
  G[15:18]=k*((-(sig_2)/sig_cr2)+(tau2/tau_cr2)**2)

  # %----Point 3----%

  T2=P*a/l
  tau3_T=-T2/(h*t2)
  tau3=q/t2+tau3_T
  sig_eq3=sqrt(sig_3**2+3*tau3**2)
  sig_cr3=sig_cr1
  tau_cr3=tau_cr1
  G[18:21]=k*sig_eq3
  G[21:24]=k*(((sig_3)/sig_cr3)+(tau3/tau_cr3)**2)
  G[24:27]=k*(((-sig_3)/sig_cr3)+(tau3/tau_cr3)**2)

  # %----Point 4----%

  tau4=-q/t2+tau1_T
  sig_eq4=sqrt(sig_4**2+3*tau4**2)
  G[27:30]=k*sig_eq4

  # %----Point 5----%

  tau5=q/t3
  sig_eq5=sqrt(sig_5**2+3*tau5**2)
  sig_cr5=((pi**2)*E*4/(12*(1-nu**2)))*(teq3/l)**2
  tau_cr5=((pi**2)*E*5.5/(12*(1-nu**2)))*(teq3/l)**2
  G[30:33]=k*sig_eq5
  G[33:36]=k*(((sig_5)/sig_cr5)+(tau5/tau_cr5)**2)
  G[36:39]=k*(((-sig_5)/sig_cr5)+(tau5/tau_cr5)**2)

  # %----Point 6----%

  tau6=-q/t2+tau3_T
  sig_eq6=sqrt(sig_6**2+3*tau6**2)
  G[39:42]=k*sig_eq6

  # %-----Constraints-----%

  Sig_C=65000*144
  Sig_T=65000*144


  G1 = np.zeros((72))

  G1[0:3]=((G[0:3])/Sig_C)-1
  G1[54:57]=-(G[0:3])/Sig_C-1
  G1[3:9]=G[3:9]-1

  G1[9:12]=(G[9:12])/Sig_C-1
  G1[57:60]=-(G[9:12])/Sig_C-1
  G1[12:18]=G[12:18]-1

  G1[18:21]=(G[18:21])/Sig_C-1
  G1[60:63]=-(G[18:21])/Sig_C-1
  G1[21:27]=G[21:27]-1

  G1[27:30]=(G[27:30])/Sig_T-1
  G1[63:66]=-(G[27:30])/Sig_T-1

  G1[30:33]=(G[30:33])/Sig_T-1
  G1[66:69]=-(G[30:33])/Sig_T-1
  G1[33:39]=G[33:39]-1

  G1[39:42]=(G[39:42])/Sig_T-1
  G1[69:72]=-(G[39:42])/Sig_T-1

  G1[42:45]=(1/2)*(ts1+ts3)/h-1
  G1[45:48]=t1/(ts1-.1*t1)-1
  G1[48:51]=t2/(ts2-.1*t2)-1
  G1[51:54]=t3/(ts3-.1*t3)-1

  constraints = [xx for xx in G1]

  return[0, constraints] # you have 72 constraints!!

def SBJ():
  """ Supersonic Business Jet aircraft conceptual design"""
  # Variables definition
  s  = COUPLING_TYPE.SHARED
  ff = COUPLING_TYPE.FEEDFORWARD
  fb = COUPLING_TYPE.FEEDBACK
  un = COUPLING_TYPE.UNCOUPLED
  dum = COUPLING_TYPE.DUMMY
  v = {
    "var1":     {"index": 1     , "sp_index": 1 , f"name": "SFC"        , "coupling_type": fb   , "link": 2     , "lb": 1     , "ub": 4       , "baseline": 1.0},   
    "var2":     {"index": 2     , "sp_index": 1 , f"name": "We"         , "coupling_type": fb   , "link": 2     , "lb": 100   , "ub": 30000   , "baseline": 1.0},   
    "var3":     {"index": 3     , "sp_index": 1 , f"name": "Wt"         , "coupling_type": ff   , "link": 3     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
    "var4":     {"index": 4     , "sp_index": 1 , f"name": "LD"         , "coupling_type": fb   , "link": 3     , "lb": 0.1   , "ub": 10      , "baseline": 1.0},   
    "var5":     {"index": 5     , "sp_index": 1 , f"name": "Ws"         , "coupling_type": fb   , "link": 4     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
    "var6":     {"index": 6     , "sp_index": 1 , f"name": "Wf"         , "coupling_type": fb   , "link": 4     , "lb": 5000  , "ub": 100000  , "baseline": 1.0}, 
    ## subproblem 2:
    "var7":     {"index": 7     , "sp_index": 2 , f"name": "SFC"        , "coupling_type": ff   , "link": 1     , "lb": 1     , "ub": 4       , "baseline": 1.0},
    "var8":     {"index": 8     , "sp_index": 2 , f"name": "We"         , "coupling_type": ff   , "link": 1     , "lb": 100   , "ub": 30000   , "baseline": 1.0},   
    "var9":     {"index": 9     , "sp_index": 2 , f"name": "D"          , "coupling_type": fb   , "link": 3     , "lb": 1000  , "ub": 70000   , "baseline": 1.0},   
    "var10":    {"index": 10    , "sp_index": 2 , f"name": "ESF"        , "coupling_type": ff   , "link": 3     , "lb": 0.5   , "ub": 1.5     , "baseline": 1.0}, 
    # local variable:
    "var11":    {"index": 11    , "sp_index": 2 , f"name": "T"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 1.0     , "baseline": 1.0},
    ## subproblem 3:
    "var12":    {"index": 12    , "sp_index": 3 , f"name": "D"          , "coupling_type": ff   , "link": 2     , "lb": 1000  , "ub": 70000   , "baseline": 1.0},
    "var13":    {"index": 13    , "sp_index": 3 , f"name": "ESF"        , "coupling_type": fb   , "link": 2     , "lb": 0.5   , "ub": 1.5     , "baseline": 1.0},   
    "var14":    {"index": 14    , "sp_index": 3 , f"name": "Wt"         , "coupling_type": fb   , "link": 1     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
    "var15":    {"index": 15    , "sp_index": 3 , f"name": "LD"         , "coupling_type": ff   , "link": 1     , "lb": 0.1   , "ub": 10      , "baseline": 1.0},   
    "var16":    {"index": 16    , "sp_index": 3 , f"name": "theta"      , "coupling_type": fb   , "link": 4     , "lb": 0.2   , "ub": 50      , "baseline": 1.0},   
    "var17":    {"index": 17    , "sp_index": 3 , f"name": "L"          , "coupling_type": ff   , "link": 4     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},  
    # shared variables:
    "var18":    {"index": 18    , "sp_index": 3 , f"name": "tc"         , "coupling_type": s    , "link": 4     , "lb": 0.01  , "ub": 0.1     , "baseline": 1.0},
    "var19":    {"index": 19    , "sp_index": 3 , f"name": "ARw"        , "coupling_type": s    , "link": 4     , "lb": 2.5   , "ub": 8.0     , "baseline": 1.0},   
    "var20":    {"index": 20    , "sp_index": 3 , f"name": "LAMBDAw"    , "coupling_type": s    , "link": 4     , "lb": 40.   , "ub": 70.     , "baseline": 1.0},   
    "var21":    {"index": 21    , "sp_index": 3 , f"name": "Sref"       , "coupling_type": s    , "link": 4     , "lb": 200.  , "ub": 800.    , "baseline": 1.0},   
    "var22":    {"index": 22    , "sp_index": 3 , f"name": "Sht"        , "coupling_type": s    , "link": 4     , "lb": 50    , "ub": 148.9   , "baseline": 1.0},   
    "var23":    {"index": 23    , "sp_index": 3 , f"name": "ARht"       , "coupling_type": s    , "link": 4     , "lb": 2.5   , "ub": 8.5     , "baseline": 1.0},  
    # local variables: 
    "var24":    {"index": 24    , "sp_index": 3 , f"name": "LAMBDAht"   , "coupling_type": un   , "link": None  , "lb": 40.   , "ub": 70.     , "baseline": 1.0},
    "var25":    {"index": 25    , "sp_index": 3 , f"name": "Lw"         , "coupling_type": un   , "link": None  , "lb": 0.01  , "ub": 0.2     , "baseline": 1.0},   
    "var26":    {"index": 26    , "sp_index": 3 , f"name": "Lht"        , "coupling_type": un   , "link": None  , "lb": 1     , "ub": 3.5     , "baseline": 1.0},   
    ## subproblem 4:
    "var27":    {"index": 27    , "sp_index": 4 , f"name": "theta"      , "coupling_type": ff   , "link": 3     , "lb": 0.2   , "ub": 50      , "baseline": 1.0},
    "var28":    {"index": 28    , "sp_index": 4 , f"name": "L"          , "coupling_type": fb   , "link": 3     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
    "var29":    {"index": 29    , "sp_index": 4 , f"name": "Ws"         , "coupling_type": ff   , "link": 1     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
    "var30":    {"index": 30    , "sp_index": 4 , f"name": "Wf"         , "coupling_type": ff   , "link": 1     , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
    # shared variables:
    "var31":    {"index": 31    , "sp_index": 4 , f"name": "tc"         , "coupling_type": s    , "link": 3     , "lb": 0.01  , "ub": 0.1     , "baseline": 1.0},
    "var32":    {"index": 32    , "sp_index": 4 , f"name": "ARw"        , "coupling_type": s    , "link": 3     , "lb": 2.5   , "ub": 8.0     , "baseline": 1.0},   
    "var33":    {"index": 33    , "sp_index": 4 , f"name": "LAMBDAw"    , "coupling_type": s    , "link": 3     , "lb": 40.   , "ub": 70.     , "baseline": 1.0},   
    "var34":    {"index": 34    , "sp_index": 4 , f"name": "Sref"       , "coupling_type": s    , "link": 3     , "lb": 200.  , "ub": 800.    , "baseline": 1.0},   
    "var35":    {"index": 35    , "sp_index": 4 , f"name": "Sht"        , "coupling_type": s    , "link": 3     , "lb": 50    , "ub": 148.9   , "baseline": 1.0},   
    "var36":    {"index": 36    , "sp_index": 4 , f"name": "ARht"       , "coupling_type": s    , "link": 3     , "lb": 2.5   , "ub": 8.5     , "baseline": 1.0},   
    # local variables:
    "var37":    {"index": 37    , "sp_index": 4 , f"name": "lambda"     , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 0.4     , "baseline": 1.0},
    "var38":    {"index": 38    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
    "var39":    {"index": 39    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
    "var40":    {"index": 40    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
    "var41":    {"index": 41    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
    "var42":    {"index": 42    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
    "var43":    {"index": 43    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
    "var44":    {"index": 44    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
    "var45":    {"index": 45    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
    "var46":    {"index": 46    , "sp_index": 4 , f"name": "t"          , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
    "var47":    {"index": 47    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var48":    {"index": 48    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var49":    {"index": 49    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var50":    {"index": 50    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var51":    {"index": 51    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var52":    {"index": 52    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var53":    {"index": 53    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var54":    {"index": 54    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    "var55":    {"index": 55    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un   , "link": None  , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
    # INTERMEDIATE PARAMETERS (NOT TO BE OPTIMIZED)
    "var56":  {"index": 56      , "sp_index": 1 ,  f"name": "range"     , "coupling_type": dum  ,  "link": None , "lb": 0.0   , "ub": 0.0     , "baseline": 1.0},
    "var57":  {"index": 57      , "sp_index": 2 ,  f"name": "Temp_E"    , "coupling_type": dum  ,  "link": None , "lb": 0.0   , "ub": 0.0     , "baseline": 1.0},
    "var58":  {"index": 58      , "sp_index": 2 ,  f"name": "Throttle_uA", "coupling_type": dum ,  "link": None , "lb": 0.0   , "ub": 0.0     , "baseline": 1.0},
    "var59":  {"index": 59      , "sp_index": 3 ,  f"name": "Pg"        , "coupling_type": dum  ,  "link": None , "lb": 0.0   , "ub": 0.0     , "baseline": 1.0},
    "var60":  {"index": 60      , "sp_index": 3 ,  f"name": "CLo1"      , "coupling_type": dum  ,  "link": None , "lb": 0.0   , "ub": 0.0     , "baseline": 1.0},
    "var61":  {"index": 61      , "sp_index": 3 ,  f"name": "CLo2"      , "coupling_type": dum  ,  "link": None , "lb": 0.0   , "ub": 0.0     , "baseline": 1.0},
  }

  Qscaling = []
  for key,value in v.items():
    scaling = v[key]["ub"] - v[key]["lb"]
    v[key]["scaling"] = scaling
    v[key]["dim"] = 1
    v[key]["value"] = v[key]["baseline"]
    Qscaling.append(.1/scaling if scaling != 0.0 and scaling != 0.0 else 1.)

  V: List[variableData] = []
  for i in range(len(v)):
    V.append(variableData(**v[f"var{i+1}"]))

  # Analyses setup; construct disciplinary analyses
  X1 = get_inputs(V,["SFC","We","LD","Ws","Wf"],sp_index=1)
  Y1 = get_outputs(V,["Wt","range"],sp_index=1)
  DA1: process = DA(inputs=X1,  
  outputs=Y1,
  blackbox=SBJ_A1,
  links=[2,3,4],
  coupling_type=COUPLING_TYPE.FEEDFORWARD)
  
  X2 = get_inputs(V,["D","T"],sp_index=2) # inputs
  Y2 = get_outputs(V,["SFC","We","ESF","Temp_E","Throttle_uA"],sp_index=2)
  DA2: process = DA(inputs=X2,
  outputs=Y2,
  blackbox=SBJ_A2,
  links=[1,3],
  coupling_type=[COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD]
  )

  X3 = get_inputs(V,["Wt","ESF","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],sp_index=3) # inputs
  Y3 = get_outputs(V,["L","D","LD","Pg","CLo1","CLo2"],sp_index=3)
  DA3: process = DA(inputs=X3,
  outputs=Y3,
  blackbox=SBJ_A3,
  links=[1,2,4],
  coupling_type=[COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD]
  ) 

  X4 = get_inputs(V,["L","tc","ARw","LAMBDAw","Sref","Sht","ARht","lambda","t","ts"],sp_index=4)
  Y4 = get_outputs(V,["Ws","Wf","theta"],sp_index=4)
  DA4: process = DA(inputs=X4,
  outputs=Y4,
  blackbox=SBJ_A4,
  links=[1,3],
  coupling_type=[COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD]
  )

  # MDA setup; construct subproblems MDA
  sp1_MDA: process = MDA(nAnalyses=1, analyses = [DA1], variables=X1, responses=Y1)
  sp2_MDA: process = MDA(nAnalyses=1, analyses = [DA2], variables=X2, responses=Y2)
  sp3_MDA: process = MDA(nAnalyses=1, analyses = [DA3], variables=X3, responses=Y3)
  sp4_MDA: process = MDA(nAnalyses=1, analyses = [DA4], variables=X4, responses=Y4)

  # Construct the coordinator
  coord = ADMM(beta = 1.3,
  nsp=4,
  budget = 500,
  index_of_master_SP=1,
  display = True,
  scaling = Qscaling,
  mode = "serial",
  M_update_scheme= w_scheme.MEDIAN,
  store_q_io=True
  )

  # Construct subproblems
  V1 = get_inputs(V,["SFC","We","LD","Ws","Wf"],sp_index=1) # inputs
  R1 = get_outputs(V,["Wt",],sp_index=1)

  sp1 = SubProblem(nv = len(V1),
  index = 1,
  vars = V1,
  resps = R1,
  is_main=1,
  analysis= sp1_MDA,
  coordination=coord,
  opt=SBJ_opt1,
  fmin_nop=np.inf,
  budget=20,
  display=False,
  psize = 1.,
  pupdate=PSIZE_UPDATE.LAST,
  freal=30000)

  V2 = get_inputs(V,["D","T"],sp_index=2) # inputs
  R2 = get_outputs(V,["SFC","We","ESF"],sp_index=2)

  sp2 = SubProblem(nv = len(V2),
  index = 2,
  vars = V2,
  resps = R2,
  is_main=0,
  analysis= sp2_MDA,
  coordination=coord,
  opt=SBJ_opt2,
  fmin_nop=np.inf,
  budget=20,
  display=False,
  psize = 1.,
  pupdate=PSIZE_UPDATE.LAST
  )

  V3 = get_inputs(V,["Wt","ESF","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],sp_index=3) # inputs
  R3 = get_outputs(V,["L","D","LD"],sp_index=3)

  sp3 = SubProblem(nv = len(V3),
  index = 3,
  vars = V3,
  resps = R3,
  is_main=0,
  analysis= sp3_MDA,
  coordination=coord,
  opt=SBJ_opt3,
  fmin_nop=np.inf,
  budget=20,
  display=False,
  psize = 1.,
  pupdate=PSIZE_UPDATE.LAST)

  V4 = get_inputs(V,["L","tc","ARw","LAMBDAw","Sref","Sht","ARht","lambda","t","ts"],sp_index=4) # inputs
  R4 = get_outputs(V,["Ws","Wf","theta"],sp_index=4)

  sp4 = SubProblem(nv = len(V4),
  index = 4,
  vars = V4,
  resps = R4,
  is_main=0,
  analysis= sp4_MDA,
  coordination=coord,
  opt=SBJ_opt4,
  fmin_nop=np.inf,
  budget=20,
  display=False,
  psize = 1.,
  pupdate=PSIZE_UPDATE.LAST)

  # Construct MDO workflow
  MDAO: MDO = MDO(
  Architecture = MDO_ARCHITECTURE.IDF,
  Coordinator = coord,
  subProblems = [sp1, sp2, sp3, sp4],
  variables = copy.deepcopy(V),
  responses = R1+R2+R3+R4,
  fmin = np.inf,
  hmin = np.inf,
  display = True,
  inc_stop = 1E-9,
  stop = "Iteration budget exhausted",
  tab_inc = [],
  noprogress_stop = 500
  )

  user = USER
  # Altitude
  user.h = 55000
  # Mach number
  user.M = 1.4

  return MDAO
  
if __name__ == "__main__":
  MDAO = SBJ()
  out = MDAO.run()
  print(out)
