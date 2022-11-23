import os
from DMDO import *
import math
import numpy as np
from numpy import cos, exp, pi, prod, sin, sqrt, subtract, inf
import copy
from typing import List, Dict, Any, Callable, Protocol, Optional

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
  
def SBJ_A1(x):
  # SBJ_obj_range
  inputs = {}
  for input_name,value in zip(["SFC","We","LD","Ws","Wf"],x): # these are your targets
    inputs[input_name] = value

  Wt = inputs["We"]+inputs["Wf"]+inputs["Ws"] # total weight
  # ["Wt"] # these are your responses
  y = Wt
  return y

def SBJ_opt1(x, y):
  inputs = {}
  for input_name,value in zip(["SFC","We","LD","Ws","Wf"],x): # these are your targets
    inputs[input_name] = value

  outputs = {}
  for outputs_name,value in zip(["Wt"],y): # these are your responses
    outputs[outputs_name] = value

  # SBJ_constraint_range
  if user.h <36089:
    theta_r = 1-0.000006875*user.h
  else:
    theta_r = 0.7519

  r = user.M * inputs["LD"] * 661 * np.sqrt(theta_r/inputs["SFC"])*math.log(outputs["Wt"]/(outputs["Wt"]-inputs["Wf"]))
  G = -r/2000 +1

  return [outputs["Wt"], [G]] # one objective, one constraint

def SBJ_A2(x):
  # SBJ_obj_power
  inputs = {}
  for input_name,value in zip(["D","T"],x): # these are your targets
    inputs[input_name] = value

  #  THIS SECTION COMPUTES SFC, ESF, AND ENGINE WEIGHT
  C=[500.0, 16000.0, 4.0 , 4360.0,  0.01375,  1.0]
  Thrust = inputs["D"]
  Dim_Throttle = inputs["T"]*16168
  s=[1.13238425638512, 1.53436586044561, -0.00003295564466, -0.00016378694115, -0.31623315541888, 0.00000410691343, -0.00005248000590, -0.00000000008574, 0.00000000190214, 0.00000001059951]
  SFC = s[0]+s[1]*user.M+s[2]*user.h+s[3]*Dim_Throttle+s[4]*user.M**2+2*user.h*user.M*s[5]+2*Dim_Throttle*user.M*s[6]+s[7]*user.h**2+2*Dim_Throttle*user.h*s[8]+s[9]*Dim_Throttle**2

  ESF = (Thrust/2)/Dim_Throttle

  We = C[3]*(ESF**1.05)*2

  # ["SFC","We","ESF"] # these are your responses
  y = [SFC,We,ESF]

def SBJ_opt2(x, y):
  # SBJ_constraint_power
  inputs = {}
  for input_name,value in zip(["D","T"],x): # these are your targets
    inputs[input_name] = value

  outputs = {}
  for outputs_name,value in zip(["SFC","We","ESF"],y): # these are your responses
    outputs[outputs_name] = value

  # -----THIS SECTION COMPUTES POLYNOMIAL CONSTRAINT FUNCTIONS-----
  Dim_Throttle = inputs["T"]*16168
  S_initial1=[user.M,user.h,inputs["D"]]
  S1=[user.M,user.h,inputs["T"]]
  flag1 = [2,4,2]
  bound1 = [.25,.25,.25]
  Temp_uA=1.02
  g1 = polyApprox(S_initial1, S1, flag1, bound1)
  g1 = g1 /Temp_uA-1
  p=[11483.7822254806, 10856.2163466548, -0.5080237941, 3200.157926969, -0.1466251679, 0.0000068572]
  Throttle_uA=p[0]+p[1]*user.M+p[2]*user.h+p[3]*user.M**2+2*p[4]*user.M*user.h+p[5]*user.h**2

  return[0, [g1, Dim_Throttle/Throttle_uA-1]] # one objective, two constraints

def SBJ_A3(x):
  # SBJ_obj_dragpolar
  inputs = {}
  for input_name,value in zip(["ESF","Wt","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],x): # these are your targets
    inputs[input_name] = value

  # %----Inputs----%
  Z = [inputs["tc"],user.h,user.M,inputs["ARw"],inputs["LAMBDAw"],inputs["Sref"],inputs["Sht"],inputs["ARht"]]
  C = [500.0, 16000.0,  4.0,  4360.0,  0.01375,  1.0]
  Z = [1,	55000,	1.40000000000000,	1,	1,	1,	1,	1]
  ARht=Z[7]
  S_ht=Z[6]    
  Nh=C[5]

  # %-----Drag computations----%

  if Z[1]<36089:
     V = Z[2]*(1116.39*sqrt(1-(6.875e-06*Z[1])))
     rho = (2.377e-03)*(1-(6.875e-06*Z[1]))**4.2561
  else:
     V = Z[2]*968.1
     rho = (2.377e-03)*(.2971)*np.exp(-(Z[1]-36089)/20806.7)
  q=.5*rho*(V**2)

  # ### Modified by S. Tosserams:
  # # scale coefficients for proper conditioning of matrix A 
  a=q*Z[5]/1e5
  b=Nh*q*S_ht/1e5
  # -------------------------------
  c= inputs["Lw"] #Lw
  d=(inputs["Lht"])*Nh*(S_ht/Z[5])
  A= np.array([[a, b], [c, d]])
  # ---- Modified by S. Tosserams:
  # ---- scale coefficient Wt for proper conditioning of matrix A 
  B=np.array([inputs["Wt"]/1e5, 0])
  # -----------------------
  try:
    CLo=np.linalg.solve(A, B)
  except:
    CLo = np.array([-np.inf, np.inf])
  delta_L=inputs["theta"]*q
  Lw1=CLo[0]*q*Z[5]-delta_L
  CLw1=Lw1/(q*Z[5])
  CLht1=-CLw1*c/d
  #  Modified by S. Tosserams:
  #  scale first coefficient of D for proper conditioning of matrix A 
  D=np.array([(inputs["Wt"]-CLw1*a-CLht1*b)/1e5, -CLw1*c-CLht1*d])
  # -----------------
  try:
    DCL = np.linalg.solve(A, D)
  except:
    DCL = np.array([np.nan,np.nan])

  if Z[2] >= 1:
    kw = Z[3] * (Z[2]**2-1) * np.cos(Z[4]*np.pi/180)/(4*Z[3]*np.sqrt(Z[2]**2-1)-2)
    kht = ARht * (Z[2]**2-1)*np.cos(inputs["LAMBDAht"]*pi/180)/(4*ARht*np.sqrt(Z[2]**2-1)-2)
  else:
    kw = 1/(np.pi*0.8*Z[3])
    kht= 1/(np.pi*0.8*ARht)
  
  S_initial1 = copy.deepcopy(inputs["ESF"])
  S1 = copy.deepcopy(inputs["ESF"])
  flag1 = 1
  bound1 = 0.25
  Fo1 = polyApprox(S_initial1 if isinstance(S_initial1, list) else [S_initial1], S1 if isinstance(S1, list) else [S1], flag1 if isinstance(flag1, list) else [flag1], bound1 if isinstance(bound1, list) else [bound1])

  CDmin = C[4]*Fo1 + 3.05*(Z[0]**(5/3))*((cos(Z[4]*np.pi/180))**(3/2))

  CDw=CDmin+kw*(CLo[0]**2)+kw*(DCL[0]**2)
  CDht=kht*(CLo[1]**2)+kht*(DCL[1]**2)
  CD=CDw+CDht
  CL=CLo[0]+CLo[1]
  L = inputs["Wt"]
  D = q*CDw*Z[5]+q*CDht*Z[6]
  LD = CL/CD

  # ["L","D","LD"] # these are your responses

  return [L, D, LD]

def SBJ_opt3(x, y):
  inputs = {}
  for input_name,value in zip(["ESF","Wt","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],x): # these are your targets
    inputs[input_name] = value

  # SBJ_constraint_dragpolar
  # %----Inputs----%
  Z = [inputs["tc"],user.h,user.M,inputs["ARw"],inputs["LAMBDAw"],inputs["Sref"],inputs["Sht"],inputs["ARht"]]
  C = [500.0, 16000.0,  4.0,  4360.0,  0.01375,  1.0]
  S_ht=Z[6]
  Nh=C[5]

  # %-----Drag computations----%

  if Z[1]<36089:
     V = Z[2]*(1116.39*sqrt(1-(6.875e-06*Z[1])))
     rho = (2.377e-03)*(1-(6.875e-06*Z[1]))**4.2561
  else:
     V = Z[2]*968.1
     rho = (2.377e-03)*(.2971)*exp(-(Z[1]-36089)/20806.7)
  q=.5*rho*(V**2)

  # ### Modified by S. Tosserams:
  # # scale coefficients for proper conditioning of matrix A 
  a=q*Z[5]/1e5
  b=Nh*q*S_ht/1e5
  # -------------------------------
  c= inputs["Lw"] #Lw
  d=(inputs["Lht"])*Nh*(S_ht/Z[5])
  A= np.array([[a, b], [c, d]])
  # ---- Modified by S. Tosserams:
  # ---- scale coefficient Wt for proper conditioning of matrix A 
  B=[inputs["Wt"]/1e5, 0]
  # -----------------------
  try:
    CLo=np.linalg.solve(A, B)
  except:
    CLo = np.array([-np.inf, np.inf])
  
  S_initial2 = copy.deepcopy(inputs["tc"])
  S2 = copy.deepcopy(Z[0])
  flag1 = [1]
  bound1 = [0.25]
  
  g1 = polyApprox(S_initial2 if isinstance(S_initial2, list) else [S_initial2], S2 if isinstance(S2, list) else [S2], flag1 if isinstance(flag1, list) else [flag1], bound1 if isinstance(bound1, list) else [bound1])


  # %------Constraints------%

  Pg_uA=1.1
  g1=g1/Pg_uA-1
  if CLo[0] > 0:
      g2=(2*(CLo[1]))-(CLo[0])
      g3=(2*(-CLo[1]))-(CLo[0])
  else:
      g2=(2*(-CLo[1]))-(CLo[0])
      g3=(2*(CLo[1]))-(CLo[0])

  return [0, [g1, g2, g3]] # no objective, three constraints

def SBJ_A4(x):
  # SBJ_obj_weight
  inputs = {}
  for input_name,value in zip(["L","tc","ARw","LAMBDAw","Sref","Sht","ARht","lambda"],x): # these are your targets
    inputs[input_name] = value

  index_t = len(inputs) # get second last 9 entries
  index_ts = len(inputs)+9 # get last 9 entries

  Z = [inputs["tc"],user.h,user.M,inputs["ARw"],inputs["LAMBDAw"],inputs["Sref"],inputs["Sht"],inputs["ARht"]]
  L = inputs["L"]
  # %----Inputs----%
  t= np.divide(x[index_t:index_t+9], 12.)#convert to feet
  ts= np.divide(x[index_ts:], 12.);#convert to feet
  LAMBDA = inputs["lambda"]
  C=[500.0, 16000.0,  4.0,  4360.0,  0.01375,  1.0]

  t1= [t[i] for i in range(3)] 
  t2= [t[i] for i in range(3, 6)]  
  t3= [t[i] for i in range(6, 9)]  
  ts1=[ts[i] for i in range(3)] 
  ts2=[ts[i] for i in range(3, 6)]  
  ts3=[ts[i] for i in range(6, 9)]  
  G=4000000*144
  E=10600000*144
  rho_alum=0.1*144
  rho_core=0.1*144/10
  rho_fuel=6.5*7.4805
  Fw_at_t=5

  # %----Inputs----%
  beta=.9
  c,c_box,Sweep_40,D_mx,b,_ = Wing_Mod(Z,LAMBDA)
  l=0.6*c_box
  h=(np.multiply([c[i] for i in range(3)], beta*float(Z[0])))-np.multiply((0.5),np.add(ts1,ts3))
  A_top=(np.multiply(t1, 0.5*l))+(np.multiply(t2, h/6))
  A_bottom=(np.multiply(t3, 0.5*l))+(np.multiply(t2, h/6))
  Y_bar=np.multiply(h, np.divide((2*A_top), (2*A_top+2*A_bottom)))
  Izz=np.multiply(2, np.multiply(A_top, np.power((h-Y_bar), 2)))+np.multiply(2, np.multiply(A_bottom,np.power((-Y_bar), 2)))
  P,Mz,Mx,bend_twist,Spanel=loads(b,c,Sweep_40,D_mx,L,Izz,Z,E)

  Phi=(Mx/(4*G*(l*h)**2))*(l/t1+2*h/t2+l/t3)
  aa=len(bend_twist)
  twist = np.array([0] * aa)
  twist[0:int(aa/3)]=bend_twist[0:int(aa/3)]+Phi[0]*180/pi
  twist[int(aa/3):int(aa*2/3)]=bend_twist[int(aa/3):int(aa*2/3)]+Phi[1]*180/pi
  twist[int(aa*2/3):aa]=bend_twist[int(aa*2/3):aa]+Phi[2]*180/pi
  deltaL_divby_q=sum(twist*Spanel*0.1*2)


  # %-----THIS SECTION COMPUTES THE TOTAL WEIGHT OF A/C-----%

  Wtop_alum=(b/4)*(c[0]+c[3])*np.mean(t1)*rho_alum
  Wbottom_alum=(b/4)*(c[0]+c[3])*np.mean(t3)*rho_alum
  Wside_alum=(b/2)*np.mean(h)*np.mean(t2)*rho_alum
  Wtop_core=(b/4)*(c[0]+c[3])*np.mean(np.subtract(ts1, t1))*rho_core
  Wbottom_core=(b/4)*(c[0]+c[3])*np.mean(np.subtract(ts3,t3))*rho_core
  Wside_core=(b/2)*np.mean(h)*np.mean(np.subtract(ts2,t2))*rho_core
  W_wingstruct=Wtop_alum+Wbottom_alum+Wside_alum+Wtop_core+Wbottom_core+Wside_core
  W_fuel_wing=np.mean(h*l)*(b/3)*(2)*rho_fuel
  Bh=sqrt(Z[7]*Z[6])
  W_ht=3.316*((1+(Fw_at_t/Bh))**-2.0)*((L*C[2]/1000)**0.260)*(Z[6]**0.806)
  Wf = C[0] + W_fuel_wing
  Ws = C[1] + W_ht + 2*W_wingstruct
  theta = deltaL_divby_q
  
  return [Ws, Wf, theta] # these are your responses

def Wing_Mod(Z, LAMBDA):
  c = [0, 0, 0, 0]
  x = [0]*8
  y = [0] * 8
  b=max(2,np.real(sqrt(Z[3]*Z[5])))
  c[0]=2*Z[5]/((1+LAMBDA)*b)
  c[3]=LAMBDA*c[0]
  x[0]=0
  y[0]=0
  x[1]=c[0]
  y[1]=0
  x[6]=(b/2)*np.tan(Z[4]*np.pi/180)
  y[6]=b/2
  x[7]=x[6]+c[3]
  y[7]=b/2
  y[2]=b/6
  x[2]=(x[6]/y[6])*y[2]
  y[4]=b/3
  x[4]=(x[6]/y[6])*y[4]
  x[5]=x[7]+((x[1]-x[7])/y[7])*(y[7]-y[4])
  y[5]=y[4]
  x[3]=x[7]+((x[1]-x[7])/y[7])*(y[7]-y[2])
  y[3]=y[2]
  c[1]=x[3]-x[2]
  c[2]=x[5]-x[4]
  TE_sweep=(np.arctan((x[7]-x[1])/y[7]))*180/np.pi
  Sweep_40=(np.arctan(((x[7]-0.6*(x[7]-x[6]))-0.4*x[1])/y[7]))*180/np.pi

  l=np.multiply([c[i] for i in range(3)], .4*np.cos(Z[4]*np.pi/180))
  k=np.multiply([c[i] for i in range(3)], .6*np.sin((90-TE_sweep)*np.pi/180)/sin((90+TE_sweep-Z[4])*np.pi/180))
  c_box=np.add(l, k)
  D_mx=np.subtract(l, np.multiply(0.407, c_box))


  return c,c_box,Sweep_40,D_mx,b,l

def loads(b,c,Sweep_40,D_mx,L,Izz,Z,E):
  NP=9 #   %----number of panels per halfspan
  n=90
  rn = int(n/NP)

  h=(b/2)/n
  x=np.linspace(0,b/2-h, n)
  x1=np.linspace(h, b/2, n)

  #%----Calculate Mx, Mz, and P----%

  l=np.linspace(0, (b/2)-(b/2)/NP, NP)


  c1mc4 = c[0]-c[3]
  f_all  =np.multiply((3*b/10), np.sqrt(np.subtract(1, ( np.power(x, 2))/(np.power(np.divide(b,2),2)))))
  f1_all =np.multiply((3*b/10), np.sqrt(np.subtract(1, ( np.power(x1, 2))/(np.power(np.divide(b,2),2)))))
  C= c[3] + 2*( (b/2- x)/b )*c1mc4
  C1=c[3] + 2*( (b/2-x1)/b )*c1mc4
  A_Tot: np.ndarray =np.multiply((h/4)*(C+C1), (np.add(f_all, f1_all)))
  Area = np.sum(A_Tot.reshape((NP,rn)), axis=1)
  Spanel = np.multiply((h*(rn)/2), (np.add([C[int(i)] for i in np.linspace(0, n-10, 9)], [C[int(i)] for i in np.linspace(9, n-1, 9)])))

    # % cos, tan, and cos**-1 of Sweep
  cosSweep = np.cos(Sweep_40*pi/180)
  cosInvSweep = 1/cosSweep
  tanCos2Sweep = np.tan(Sweep_40*pi/180)*cosSweep*cosSweep



  p=np.divide(L*Area, sum(Area))
  


  # % Replace T by:
  Tcsp = np.cumsum(p)
  Tsp = Tcsp[-1]
  temp = [0]
  temp = temp +[Tcsp[i] for i in range(len(Tcsp)-1)]
  T = np.subtract(Tsp, temp)
  pl = np.multiply(p, l)
  Tcspl = np.cumsum(pl)
  Tspl = Tcspl[-1]
  Mb = np.multiply(np.subtract( np.subtract(Tspl, Tcspl), np.multiply(l, np.subtract(Tsp, Tcsp)) ),cosInvSweep)


  P=[T[int(i)] for i in np.arange(0,NP-1, int(NP/3))]
  Mx=np.multiply(P, D_mx)
  Mz=[Mb[int(i)] for i in np.arange(0,NP-1, int(NP/3))]

  # %----Calculate Wing Twist due to Bending----%
  I = np.zeros((NP))
  chord=c[3]+ (np.divide(2*(b/2-l), b))*c1mc4
  y = np.zeros((2,9))
  y[0,:]=(l-.4*chord*tanCos2Sweep)*cosInvSweep
  y[1,:]=(l+.6*chord*tanCos2Sweep)*cosInvSweep
  y[1,0]=0
  I[0:int(NP/3)]=sqrt((Izz[0]**2+Izz[1]**2)/2)
  I[int(NP/3):int(2*NP/3)]=sqrt((Izz[1]**2+Izz[2]**2)/2)
  I[int(2*NP/3):int(NP)]=sqrt((Izz[2]**2)/2)

  La=y[0,1:NP]-y[0,0:NP-1]
  La = np.append(0, La)
  Lb=y[1,1:NP]-y[1,0:NP-1]
  Lb=np.append(0, Lb)
  A=T*La**3/(3*E*I)+Mb*La**2./(2*E*I)
  B=T*Lb**3/(3*E*I)+Mb*Lb**2./(2*E*I)
  Slope_A=T*La**2./(2*E*I)+Mb*La/(E*I)
  Slope_B=T*Lb**2./(2*E*I)+Mb*Lb/(E*I)
  for i in range(NP-1):
    Slope_A[i+1]=Slope_A[i]+Slope_A[i+1]
    Slope_B[i+1]=Slope_B[i]+Slope_B[i+1]
    A[i+1]=A[i]+Slope_A[i]*La[i+1]+A[i+1]
    B[i+1]=B[i]+Slope_B[i]*Lb[i+1]+B[i+1]

  bend_twist=((B-A)/chord)*180/pi
  for i in range(1, len(bend_twist)):
    if bend_twist[i]<bend_twist[i-1]:
        bend_twist[i]=bend_twist[i-1]

  return P,Mz,Mx,bend_twist,Spanel

def SBJ_opt4(x, y):
  inputs = {}
  for input_name,value in zip(["L","tc","ARw","LAMBDAw","Sref","Sht","ARht","lambda"],x): # these are your targets
    inputs[input_name] = value
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


  return[0, [xx for xx in G1]] # you have 72 constraints!!

def polyApprox(S, S_new, flag, S_bound):
  S_norm = []
  S_shifted = []
  Ai = []
  Aij = np.zeros((len(S),len(S)))
  for i in range(len(S)):
    S_norm.append(S_new[int(i/S[i])])
    if S_norm[i]>1.25:
      S_norm[i]=1.25
    elif S_norm[i]<0.75:
        S_norm[i]=0.75
    S_shifted.append(S_norm[i] - 1)
    a = 0.1
    b = a


    if flag[i]==5:
      # CALCULATE POLYNOMIAL COEFFICIENTS (S-ABOUT ORIGIN)
      So=0
      Sl=So-S_bound[i]
      Su=So+S_bound[i]
      Mtx_shifted = np.array([[1, Sl, Sl**2], [1, So, So**2], [1, Su, Su**2]])

      F_bound = np.array([1+(.5*a)**2, 1, 1+(.5*b)**2])
      A = np.linalg.solve(Mtx_shifted, F_bound)
      Ao = A[0]
      Ai.append(A[1])
      Aij[i,i] = A[2]

      # CALCULATE POLYNOMIAL COEFFICIENTS
    else:
      if flag[i] == 0:
        S_shifted.append(0)
      elif flag[i]==3:
        a *= -1.
        b=copy.deepcopy(a)
      elif flag[i]==2:
        b = 2 * a
      elif flag[i] == 4:
        a *= -1
        b = 2*a
      # DETERMINE BOUNDS ON FF DEPENDING ON SLOPE-SHAPE
      #  CALCULATE POLYNOMIAL COEFFICIENTS (S-ABOUT ORIGIN)
      So=0
      Sl=So-S_bound[i]
      Su=So+S_bound[i]
      Mtx_shifted = np.array([[1, Sl, Sl**2], [1, So, So**2], [1, Su, Su**2]])
      F_bound = np.array([1-.5*a, 1, 1+.5*b])
      A = np.linalg.solve(Mtx_shifted, F_bound)
      Ao = A[0]
      Ai.append(A[1])
      Aij[i,i] = A[2]
    
      #  CALCULATE POLYNOMIAL COEFFICIENTS
  R = np.array([[0.2736,    0.3970,    0.8152,    0.9230,    0.1108], 
                [0.4252,    0.4415,    0.6357,    0.7435,    0.1138],
                [0.0329,    0.8856,    0.8390,    0.3657,    0.0019],
                [0.0878,    0.7248,    0.1978,    0.0200,    0.0169],
                [0.8955,    0.4568,    0.8075,    0.9239,    0.2525]])
  
  for i in range(len(S)):
    for j in range(i+1, len(S)):
      Aij[i, j] = Aij[i,i] * R[i,j]
      Aij[j, i] = Aij[i, j]

  S_shifted = np.array(S_shifted)
  
  FF = Ao + np.dot(Ai, (np.transpose(S_shifted))) + (1/2)*np.dot(np.dot(S_shifted, Aij), np.transpose(S_shifted))

  return FF

def SBJ():
  """ Supersonic Business Jet aircraft conceptual design"""
  # Variables definition
  s  = COUPLING_TYPE.SHARED
  ff = COUPLING_TYPE.FEEDFORWARD
  fb = COUPLING_TYPE.FEEDBACK
  un = COUPLING_TYPE.UNCOUPLED
  dum = COUPLING_TYPE.DUMMY
  v = {
      "var1":     {"index": 1     , "sp_index": 1 , f"name": "SFC"        , "coupling_type": fb  , "link": 2       , "lb": 1     , "ub": 4       , "baseline": 1.0},   
      "var2":     {"index": 2     , "sp_index": 1 , f"name": "We"         , "coupling_type": fb  , "link": 2       , "lb": 100   , "ub": 30000   , "baseline": 1.0},   
      "var3":     {"index": 3     , "sp_index": 1 , f"name": "Wt"         , "coupling_type": ff  , "link": 3       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
      "var4":     {"index": 4     , "sp_index": 1 , f"name": "LD"         , "coupling_type": fb  , "link": 3       , "lb": 0.1   , "ub": 10      , "baseline": 1.0},   
      "var5":     {"index": 5     , "sp_index": 1 , f"name": "Ws"         , "coupling_type": fb  , "link": 4       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
      "var6":     {"index": 6     , "sp_index": 1 , f"name": "Wf"         , "coupling_type": fb  , "link": 4       , "lb": 5000  , "ub": 100000  , "baseline": 1.0}, 
      ## subproblem 2:
      "var7":     {"index": 7     , "sp_index": 2 , f"name": "SFC"        , "coupling_type": ff  , "link": 1       , "lb": 1     , "ub": 4       , "baseline": 1.0},
      "var8":     {"index": 8     , "sp_index": 2 , f"name": "We"         , "coupling_type": ff  , "link": 1       , "lb": 100   , "ub": 30000   , "baseline": 1.0},   
      "var9":     {"index": 9     , "sp_index": 2 , f"name": "D"          , "coupling_type": fb  , "link": 3       , "lb": 1000  , "ub": 70000   , "baseline": 1.0},   
      "var10":    {"index": 10    , "sp_index": 2 , f"name": "ESF"        , "coupling_type": ff  , "link": 3       , "lb": 0.5   , "ub": 1.5     , "baseline": 1.0}, 
      # local variable:
      "var11":    {"index": 11    , "sp_index": 2 , f"name": "T"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 1.0     , "baseline": 1.0},
      ## subproblem 3:
      "var12":    {"index": 12    , "sp_index": 3 , f"name": "D"          , "coupling_type": ff  , "link": 2       , "lb": 1000  , "ub": 70000   , "baseline": 1.0},
      "var13":    {"index": 13    , "sp_index": 3 , f"name": "ESF"        , "coupling_type": fb  , "link": 2       , "lb": 0.5   , "ub": 1.5     , "baseline": 1.0},   
      "var14":    {"index": 14    , "sp_index": 3 , f"name": "Wt"         , "coupling_type": fb  , "link": 1       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
      "var15":    {"index": 15    , "sp_index": 3 , f"name": "LD"         , "coupling_type": ff  , "link": 1       , "lb": 0.1   , "ub": 10      , "baseline": 1.0},   
      "var16":    {"index": 16    , "sp_index": 3 , f"name": "theta"      , "coupling_type": fb  , "link": 4       , "lb": 0.2   , "ub": 50      , "baseline": 1.0},   
      "var17":    {"index": 17    , "sp_index": 3 , f"name": "L"          , "coupling_type": ff  , "link": 4       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},  
      # shared variables:
      "var18":    {"index": 18    , "sp_index": 3 , f"name": "tc"         , "coupling_type": s   , "link": 4       , "lb": 0.01  , "ub": 0.1     , "baseline": 1.0},
      "var19":    {"index": 19    , "sp_index": 3 , f"name": "ARw"        , "coupling_type": s   , "link": 4       , "lb": 2.5   , "ub": 8.0     , "baseline": 1.0},   
      "var20":    {"index": 20    , "sp_index": 3 , f"name": "LAMBDAw"    , "coupling_type": s   , "link": 4       , "lb": 40.   , "ub": 70.     , "baseline": 1.0},   
      "var21":    {"index": 21    , "sp_index": 3 , f"name": "Sref"       , "coupling_type": s   , "link": 4       , "lb": 200.  , "ub": 800.    , "baseline": 1.0},   
      "var22":    {"index": 22    , "sp_index": 3 , f"name": "Sht"        , "coupling_type": s   , "link": 4       , "lb": 50    , "ub": 148.9   , "baseline": 1.0},   
      "var23":    {"index": 23    , "sp_index": 3 , f"name": "ARht"       , "coupling_type": s   , "link": 4       , "lb": 2.5   , "ub": 8.5     , "baseline": 1.0},  
      # local variables: 
      "var24":    {"index": 24    , "sp_index": 3 , f"name": "LAMBDAht"   , "coupling_type": un  , "link": None    , "lb": 40.   , "ub": 70.     , "baseline": 1.0},
      "var25":    {"index": 25    , "sp_index": 3 , f"name": "Lw"         , "coupling_type": un  , "link": None    , "lb": 0.01  , "ub": 0.2     , "baseline": 1.0},   
      "var26":    {"index": 26    , "sp_index": 3 , f"name": "Lht"        , "coupling_type": un  , "link": None    , "lb": 1     , "ub": 3.5     , "baseline": 1.0},   
      ## subproblem 4:
      "var27":    {"index": 27    , "sp_index": 4 , f"name": "theta"      , "coupling_type": ff  , "link": 3       , "lb": 0.2   , "ub": 50      , "baseline": 1.0},
      "var28":    {"index": 28    , "sp_index": 4 , f"name": "L"          , "coupling_type": fb  , "link": 3       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
      "var29":    {"index": 29    , "sp_index": 4 , f"name": "Ws"         , "coupling_type": ff  , "link": 1       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
      "var30":    {"index": 30    , "sp_index": 4 , f"name": "Wf"         , "coupling_type": ff  , "link": 1       , "lb": 5000  , "ub": 100000  , "baseline": 1.0},   
      # shared variables:
      "var31":    {"index": 31    , "sp_index": 4 , f"name": "tc"         , "coupling_type": s   , "link": 3       , "lb": 0.01  , "ub": 0.1     , "baseline": 1.0},
      "var32":    {"index": 32    , "sp_index": 4 , f"name": "ARw"        , "coupling_type": s   , "link": 3       , "lb": 2.5   , "ub": 8.0     , "baseline": 1.0},   
      "var33":    {"index": 33    , "sp_index": 4 , f"name": "LAMBDAw"    , "coupling_type": s   , "link": 3       , "lb": 40.   , "ub": 70.     , "baseline": 1.0},   
      "var34":    {"index": 34    , "sp_index": 4 , f"name": "Sref"       , "coupling_type": s   , "link": 3       , "lb": 200.  , "ub": 800.    , "baseline": 1.0},   
      "var35":    {"index": 35    , "sp_index": 4 , f"name": "Sht"        , "coupling_type": s   , "link": 3       , "lb": 50    , "ub": 148.9   , "baseline": 1.0},   
      "var36":    {"index": 36    , "sp_index": 4 , f"name": "ARht"       , "coupling_type": s   , "link": 3       , "lb": 2.5   , "ub": 8.5     , "baseline": 1.0},   
      # local variables:
      "var37":    {"index": 37    , "sp_index": 4 , f"name": "lambda"     , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 0.4     , "baseline": 1.0},
      "var38":    {"index": 38    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
      "var39":    {"index": 39    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
      "var40":    {"index": 40    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
      "var41":    {"index": 41    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},   
      "var42":    {"index": 42    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
      "var43":    {"index": 43    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
      "var44":    {"index": 44    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
      "var45":    {"index": 45    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
      "var46":    {"index": 46    , "sp_index": 4 , f"name": "t"          , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 4.0     , "baseline": 1.0},
      "var47":    {"index": 47    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var48":    {"index": 48    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var49":    {"index": 49    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var50":    {"index": 50    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var51":    {"index": 51    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var52":    {"index": 52    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var53":    {"index": 53    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var54":    {"index": 54    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
      "var55":    {"index": 55    , "sp_index": 4 , f"name": "ts"         , "coupling_type": un  , "link": None    , "lb": 0.1   , "ub": 9.0     , "baseline": 1.0},
  }

  fstar = 33600
  frealistic = 30000
  altitude = 55000
  Mach = 1.4

  Qscaling = []
  for key,value in v.items():
    scaling = v[key]["ub"] - v[key]["lb"]
    v[key]["scaling"] = scaling
    v[key]["dim"] = 1
    v[key]["value"] = v[key]["baseline"]
    Qscaling.append(.1/scaling if .1/scaling != np.inf and .1/scaling != np.nan else 1.)

  V: List[variableData] = []
  for i in range(len(v)):
    V.append(variableData(**v[f"var{i+1}"]))

  # Analyses setup; construct disciplinary analyses
  v1 = get_inputs(V,["SFC","We","LD","Ws","Wf"],sp_index=1) # inputs
  y1 = get_outputs(V,["Wt"],sp_index=1) # outputs
  V1 = copy.deepcopy(v1)
  Y1 = copy.deepcopy(y1)
  DA1: process = DA(inputs=V1,
  outputs=Y1,
  blackbox=SBJ_A1,
  links=[2,3,4],
  coupling_type=COUPLING_TYPE.FEEDFORWARD)
  
  v2 = get_inputs(V,["D","T"],sp_index=2) # inputs
  y2 = get_outputs(V,["SFC","We","ESF"],sp_index=2) # outputs
  V2 = copy.deepcopy(v2)
  Y2 = copy.deepcopy(y2)
  DA2: process = DA(inputs=V2,
  outputs=Y2,
  blackbox=SBJ_A2,
  links=[1,3],
  coupling_type=[COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD]
  )

  v3 = get_inputs(V,["ESF","Wt","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],sp_index=3) # inputs
  y3 = get_outputs(V,["L","D","LD"],sp_index=3) # outputs
  V3 = copy.deepcopy(v3)
  Y3 = copy.deepcopy(y3)
  DA3: process = DA(inputs=V3,
  outputs=Y3,
  blackbox=SBJ_A3,
  links=[1,2,4],
  coupling_type=[COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD]
  )

  v4 = get_inputs(V,["L","tc","ARw","LAMBDAw","Sref","Sht","ARht","lambda","t","ts"],sp_index=4) # inputs
  y4 = get_outputs(V,["Ws","Wf","theta"],sp_index=4) # outputs
  V4 = copy.deepcopy(v4)
  Y4 = copy.deepcopy(y4)
  DA4: process = DA(inputs=V4,
  outputs=Y4,
  blackbox=SBJ_A4,
  links=[1,3],
  coupling_type=[COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD, COUPLING_TYPE.FEEDFORWARD]
  )

# MDA setup; construct subproblems MDA
  sp1_MDA: process = MDA(nAnalyses=1, analyses = [DA1], variables=V1, responses=Y1)
  sp2_MDA: process = MDA(nAnalyses=1, analyses = [DA2], variables=V2, responses=Y2)
  sp3_MDA: process = MDA(nAnalyses=1, analyses = [DA3], variables=V3, responses=Y3)
  sp4_MDA: process = MDA(nAnalyses=1, analyses = [DA4], variables=V4, responses=Y4)

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
  sp1 = SubProblem(nv = len(V1),
  index = 1,
  vars = V1,
  resps = Y1,
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

  sp2 = SubProblem(nv = len(V2),
  index = 2,
  vars = V2,
  resps = Y2,
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

  sp3 = SubProblem(nv = len(V3),
  index = 3,
  vars = V3,
  resps = Y3,
  is_main=0,
  analysis= sp3_MDA,
  coordination=coord,
  opt=SBJ_opt3,
  fmin_nop=np.inf,
  budget=20,
  display=False,
  psize = 1.,
  pupdate=PSIZE_UPDATE.LAST)

  sp4 = SubProblem(nv = len(V4),
  index = 4,
  vars = V4,
  resps = Y4,
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
  variables = V,
  responses = Y1+Y2+Y3+Y4,
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
  
  out = MDAO.run()
  print(out)
  
if __name__ == "__main__":
  SBJ()
