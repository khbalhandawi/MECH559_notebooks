import os
import math
import numpy as np
import copy
from numpy import cos, exp, pi, prod, sin, sqrt, subtract, inf
from DMDO import USER

def SBJ_A1(x):
  inputs = {}
  for input_name,value in zip(["SFC","We","LD","Ws","Wf"],x): # these are your targets
    inputs[input_name] = value

  Wt = inputs["We"]+inputs["Wf"]+inputs["Ws"] # total weight

  # SBJ_constraint_range
  if user.h <36089:
    theta_r = 1-0.000006875*user.h
  else:
    theta_r = 0.7519

  range_aircraft = user.M * inputs["LD"] * 661 * np.sqrt(theta_r/inputs["SFC"])*math.log(Wt/(Wt-inputs["Wf"]))

  return [Wt, range_aircraft] # these are your responses

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

  S_initial1=[user.M,user.h,inputs["D"]]
  S1=[user.M,user.h,inputs["T"]]
  flag1 = [2,4,2]
  bound1 = [.25,.25,.25]
  Temp_E = polyApprox(S_initial1, S1, flag1, bound1)
  p=[11483.7822254806, 10856.2163466548, -0.5080237941, 3200.157926969, -0.1466251679, 0.0000068572]
  Throttle_uA=p[0]+p[1]*user.M+p[2]*user.h+p[3]*user.M**2+2*p[4]*user.M*user.h+p[5]*user.h**2
  
  return [SFC,We,ESF,Temp_E,Throttle_uA]

def SBJ_A3(x):
  # SBJ_obj_dragpolar
  inputs = {}
  for input_name,value in zip(["Wt","ESF","theta","tc","ARw","LAMBDAw","Sref","Sht","ARht","LAMBDAht","Lw","Lht"],x): # these are your targets
    inputs[input_name] = value

  # %----Inputs----%
  Z = [inputs["tc"],user.h,user.M,inputs["ARw"],inputs["LAMBDAw"],inputs["Sref"],inputs["Sht"],inputs["ARht"]]
  C = [500.0, 16000.0,  4.0,  4360.0,  0.01375,  1.0]
  # Z = [1,	55000,	1.40000000000000,	1,	1,	1,	1,	1]
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

  # %------Constraints------%
  # SBJ_constraint_dragpolar
  S_initial2 = copy.deepcopy(inputs["tc"])
  S2 = copy.deepcopy(Z[0])
  flag1 = [1]
  bound1 = [0.25]
  
  Pg = polyApprox(S_initial2 if isinstance(S_initial2, list) else [S_initial2], S2 if isinstance(S2, list) else [S2], flag1 if isinstance(flag1, list) else [flag1], bound1 if isinstance(bound1, list) else [bound1])

  return [L, D, LD, Pg, CLo[0], CLo[1]]

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

# %% Constants
user = USER
# Altitude
user.h = 55000
# Mach number
user.M = 1.4