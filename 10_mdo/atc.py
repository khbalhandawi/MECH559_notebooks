from typing import List, Dict, Any, Union
import OMADS
import numpy as np
import json
from copy import deepcopy

def get_by_address(address):
    return [x for x in globals().values() if id(x)==address]

class ATC():

    def __init__(self,disciplines:List[callable],P_global:callable,P_local:List[callable],
        variables:List[Dict[str,Union[str,float,int]]],gamma:float=0.1,beta:float=1.3,
        w0:float=1.0,v0:float=0.0,opt_options:Dict[str,Any]=None):

        self.disciplines = disciplines
        self.P_global = P_global
        assert len(P_local) == len(disciplines) # make sure there are as many disciplines as subproblems
        self.P_local = P_local
        self.n_disciplines = len(disciplines)
        self._variables = variables

        # make shared targets
        for i in range(self.n_disciplines):
            x0t_dicts = deepcopy(self.get_x("shared"))
            for d in x0t_dicts:
                d["target"] = True
                d["subproblem"] = i+1
                self._variables += [d]

        # make coupling targets
        for i in range(self.n_disciplines):
            ut_dicts = deepcopy(self.get_x("coupling"))
            for d in ut_dicts:
                d["target"] = True
                d["subproblem"] = i+1
                self._variables += [d]

        # set initial weights
        self.set_target(w0,"shared","w",subproblem=[1,2])
        self.set_target(v0,"shared","v",subproblem=[1,2])
        self.set_target(w0,"coupling","w",subproblem=[1,2])
        self.set_target(v0,"coupling","v",subproblem=[1,2])

        # default penalty weights
        self.gamma = gamma
        self.beta = beta
        
        self.cu = [np.zeros_like(ui_init) for ui_init in u_init] # copy state variables
        self.cx0 = [np.zeros_like(x0_init),]*self.n_disciplines # copy shared variables
        self.cu_old = [np.zeros_like(ui_init) for ui_init in u_init] # copy state variables
        self.cx0_old = [np.zeros_like(x0_init),]*self.n_disciplines # copy shared variables

        # penalty values
        self.n_evaluations = [0,]*(self.n_disciplines+1)

        self.opt_options = {
            "seed": 0, "budget": 100000, "tol": 1e-12, 
            "display": False, 
            "opportunistic": False
        } if opt_options == None else opt_options
        self.optimizer = OMADS.main

    def get_x(self,var_type:str,var_attribute:str=None) -> Union[np.ndarray,List[Any]]:
        """get shared variables"""

        variables = [var for var in self._variables if var["type"] == var_type and not var["target"]]

        if var_attribute in ["value","w","v"]:
            return np.array([var[var_attribute] for var in variables])
        elif var_attribute == None:
            return variables
        else:
            return [var[var_attribute] for var in variables]

    def set_x(self,x:Union[List,np.ndarray,float,int],var_type:str,var_attribute:str):
        """update shared variables"""
        var_index = [i for i,var in enumerate(self._variables) if var["type"] == var_type and not var["target"]]
        if isinstance(x,(list,np.ndarray)):
            assert len(x) == len(var_index)
            for i,xi in zip(var_index,x):
                self._variables[i][var_attribute] = xi
        elif isinstance(x,(float,int)):
            for i in var_index:
                self._variables[i][var_attribute] = x

    def get_target(self,var_type:str,var_attribute:str=None,subproblem:Union[int,List[int]]=None) -> Union[np.ndarray,List[Any]]:
        """get shared variables"""
        
        if subproblem is None:
            variables = [var for var in self._variables if var["type"] == var_type and var["target"]]
        elif isinstance(subproblem,int):
            variables = [var for var in self._variables if var["type"] == var_type and var["subproblem"] == subproblem and var["target"]]
        elif isinstance(subproblem,list):
            variables = [var for var in self._variables if var["type"] == var_type and var["subproblem"] in subproblem and var["target"]]

        if var_attribute in ["value","w","v"]:
            return np.array([var[var_attribute] for var in variables])
        elif var_attribute == None:
            return variables
        else:
            return [var[var_attribute] for var in variables]

    def set_target(self,x:Union[List,np.ndarray,float,int],var_type:str,var_attribute:str,subproblem:Union[int,List[int]]=None):
        """update shared variables"""

        if subproblem is None:
            var_index = [i for i,var in enumerate(self._variables) if var["type"] == var_type and var["target"]]
        elif isinstance(subproblem,int):
            var_index = [i for i,var in enumerate(self._variables) if var["type"] == var_type and var["subproblem"] == subproblem and var["target"]]
        elif isinstance(subproblem,list):
            var_index = [i for i,var in enumerate(self._variables) if var["type"] == var_type and var["subproblem"] in subproblem and var["target"]]

        if isinstance(x,(list,np.ndarray)):
            assert len(x) == len(var_index)
            for i,xi in zip(var_index,x):
                self._variables[i][var_attribute] = xi
        elif isinstance(x,(float,int)):
            for i in var_index:
                self._variables[i][var_attribute] = x

    def _flatten(self,l:List[List[float]]) -> List[float]:
        return [item for sublist in l for item in sublist]

    def _unflatten(self,list_flat:List[float],indices:List[int]) -> List[np.ndarray]:
        # un-flatten ut
        l = []; n_start = 0
        for n in indices:
            l += [np.array(list_flat[n_start:n+n_start])]
            n_start += n
        return l

    def penalty(self,ci:np.ndarray,v,w):
        return v.T @ np.abs(ci) + np.linalg.norm(w*ci)**2

    def bb_sys(self,vars):

        x0 = np.array(vars[:self.n_shared]) # shared variables (x_0) and state targets (u^t)
        ut_flat = vars[self.n_shared:] # shared variables (x_0) and state targets (u^t)
        ut = self._unflatten(ut_flat,self.n_states) # un-flatten ut

        # retrieve shared and local variables (x_0,x_i), and state variables (u) and compute global objective and constraints
        if callable(self.f_global):
            f = self.f_global(x0,self.xi,ut).astype("float64") # to avoid float128 dtype errors
        else:
            f = self.f_global

        g = []
        for gi in self.g_global:
            if callable(gi):
                g += [gi(x0,self.xi,ut).astype("float64")] # to avoid float128 dtype errors
            else:
                g += [gi]

        for i,discipline in enumerate(self.disciplines):
            ui = discipline(x0,self.xi,ut) # estimate coupling variables
            # find inconsistency
            phi_i = self.penalty(self.x0t[i] - x0,self.vx0[i],self.wx0[i]) + \
                self.penalty(ut[i] - ui, self.vu[i], self.wu[i])
            f += phi_i

        self.n_evaluations[0] += 1 # update iteration counter

        return [f,g]

    def bb_i(self,vars,params):

        i = params[0]
        # self:ATC = get_by_address(params[1])[0]
        # i is the discipline index
        x0it = np.array(vars[:self.n_shared]) # shared targets (x_0,i^t)
        xi = np.array(vars[self.n_shared:]) # local variables (x_i)

        # solve for coupling variables
        discipline_i = self.disciplines[i]
        u_hat = discipline_i(x0it,xi,self.ut) # estimate coupling variables

        f = 0
        # local objective
        if callable(self.f_local[i]):
            f += self.f_local[i](x0it,xi,u_hat).astype("float64") # to avoid float128 dtype errors
        else:
            f += self.f_local[i]

        # local constraints
        g = []
        for gi in self.g_local[i]:
            if callable(gi):
                g += [gi(x0it,xi,u_hat).astype("float64")]
            else:
                g += [gi]

        # find inconsistency
        phi_i = self.penalty(x0it - self.x0,self.vx0[i],self.wx0[i]) + \
            self.penalty(self.ut[i] - u_hat,self.vu[i],self.wu[i])
        f += phi_i

        self.n_evaluations[i+1] += 1 # update iteration counter

        return [f,g]

    def opt_bb_discipline(self,i,opt_options):

        # Optimization setup
        x0it_init = self.x0_init
        xi_init = self.xi_init[i]
        x_init = list(x0it_init) + list(xi_init) # initial local variables (x_i) and shared targets (x_0,i^t)
        
        lb = self.lbs + self.lbi[i] # lower bound of shared and local variable
        ub = self.ubs + self.ubi[i] # upper bound of shared and local variable

        eval = {"blackbox": self.bb_i, "constants": [i,]}
        x0it_labels = ["x_{0,%i,%i}^t"%(i+1,k+1) for k in range(len(x0it_init))]
        xi_labels = ["x_{%i,%i}"%(i+1,k+1) for k in range(len(xi_init))]

        param = {"baseline": x_init,
                "lb": lb,
                "ub": ub,
                "var_names": x0it_labels + xi_labels,
                "scaling": 10.0,
                "post_dir": "./post"} # these are OMADS specific options. You can modify them as necessary
        data = {"evaluator": eval, "param": param, "options":opt_options}

        out = self.optimizer(data)

        x0it = np.array(out["xmin"][0:self.n_shared])
        xi = np.array(out["xmin"][self.n_shared:])
        # update local variables (x_i) and shared targets (x_0,i^t)
        self.x0t[i] = x0it
        self.xi[i] = xi

        # solve for coupling variables
        discipline_i = self.disciplines[i]
        u_hat = discipline_i(x0it,xi,self.ut) # estimate coupling variables

        self.u[i] = u_hat # update state variable estimate globally

    def opt_sys(self,opt_options):

        # Optimization setup
        ut_init = self._flatten(self.u_init) # flatten states
        x0_init = self.x0_init
        x_init = list(x0_init) + list(ut_init) # initial shared variables (x_0) and state targets (u^t)
        
        flat_lbu = self._flatten(self.lbu) # flatten state bounds
        flat_ubu = self._flatten(self.ubu) # flatten state bounds
        lb = self.lbs + flat_lbu # lower bound of shared and state variables
        ub = self.ubs + flat_ubu # upper bound of shared and state variables

        eval = {"blackbox": self.bb_sys}
        x0_labels = ["x_{0,%i}"%(k+1) for k in range(len(x0_init))]
        u_labels =  ["u_%i,%i^t" %(i+1,k+1) for i,sublist in enumerate(self.u_init) for k,item in enumerate(sublist)]
        param = {"baseline": x_init,
                "lb": lb,
                "ub": ub,
                "var_names": x0_labels + u_labels,
                "scaling": 10.0,
                "post_dir": "./post"} # these are OMADS specific options. You can modify them as necessary
        data = {"evaluator": eval, "param": param, "options":opt_options}

        out = self.optimizer(data)

        x0 = np.array(out["xmin"][:self.n_shared])
        ut_flat = out["xmin"][self.n_shared:]

        # update shared (x_0) and state targets (u^t)
        self.x0 = x0
        self.ut = self._unflatten(ut_flat,self.n_states) # un-flatten ut

        # update penalties
        for i in range(self.n_disciplines):
            # find inconsistency
            phi_i = self.penalty(self.x0t[i] - self.x0,self.vx0[i],self.wx0[i]) + \
                self.penalty(self.ut[i] - self.u[i], self.vu[i], self.wu[i])
            self.phi_old[i] = self.phi[i]
            self.phi[i] = phi_i

            # update inconsistencies
            self.cu[i] = self.ut[i] - self.u[i]
            self.cx0[i] = self.x0t[i] - self.x0

        return float(out["fmin"])

    def run_atc(self,k_max=100,tol=1e-6):
        
        k_outer = 0
        while True:

            # run sub-system optimization first
            for i in range(self.n_disciplines):
                self.opt_bb_discipline(i,self.opt_options)

            # run the system optimization loop
            self.opt_sys(self.opt_options)

            # check outer loop convergence
            if k_outer > 1:
                c_vector = self._flatten(self.cu)+self._flatten(self.cx0) # flatten states
                if k_outer > k_max or max([abs(ci) for ci in c_vector]) <= tol:
                    break

            # update penalty weights
            for i in range(self.n_disciplines):
                self.vx0[i] += 2*((self.wx0[i])**2)*np.abs(self.cx0[i])
                cond = np.abs(self.cx0[i]) > self.gamma*np.abs(self.cx0_old[i])
                self.wx0[i][cond] *= self.beta
                self.cx0_old[i] = self.cx0[i]

                self.vu[i] += 2*((self.wu[i])**2)*np.abs(self.cu[i])
                cond = np.abs(self.cu[i]) > self.gamma*np.abs(self.cu_old[i])
                self.wu[i][cond] *= self.beta
                self.cu_old[i] = self.cu[i]

            k_outer += 1 # update iteration counter

        # retrieve local variables (x_i), shared variables (x_0), and state variables (u)
        f = self.f_global(self.x0,self.xi,self.u) # calculate global objective

        print("ATC terminated")

        return self.x0,self.cx0,self.xi,self.u,self.cu,self.phi,f

    def __deepcopy__(self, memo): # memo is a dict of id's to copies
        """
        creates a deep independent copy of the class instance self.
        https://stackoverflow.com/a/15774013
        """
        id_self = id(self)
        _copy = memo.get(id_self) # memoization avoids unnecessary recursion
        if _copy is None:
            _copy = type(self)(self.disciplines,self.f_global,self.g_global,self.f_local,
            self.g_local,self.x0_init,self.xi_init,self.u_init,
            self.lbs,self.ubs,self.lbi,self.ubi,self.lbu,self.ubu,
            self.gamma,self.beta,opt_options=self.opt_options)

            # initialization
            _copy.x0 = self.x0
            _copy.xi = self.xi
            _copy.u = self.u
            _copy.ut = self.ut
            _copy.x0t = self.x0t
            
            # inconsistency to updated during coordination
            _copy.wu = self.wu
            _copy.wx0 = self.wx0
            _copy.vu = self.vu
            _copy.vx0 = self.vx0
            _copy.cu = self.cu
            _copy.cx0 = self.cx0
            _copy.cu_old = self.cu_old
            _copy.cx0_old = self.cx0_old

            # penalty values
            _copy.F_total_old = self.F_total_old
            _copy.F_total = self.F_total
            _copy.phi_old = self.phi_old
            _copy.phi = self.phi
            _copy.n_evaluations = self.n_evaluations

            memo[id_self] = _copy 
        return _copy

if __name__ == "__main__":

    def prob0(x0,x,u):
        """This is the system level problem"""
        u1 = u[0] # coupling variable from SP1
        u2 = u[1] # coupling variable from SP2
        x1 = x[0] # local variable from SP1
        x2 = x[1] # local variable from SP2 (unused)

        f = x0[0] + x1[0] + u1[0] + u2[0] # global objective
        g = [0,] # global constraints (none)
        return [f,g]

    def prob1(x0,xi,ui):
        """This is the subsystem 1 problem"""
        f = 0 # local objective (none)
        g = [0,] # local constraints (none)
        return [f,g]

    def prob2(x0,xi,ui):
        """This is the subsystem 2 problem"""
        f = 0 # local objective (none)
        g = [xi[0] + ui[0] - 10,] # local constraints
        return [f,g]

    # these are your coupled disciplines
    def discpline_1(x0,xi,u):
        """This is the discpline 1 analysis"""
        u2 = u[1] # target variable coming from other discipline
        u1 = [np.log(x0[0]+1e-6) + np.log(xi[0]+1e-6) + np.log(u2[0]+1e-6),] # 1e-6 to avoid log(0) errors
        return np.array(u1)

    def discpline_2(x0,xi,u):
        """This is the discpline 2 analysis"""
        u1 = u[0] # target variable coming from other discipline
        u2 = [1/(x0[0]+1e-6) + 1/(xi[0]+1e-6) + 1/(u1[0]+1e-6),] # 1e-6 to avoid division by zero errors
        return np.array(u2)

    # Variable definitions
    variables = [
        {"name": "u", "type": "shared", "subproblem": 0, "target": False, "value": 1.0, "lb": 0.0, "ub": 10.0},
        {"name": "v", "type": "local", "subproblem": 1, "target": False,  "value": 1.0, "lb": 0.0, "ub": 10.0},
        {"name": "w", "type": "local", "subproblem": 2, "target": False,  "value": 1.0, "lb": 0.0, "ub": 10.0},
        {"name": "a", "type": "coupling", "subproblem": 1, "target": False,  "value": 1.0, "lb": 0.0, "ub": 10.0},
        {"name": "b", "type": "coupling", "subproblem": 2, "target": False,  "value": 1.0, "lb": 0.0, "ub": 10.0},
    ]
    # optimization algorithm options
    opt_options = {
        "seed": 0, 
        "budget": 1000, 
        "tol": 1e-12, 
        "display": False, 
        "opportunistic": False
    }
    my_atc = ATC(disciplines=[discpline_1,discpline_2],P_global=prob0, 
        P_local=[prob1,prob2],variables=variables,opt_options=opt_options,w0=1.0,v0=0.0,
        beta=1.1,gamma=0.1)
    x0,cx0,xi,u,cu,phi,f = my_atc.run_atc(k_max=200,tol=1e-6)

    print(f)
    print(list(x0)+my_atc._flatten(xi))
    print(my_atc._flatten(u))
    print(my_atc._flatten(cx0) + my_atc._flatten(cu))
    g = g21(x0,xi,u[1])
    print(g)