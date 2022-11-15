from typing import List, Dict, Any
import OMADS
import numpy as np
import json

def get_by_address(address):
    return [x for x in globals().values() if id(x)==address]

class ATC():

    def __init__(self,disciplines:List[callable],f_global:callable,g_global:callable,f_local:List[callable],
        g_local:List[callable],xs_init:List[float],xi_init:List[float],u_init:List[List[float]],
        lbs:List[float],ubs:List[float],lbi:List[float],ubi:List[float],lbu:List[List[float]],ubu:List[List[float]],
        gamma:float=0.8,beta:float=1.3,w0:float=0.01,v0:float=0.01,
        in_opt_options:Dict[str,Any]=None):

        self.disciplines = disciplines
        self.n_disciplines = len(disciplines)
        self.f_local = f_local
        self.g_local = g_local
        self.f_global = f_global
        self.g_global = g_global
        
        # shared variables
        self.xs_init = np.array(xs_init)
        self.lbs = lbs
        self.ubs = ubs
        self.n_shared = len(xs_init)

        # local variables
        self.xi_init = [np.array(xik_init) for xik_init in xi_init]
        self.lbi = lbi
        self.ubi = ubi
        self.n_local = [len(xi) for xi in xi_init]

        # state variable
        self.u_init = u_init
        self.lbu = lbu
        self.ubu = ubu
        self.n_states = [len(ui) for ui in u_init]

        # default penalty weights
        self.gamma = gamma
        self.beta = beta
        self.w = [w0,]*self.n_disciplines
        self.v = [v0,]*self.n_disciplines

        # initialization
        self.xs = [np.array(xs_init),]*self.n_disciplines
        self.xi = [np.array(xik_init) for xik_init in xi_init]
        self.u = [np.array(ui_init) for ui_init in u_init]
        self.ut = [np.array(ui_init) for ui_init in u_init] # copy state variables
        self.x0t = [np.array(xs_init),]*self.n_disciplines # copy shared variables
        
        # inconsistency to updated during coordination
        self.cu = [np.zeros_like(ui_init) for ui_init in u_init] # copy state variables
        self.cxs = [np.zeros_like(xs_init),]*self.n_disciplines # copy shared variables

        # penalty values
        self.phi_old = [0.,]*self.n_disciplines
        self.phi = [0.,]*self.n_disciplines
        self.k_outer = 0
        self.k_inner = [0,]*self.n_disciplines

        self.in_opt_options = {
            "seed": 0, "budget": 100000, "tol": 1e-12, 
            "display": False, 
            "opportunistic": False
        } if in_opt_options == None else in_opt_options
        self.optimizer = OMADS.main

    def _flatten(self,l:List[List[float]]) -> List[float]:
        return [item for sublist in l for item in sublist]

    def _unflatten(self,list_flat:List[float],indices:List[int]) -> List[np.ndarray]:
        # un-flatten ut
        l = []; n_start = 0
        for n in indices:
            l += [np.array(list_flat[n_start:n+n_start])]
            n_start += n
        return l

    def penalty(self,hi:np.ndarray,i):
        return np.sum(self.v[i]*np.abs(hi) + ((self.w[i]*hi)**2))

    def bb_sys(self,vars,params):

        # self:ATC = get_by_address(params[0])[0]

        x0 = np.array(vars[:self.n_shared]) # shared variables (x_0) and state targets (u^t)
        ut_flat = vars[self.n_shared:] # shared variables (x_0) and state targets (u^t)

        # update shared (x_0) and state targets (u^t)
        self.xs = x0
        self.ut = self._unflatten(ut_flat,self.n_states) # un-flatten ut

        # run sub-system optimization first
        for i in range(self.n_disciplines):
            self.opt_bb_discipline(i,self.in_opt_options)

        # retrieve shared and local variables (x_0,x_i), and state variables (u) and compute global objective and constraints
        if callable(self.f_global):
            f = self.f_global(self.xs,self.xi,self.u).astype("float64") # to avoid float128 dtype errors
        else:
            f = self.f_global

        g = []
        for gi in self.g_global:
            if callable(gi):
                g += [gi(self.xs,self.xi,self.u).astype("float64")] # to avoid float128 dtype errors
            else:
                g += [gi]

        for i in range(self.n_disciplines):
            f += self.phi[i]

        # update penalty
        for i in range(len(self.v)):
            self.v[i] = min(self.v[i] + min(2*((self.w[i])**2)*self.phi[i],1),1e200)
            if np.abs(self.phi[i]) > self.gamma*np.abs(self.phi_old[i]):
                self.w[i] = min(self.w[i]*self.beta,1e200)

        self.k_outer += 1 # update iteration counter

        dict_outer = {
            "k_outer" : self.k_outer,
            "v" : [float(vi) for vi in self.v],
            "w": [float(wi) for wi in self.w]
        }
        # write to log file
        with open("outer_loop.json", "w") as file:
            json.dump(dict_outer,file,indent=4)

        return [f,g]

    def bb_i(self,vars,params):

        i = params[0]
        # self:ATC = get_by_address(params[1])[0]
        # i is the discipline index
        x0it = np.array(vars[:self.n_shared]) # shared targets (x_0,i^t)
        xi = np.array(vars[self.n_shared:]) # local variables (x_i)

        # retrieve shared (x_0) and state targets (u^t)
        x0 = np.array(self.xs)
        ut = self.ut

        # solve for coupling variables
        discipline_i = self.disciplines[i]
        u_hat = discipline_i(x0,xi,self.ut) # estimate coupling variables

        # find inconsistency
        phi_i = self.penalty(x0it - x0,i) + self.penalty(ut[i] - u_hat,i)

        if callable(self.f_local[i]):
            f = self.f_local[i](x0it,xi,u_hat).astype("float64") # to avoid float128 dtype errors
        else:
            f = self.f_local[i]

        g = []
        for gi in self.g_local[i]:
            if callable(gi):
                g += [gi(x0it,xi,u_hat).astype("float64")]
            else:
                g += [gi]

        # local objective
        f += phi_i
        # f += self.penalty(np.array(g),i) # to avoid optimization issues with constraint

        self.k_inner[i] += 1 # update iteration counter

        return [f,g] # to avoid float128 dtype errors

    def opt_bb_discipline(self,i,opt_options):

        # Optimization setup
        x0it_init = self.xs_init
        xi_init = self.xi_init[i]
        x0 = list(x0it_init) + list(xi_init) # initial local variables (x_i) and shared targets (x_0,i^t)
        
        lb = self.lbs + self.lbi[i] # lower bound of shared and local variable
        ub = self.ubs + self.ubi[i] # upper bound of shared and local variable

        eval = {"blackbox": self.bb_i, "constants": [i,id(self)]}
        x0it_labels = ["x_{0,%i,%i}^t"%(i+1,k+1) for k in range(len(x0it_init))]
        xi_labels = ["x_{%i,%i}"%(i+1,k+1) for k in range(len(xi_init))]

        param = {"baseline": x0,
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

        # update penalties
        phi_i = self.penalty(x0it - self.xs,i) + self.penalty(self.ut[i] - u_hat,i)
        self.phi_old[i] = self.phi[i]
        self.phi[i] = phi_i

        # update inconsistencies
        self.cu[i] = self.ut[i] - u_hat
        self.cxs[i] = x0it - self.xs

    def opt_sys(self,opt_options:Dict[str,Any]=None,current_id=None):
        
        current_id = current_id if current_id is not None else id(self)
        opt_options = {
            "seed": 0, 
            "budget": 10000, 
            "tol": 1e-12, 
            "display": False,
            "opportunistic": False
        } if opt_options == None else opt_options

        # Optimization setup
        ut_init = self._flatten(self.u_init) # flatten states
        x0_init = self.xs_init
        x0 = list(x0_init) + list(ut_init) # initial shared variables (x_0) and state targets (u^t)
        
        flat_lbu = self._flatten(self.lbu) # flatten state bounds
        flat_ubu = self._flatten(self.ubu) # flatten state bounds
        lb = self.lbs + flat_lbu # lower bound of shared and state variables
        ub = self.ubs + flat_ubu # upper bound of shared and state variables

        eval = {"blackbox": self.bb_sys, "constants": [current_id,]}
        x0_labels = ["x_{0,%i}"%(k+1) for k in range(len(x0_init))]
        u_labels =  ["u_%i,%i^t" %(i+1,k+1) for i,sublist in enumerate(self.u_init) for k,item in enumerate(sublist)]
        param = {"baseline": x0,
                "lb": lb,
                "ub": ub,
                "var_names": x0_labels + u_labels,
                "scaling": 10.0,
                "post_dir": "./post"} # these are OMADS specific options. You can modify them as necessary
        data = {"evaluator": eval, "param": param, "options":opt_options}

        out = self.optimizer(data)

        x0 = np.array(out["xmin"][:self.n_shared])
        ut_flat = out["xmin"][self.n_shared:]

        # write to log file
        with open("outer_loop.json", "r") as file:
            dict_outer = json.load(file)

        self.k_outer = dict_outer["k_outer"]
        self.v = dict_outer["v"]
        self.w = dict_outer["w"]

        # RERUN everything at converged solution
        # update shared (x_0) and state targets (u^t)
        self.xs = x0
        self.ut = self._unflatten(ut_flat,self.n_states) # un-flatten ut

        eval_opt_options = {"seed": 0, "budget": 1000000000, "tol": 1e-12, "display": False} if opt_options == None else opt_options
        # run sub-system optimization first
        for i in range(self.n_disciplines):
            self.opt_bb_discipline(i,eval_opt_options)

        # retrieve local variables (x_i), shared targets (x_0,i^t), and state variables (u)
        f = self.f_global(self.xs,self.xi,self.u)

        print("ATC terminated")

        return self.xs,self.cxs,self.xi,self.u,self.cu,self.phi,f

    def __deepcopy__(self, memo): # memo is a dict of id's to copies
        """
        creates a deep independent copy of the class instance self.
        https://stackoverflow.com/a/15774013
        """
        id_self = id(self)
        _copy = memo.get(id_self) # memoization avoids unnecessary recursion
        if _copy is None:
            _copy = type(self)(self.disciplines,self.f_global,self.g_global,self.f_local,
            self.g_local,self.xs_init,self.xi_init,self.u_init,
            self.lbs,self.ubs,self.lbi,self.ubi,self.lbu,self.ubu,
            self.gamma,self.beta,in_opt_options=self.in_opt_options)

            # initialization
            _copy.xs = self.xs
            _copy.xi = self.xi
            _copy.u = self.u
            _copy.ut = self.ut
            _copy.x0t = self.x0t
            
            # inconsistency to updated during coordination
            _copy.cu = self.cu
            _copy.cxs = self.cxs

            # penalty values
            _copy.phi_old = self.phi_old
            _copy.phi = self.phi
            _copy.k_outer = self.k_outer
            _copy.k_inner = self.k_inner

            memo[id_self] = _copy 
        return _copy

if __name__ == "__main__":

    # these are your coupled disciplines
    U1 = lambda x0,xi,u: np.array([np.log(x0[0]+1e-6) + np.log(xi[0]+1e-6) + np.log(u[1][0]+1e-6)])
    U2 = lambda x0,xi,u: np.array([1/(x0[0]+1e-6) + 1/(xi[0]+1e-6) + 1/(u[0][0]+1e-6)])
    f = lambda xs,x,u: xs[0] + x[0][0] + u[0][0] + u[1][0]
    g = lambda xs,x,u: x[1][0] + u[1][0] - 10 # is not using the shared variable (but could in practice)
    g21 = lambda xs,xi,ui: xi[0] + ui[0] - 10 # is not using the shared variable (but could in practice)

    # system optimization problem
    f_global = f
    g_global = [0,] # no global constraints
    # sub-system 1 optimization problem
    f1 = 0 # no local objective
    g1 = [0,] # no local constraints
    # sub-system 2 optimization problem
    f2 = 0 # no local objective
    g2 = [g21,]

    f_local = [f1,f2]
    g_local = [g1,g2]
    disciplines = [U1,U2]

    u0 = 5.0; v0 =5.0; w0 = 5.0; a0 = 5.0; b0 = 5.0
    xs_init = [u0,] # initial guess for shared variable (u)
    xi_init = [[v0,],[w0,]] # initial guess for local variables (v and w)
    u_init = [[a0,],[b0,]] # initial guess for state variables (a and b)

    lbs = [0.0,]; ubs = [10.,] # bounds of shared variable (u)
    lbi = [[0.0,], [0.0,]]; ubi = [[10.,], [10.,]] # bounds of local variables (v and w)
    lbu = [[0.01,], [0.01,]]; ubu = [[20.,], [20.,]] # bounds of state variables (a and b)

    in_opt_options = {
        "seed": 0, 
        "budget": 500, 
        "tol": 1e-12, 
        "display": False, 
        "opportunistic": False
    }
    out_opt_options = {
        "seed": 0, 
        "budget": 200, 
        "tol": 1e-12, 
        "display": False,
        "opportunistic": False
    }
    my_atc = ATC(disciplines=disciplines,f_global=f_global, 
        g_global=g_global,f_local=f_local,g_local=g_local,xs_init=xs_init,
        xi_init=xi_init,u_init=u_init,lbs=lbs,ubs=ubs,lbi=lbi,ubi=ubi,
        lbu=lbu,ubu=ubu,in_opt_options=in_opt_options)
    my_atc.opt_sys(out_opt_options)