from typing import List, Dict, Any
import OMADS
import numpy as np
import json

def get_by_address(address):
    return [x for x in globals().values() if id(x)==address]

class ATC():

    def __init__(self,disciplines:List[callable],f_global:callable,g_global:callable,f_local:List[callable],
        g_local:List[callable],x0_init:List[float],xi_init:List[float],u_init:List[List[float]],
        lbs:List[float],ubs:List[float],lbi:List[float],ubi:List[float],lbu:List[List[float]],ubu:List[List[float]],
        gamma:float=0.1,beta:float=1.3,w0:float=0.5,v0:float=0.5,
        opt_options:Dict[str,Any]=None):

        self.disciplines = disciplines
        self.n_disciplines = len(disciplines)
        self.f_local = f_local
        self.g_local = g_local
        self.f_global = f_global
        self.g_global = g_global
        
        # shared variables
        self.x0_init = np.array(x0_init)
        self.lbs = lbs
        self.ubs = ubs
        self.n_shared = len(x0_init)

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

        # initialization
        self.x0 = np.array(x0_init)
        self.xi = [np.array(xik_init) for xik_init in xi_init]
        self.u = [np.array(ui_init) for ui_init in u_init]
        self.ut = [np.array(ui_init) for ui_init in u_init] # copy state variables
        self.x0t = [np.array(x0_init),]*self.n_disciplines # copy shared variables
        
        # inconsistency to updated during coordination
        self.wx0 = [w0*np.ones(len(x0_init)),]*self.n_disciplines
        self.wu = [w0*np.ones(len(ui_init)) for ui_init in u_init]
        self.vx0 = [v0*np.ones(len(x0_init)),]*self.n_disciplines
        self.vu = [v0*np.ones(len(ui_init)) for ui_init in u_init]
        
        self.cu = [np.zeros_like(ui_init) for ui_init in u_init] # copy state variables
        self.cx0 = [np.zeros_like(x0_init),]*self.n_disciplines # copy shared variables
        self.cu_old = [np.zeros_like(ui_init) for ui_init in u_init] # copy state variables
        self.cx0_old = [np.zeros_like(x0_init),]*self.n_disciplines # copy shared variables

        # penalty values
        self.F_total_old = 0.
        self.F_total = 0.
        self.phi_old = [0.,]*self.n_disciplines
        self.phi = [0.,]*self.n_disciplines
        self.k_sys = 0
        self.k_disp = [0,]*self.n_disciplines

        self.opt_options = {
            "seed": 0, "budget": 100000, "tol": 1e-12, 
            "display": False, 
            "opportunistic": False
        } if opt_options == None else opt_options
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

    def penalty(self,ci:np.ndarray,v,w):
        return v.T @ np.abs(ci) + np.linalg.norm(w*ci)**2

    def bb_sys(self,vars,params):

        self:ATC = get_by_address(params[0])[0]

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

        self.k_sys += 1

        return [f,g]

    def bb_i(self,vars,params):

        i = params[0]
        self:ATC = get_by_address(params[1])[0]
        # i is the discipline index
        x0it = np.array(vars[:self.n_shared]) # shared targets (x_0,i^t)
        xi = np.array(vars[self.n_shared:]) # local variables (x_i)

        # solve for coupling variables
        discipline_i = self.disciplines[i]
        u_hat = discipline_i(x0it,xi,self.ut) # estimate coupling variables

        # global objective
        xi_all = self.xi.copy(); xi_all[i] = xi
        u_all = self.ut.copy(); u_all[i] = u_hat
        if callable(self.f_global):
            f = self.f_global(x0it,xi_all,u_all).astype("float64") # to avoid float128 dtype errors
        else:
            f = self.f_global
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

        self.k_disp[i] += 1 # update iteration counter

        return [f,g] # to avoid float128 dtype errors

    def opt_bb_discipline(self,i,opt_options):

        # Optimization setup
        x0it_init = self.x0_init
        xi_init = self.xi_init[i]
        x_init = list(x0it_init) + list(xi_init) # initial local variables (x_i) and shared targets (x_0,i^t)
        
        lb = self.lbs + self.lbi[i] # lower bound of shared and local variable
        ub = self.ubs + self.ubi[i] # upper bound of shared and local variable

        eval = {"blackbox": self.bb_i, "constants": [i,id(self)]}
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

        # # update penalties
        # phi_i = self.penalty(x0it - self.x0,i) + self.penalty(self.ut[i] - u_hat,i)
        # self.phi_old[i] = self.phi[i]
        # self.phi[i] = phi_i

        # # update inconsistencies
        # self.cu[i] = self.ut[i] - u_hat
        # self.cx0[i] = x0it - self.x0

    def opt_sys(self,opt_options):

        # Optimization setup
        ut_init = self._flatten(self.u_init) # flatten states
        x0_init = self.x0_init
        x_init = list(x0_init) + list(ut_init) # initial shared variables (x_0) and state targets (u^t)
        
        flat_lbu = self._flatten(self.lbu) # flatten state bounds
        flat_ubu = self._flatten(self.ubu) # flatten state bounds
        lb = self.lbs + flat_lbu # lower bound of shared and state variables
        ub = self.ubs + flat_ubu # upper bound of shared and state variables

        eval = {"blackbox": self.bb_sys, "constants": [id(self),]}
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
            k_inner = 0
            # while True:
            #     # run sub-system optimization first
            #     for i in range(self.n_disciplines):
            #         self.opt_bb_discipline(i,self.opt_options)

            #     # run the system optimization loop
            #     self.F_total = self.opt_sys(self.opt_options)

            #     # check inner loop convergence
            #     if k_inner > 1:
            #         if k_inner > k_max_inner or np.abs(self.F_total - self.F_total_old) / (np.abs(self.F_total) + 1) <= tol_outer:
            #             break

            #     self.F_total_old = self.F_total
            #     k_inner += 1

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

            # update penalty
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

        # retrieve local variables (x_i), shared targets (x_0,i^t), and state variables (u)
        f = self.f_global(self.x0,self.xi,self.u)

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
            _copy.k_sys = self.k_sys
            _copy.k_disp = self.k_disp

            memo[id_self] = _copy 
        return _copy

if __name__ == "__main__":

    # these are your coupled disciplines
    U1 = lambda x0,xi,u: np.array([np.log(x0[0]+1e-6) + np.log(xi[0]+1e-6) + np.log(u[1][0]+1e-6)])
    U2 = lambda x0,xi,u: np.array([1/(x0[0]+1e-6) + 1/(xi[0]+1e-6) + 1/(u[0][0]+1e-6)])
    f = lambda x0,x,u: x0[0] + x[0][0] + u[0][0] + u[1][0]
    g = lambda x0,x,u: x[1][0] + u[1][0] - 10 # is not using the shared variable (but could in practice)
    g21 = lambda x0,xi,ui: xi[0] + ui[0] - 10 # is not using the shared variable (but could in practice)

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

    u0 = 1.0; v0 = 1.0; w0 = 1.0; a0 = 1.0; b0 = 1.0
    x0_init = [u0,] # initial guess for shared variable (u)
    xi_init = [[v0,],[w0,]] # initial guess for local variables (v and w)
    u_init = [[a0,],[b0,]] # initial guess for state variables (a and b)

    lbs = [0.0,]; ubs = [10.,] # bounds of shared variable (u)
    lbi = [[0.0,], [0.0,]]; ubi = [[10.,], [10.,]] # bounds of local variables (v and w)
    lbu = [[0.0,], [0.0,]]; ubu = [[10.,], [10.,]] # bounds of state variables (a and b)

    opt_options = {
        "seed": 0, 
        "budget": 1000, 
        "tol": 1e-12, 
        "display": False, 
        "opportunistic": False
    }
    my_atc = ATC(disciplines=disciplines,f_global=f_global, 
        g_global=g_global,f_local=f_local,g_local=g_local,x0_init=x0_init,
        xi_init=xi_init,u_init=u_init,lbs=lbs,ubs=ubs,lbi=lbi,ubi=ubi,
        lbu=lbu,ubu=ubu,opt_options=opt_options,w0=0.1,v0=0.1,
        beta=1.1,gamma=0.1)
    x0,cx0,xi,u,cu,phi,f = my_atc.run_atc(k_max=200,tol=1e-6)

    print(f)
    print(list(x0)+my_atc._flatten(xi))
    print(my_atc._flatten(u))
    print(my_atc._flatten(cx0) + my_atc._flatten(cu))
    g = g21(x0,xi,u[1])
    print(g)