from typing import List
import OMADS
import numpy as np

class ATC():

    def __init__(self,disciplines:List[callable],x_init:List[float],u_init:List[float],lb:List[float],ub:List[float],
        lb_u:List[float],ub_u:List[float],gamma:float=0.1,beta:float=1.3,w0:List[float]=None,v0:List[float]=None):
        self.disciplines = disciplines

        self.x_init = x_init
        self.lb = lb
        self.ub = ub
        self.u_init = u_init
        self.lb_u = lb_u
        self.ub_u = ub_u

        # default penalty weights
        self.gamma = gamma
        self.beta = beta
        self.w = [0.01,]*len(disciplines) if w0 is None else w0
        self.v = [0.01,]*len(disciplines) if v0 is None else v0

        # initialization
        self.x = self.x_init
        self.u = np.array(self.u_init)
        self.x0t = [self.x_init[0],]*len(disciplines) # copy shared variables
        
        self.phi_old = [1,]*len(disciplines)
        self.phi = [1,]*len(disciplines)
        self.k_outer = 0
        self.k_inner = [0,]*len(disciplines)

    def penalty(self,hi,i):
        return self.v[i]*np.abs(hi) + (self.w[i]*hi)**2

    def bb_1(self,vars):
        i = 0

        [x1,x01t] = vars # local variables (x_i) and shared targets (x_0,i^t)

        # retrieve shared (x_0) and state targets (u^t)
        x0 = self.x[0]
        ut = self.ut

        # solve for coupling variables
        discipline_1 = self.disciplines[i]
        mask = [True,]*len(self.disciplines)
        mask[i] = False # select all state variables but current one
        u_hat = discipline_1([x0,],[x1,],self.ut[mask]) # estimate coupling variables

        # find inconsistency
        phi1 = self.penalty(x01t - x0,i) + self.penalty(ut[i] - u_hat,i)

        # local objective
        f = phi1

        self.k_inner[i] += 1 # update iteration counter

        return f.astype("float64") # to avoid float128 dtype errors

    def bb_2(self,vars):
        i = 1

        [x2,x02t] = vars # local variables (x_i) and shared targets (x_0,i^t)

        # retrieve shared (x_0) and state targets (u^t)
        x0 = self.x[0]
        ut = self.ut

        # solve for coupling variables
        discipline_2 = self.disciplines[i]
        mask = [True,]*len(self.disciplines)
        mask[i] = False # select all state variables but current one
        u_hat = discipline_2([x0,],[x2,],self.ut[mask]) # estimate coupling variables

        # find inconsistency
        phi2 = self.penalty(x02t - x0,i) + self.penalty(ut[i] - u_hat,i)

        # local objective and constraint
        f = phi2
        g2 = x2 + ut[mask] - 10

        self.k_inner[i] += 1 # update iteration counter

        return [f.astype("float64"),g2.astype("float64")] # to avoid float128 dtype errors

    def bb_sys(self,vars):

        x0 = vars[0] # shared variables (x_0) and state targets (u^t)
        ut = vars[1:] # shared variables (x_0) and state targets (u^t)

        # update shared (x_0) and state targets (u^t)
        self.x[0] = x0
        self.ut = np.array(ut)

        # run sub-system optimization first
        self.opt_bb_discipline(0)
        self.opt_bb_discipline(1)

        # retrieve local variables (x_i), shared targets (x_0,i^t), and state variables (u)
        x1,x2 = self.x[1:]
        x01t,x02t = self.x0t # not required here
        u_hat = self.u

        f0 = x0 + x1 + u_hat[0] + u_hat[1]
        f = f0 + self.phi[0] + self.phi[1]

        # update penalty
        for i in range(len(self.v)):
            self.v[i] += min(2*((self.w[i])**2)*self.phi[i],0.1)
            if np.abs(self.phi[i]) > self.gamma*np.abs(self.phi_old[i]):
                self.w[i] *= self.beta

        self.k_outer += 1 # update iteration counter

        return f.astype("float64") # to avoid float128 dtype errors

    def opt_bb_discipline(self,i):

        # Optimization setup
        xi_init = self.x_init[i+1]
        x0it_init = self.x_init[0]
        x0 = [xi_init,x0it_init] # initial local variables (x_i) and shared targets (x_0,i^t)
        
        lb = [self.lb[0],self.lb[i]] # lower bound of shared and local variable
        ub = [self.ub[0],self.ub[i]] # upper bound of shared and local variable

        if i == 0:
            eval = {"blackbox": self.bb_1}
        elif i == 1:
            eval = {"blackbox": self.bb_2}

        param = {"baseline": x0,
                "lb": lb,
                "ub": ub,
                "var_names": ["x_%i"%i, "x_{0,%i}^t"%(i+1)],
                "scaling": 10.0,
                "post_dir": "./post"} # these are OMADS specific options. You can modify them as necessary
        options = {"seed": 0, "budget": 100000, "tol": 1e-12, "display": False}

        data = {"evaluator": eval, "param": param, "options":options}

        out = {}
        # out is a dictionary that will hold output data of the final solution. The out dictionary has three keys: "xmin", "fmin" and "hmin"

        out = OMADS.main(data)

        xi = out["xmin"][0]
        x0it = out["xmin"][1]
        # update local variables (x_i) and shared targets (x_0,i^t)
        self.x[i+1] = xi
        self.x0t[i] = x0it

        # solve for coupling variables
        discipline_i = self.disciplines[i]
        mask = [True,]*len(self.disciplines)
        mask[i] = False # select all state variables but current one
        u_hat = discipline_i([x0it,],[xi,],self.ut[mask]) # estimate coupling variables

        self.u[i] = u_hat # update state variable estimate globally

        # update inconsistency
        phi_i = self.penalty(x0it - self.x[0],i) + self.penalty(self.ut[i] - u_hat,i)
        self.phi_old[i] = self.phi[i]
        self.phi[i] = phi_i

    def opt_sys(self):
        
        # Optimization setup
        ut_init = self.u_init
        x0_init = self.x_init[0]
        x0 = [x0_init,] + ut_init # initial shared variables (x_0) and state targets (u^t)
        
        lb = [self.lb[0],] + self.lb_u # lower bound of shared and state variables
        ub = [self.ub[0],] + self.ub_u # upper bound of shared and state variables

        eval = {"blackbox": self.bb_sys}
        param = {"baseline": x0,
                "lb": lb,
                "ub": ub,
                "var_names": ["x_0",] +  ["u_%i^t" %(i+1) for i in range(len(self.u_init))],
                "scaling": 10.0,
                "post_dir": "./post"} # these are OMADS specific options. You can modify them as necessary
        options = {"seed": 0, "budget": 100000, "tol": 1e-12, "display": False}

        data = {"evaluator": eval, "param": param, "options":options}

        out = {}
        # out is a dictionary that will hold output data of the final solution. The out dictionary has three keys: "xmin", "fmin" and "hmin"

        out = OMADS.main(data)

        x0 = out["xmin"][0]
        ut = out["xmin"][1:]

        # RERUN everything at converged solution
        # update shared (x_0) and state targets (u^t)
        self.x[0] = x0
        self.ut = np.array(ut)

        # run sub-system optimization first
        self.opt_bb_discipline(0)
        self.opt_bb_discipline(1)

        # retrieve local variables (x_i), shared targets (x_0,i^t), and state variables (u)
        x1,x2 = self.x[1:]
        u_hat = self.u
        f0 = x0 + x1 + u_hat[0] + u_hat[1]
        f = f0 + self.phi[0] + self.phi[1]

        print("ATC terminated")

        return self.x,self.u,self.phi,f

if __name__ == "__main__":

    # these are your coupled disciplines
    U1 = lambda x0,xi,uj: np.log(x0[0]+1e-6) + np.log(xi[0]+1e-6) + np.log(uj[0]+1e-6)
    U2 = lambda x0,xi,uj: 1/(x0[0]+1e-6) + 1/(xi[0]+1e-6) + 1/(uj[0]+1e-6)

    disciplines = [U1,U2]

    u0 = 1.0; v0 = 1.0; w0 = 1.0
    x_init = [u0,v0,w0]

    a0 = 1.0; b0 = 1.0
    u_init = [a0,b0]

    lb = [0.0, 0.0, 0.0]
    ub = [10., 10., 10.]
    lb_u = [1e-1, 1e-1]
    ub_u = [10., 10.]

    my_atc = ATC(disciplines=disciplines,x_init=x_init,u_init=u_init,lb=lb,ub=ub,lb_u=lb_u,ub_u=ub_u)
    my_atc.opt_sys()