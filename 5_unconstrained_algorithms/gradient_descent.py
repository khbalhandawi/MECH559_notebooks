import numpy as np
from typing import Callable, Tuple
import warnings
from copy import copy

from numpy.linalg import LinAlgError
class Optimization():

    def __init__(self,f:Callable[[np.ndarray,], np.ndarray],grad_f:Callable[[np.ndarray,], np.ndarray],
        H_f:Callable[[np.ndarray,], np.ndarray]=None,x0:np.ndarray=np.array([[0.0,0.0]]),verbose:bool=False):

        self.f = f
        self.grad_f = grad_f
        self.H_f = H_f
        self.xk = copy(x0)
        self.ndim = x0.shape[0]
        self.verbose = verbose

        self.fk = self.f(self.xk)
        self.grad_fk = self.grad_f(self.xk).squeeze(axis=2)
        if self.H_f is not None: self.Hk = self.H_f(self.xk).squeeze(axis=2)


class GradientDesc(Optimization):

    def __init__(self,f:Callable[[np.ndarray,], np.ndarray],grad_f:Callable[[np.ndarray,], np.ndarray],
        H_f:Callable[[np.ndarray,], np.ndarray]=None,x0:np.ndarray=np.array([[0.0,0.0]]),
        epsilon:float=0.1,beta:float=0.9,direction:str="steepest",alpha:float=1.0,line_search:str="armijo",
        hessian_approximation:str="BFGS",verbose:bool=False):

        super().__init__(f,grad_f,H_f,x0,verbose)

        self.epsilon = epsilon
        self.beta = beta
        self.direction = direction
        self.alpha = alpha
        self.line_search = line_search
        self.hessian_approximation = hessian_approximation

        # initialize
        if self.direction == "newton":
            assert self.H_f is not None, "you must provide the Hessian"
            if self.hessian_approximation == "exact":
                self.Bk = np.linalg.inv(self.H_f(self.xk).squeeze(axis=2))
            elif self.hessian_approximation == "BFGS":
                self.Bk = np.eye(self.ndim)
        elif self.direction == "steepest":
            self.Bk = np.eye(self.ndim)
        self.dk = - self.Bk @ self.grad_fk

    def _line_search(self,xk,dk,epsilon=0.1,beta=0.9,alpha=0.5) -> int:
        phi = lambda alpha_k: self.f(xk + alpha_k*dk)
        phi_prime = lambda alpha_k: self.grad_f(xk + alpha_k*dk).T @ dk

        n = 0
        while True:
            if n > 100:
                return alpha
                
            if (phi(alpha) <= phi(0) - epsilon*phi_prime(0)*alpha) and \
                (phi(2*alpha) >= phi(0) - 2*epsilon*phi_prime(0)*alpha):
                return alpha

            if phi(alpha) > phi(0) - epsilon*phi_prime(0)*alpha:
                alpha = beta*alpha
            elif phi(2*alpha) < phi(0) - 2*epsilon*phi_prime(0)*alpha:
                alpha = alpha / beta
            
            # # to avoid oscillating alpha
            # beta = min(beta + beta * (n/1000),1-1e-3)
            n += 1
        
    def step(self) -> Tuple[float,float]:

        norm_df = np.linalg.norm(self.grad_f(self.xk).squeeze(axis=2))

        # determine direction
        if self.direction == "newton":
            self.Hk = self.H_f(self.xk).squeeze(axis=2)

            try:
                Bk_true = np.linalg.inv(self.Hk)
            except LinAlgError:
                Bk_true = np.eye(self.ndim)
                warnings.warn("Singular true matrix, setting Bk_true to the identity matrix")

            error_BFGS = ((self.Bk - Bk_true)**2).mean(axis=None)
        elif self.direction == "steepest":
            error_BFGS = None        

        # determine step size
        if self.line_search == "armijo":
            alpha_k = self._line_search(self.xk,self.dk,self.epsilon,self.beta,self.alpha) # perform Armijo line search
        elif self.line_search == "fixed":
            alpha_k = self.alpha
        elif self.line_search == "exact":
            assert self.H_f is not None, "you must provide the Hessian"
            alpha_k = (self.grad_fk.T @ self.grad_fk) / (self.grad_fk.T @ self.Hk @ self.grad_fk)
            if self.grad_fk.T @ self.Hk @ self.grad_fk < 0:
                warnings.warn("SOSCs failed during exact line search")
        
        sk = alpha_k*self.dk

        if self.verbose:
            print("alpha_k = %f" %alpha_k)
            print("norm of gradient = ")
            print(np.linalg.norm(self.grad_fk))
            print("xk = ")
            print(self.xk)
            if self.direction == "newton":
                print("Bk = ")
                print(self.Bk)
                print("Bk_true = ")
                print(Bk_true)

        #########################################
        # UPDATES
        self.xk = self.xk + sk
        # Update Hessian approximation
        if self.hessian_approximation == "exact" and self.direction == "newton":
            self.Bk = np.linalg.inv(self.H_f(self.xk).squeeze(axis=2))
        elif self.hessian_approximation == "BFGS" and self.direction == "newton":
            yk = self.grad_f(self.xk).squeeze(axis=2) - self.grad_fk
            self.Bk = self.Bk + (((sk.T @ yk + yk.T @ self.Bk @ yk) * (sk @ sk.T)) / ((sk.T @ yk)**2)) - \
                ((self.Bk @ yk @ sk.T + sk @ yk.T @ self.Bk) / (sk.T @ yk))
        # update f and its gradients
        self.fk = self.f(self.xk)
        self.grad_fk = self.grad_f(self.xk).squeeze(axis=2)
        if self.H_f is not None: self.Hk = self.H_f(self.xk).squeeze(axis=2)
        # update direction
        self.dk = - self.Bk @ self.grad_fk
        return norm_df,error_BFGS

    def optimize(self,tol:float=1e-6,n_iterations:int=1000) -> Tuple[np.ndarray,float]:
        k = 0
        while True:
            norm_df,_ = self.step()
            if norm_df <= tol or k >= n_iterations:
                return self.xk,self.fk
            k += 1


class ConjugateDesc(Optimization):

    def __init__(self,f:Callable[[np.ndarray,], np.ndarray],grad_f:Callable[[np.ndarray,], np.ndarray],
        H_f:Callable[[np.ndarray,], np.ndarray]=None,x0:np.ndarray=np.array([[0.0,0.0]]),verbose:bool=False):

        super().__init__(f,grad_f,H_f,x0,verbose)

        # initialize
        self.gk = self.grad_fk
        self.dk = -self.gk
        self.beta_k = 0.0
        
    def step(self) -> float:

        norm_df = np.linalg.norm(self.grad_f(self.xk).squeeze(axis=2))      

        # determine step size
        alpha_k = -(self.gk.T @ self.dk) / (self.dk.T @ self.Hk @ self.dk)

        if self.verbose:
            print("alpha_k = %f" %alpha_k)
            print("norm of gradient = ")
            print(np.linalg.norm(self.grad_fk))
            print("xk = ")
            print(self.xk)

        #########################################
        # UPDATES
        self.xk = self.xk + alpha_k*self.dk

        # update direction
        self.gk = self.grad_f(self.xk).squeeze(axis=2)
        self.beta_k = (self.gk.T @ self.Hk @ self.dk) / (self.dk.T @ self.Hk @ self.dk)
        self.dk = - self.gk + self.beta_k * self.dk

        # update f and its gradients
        self.fk = self.f(self.xk)
        self.grad_fk = self.grad_f(self.xk).squeeze(axis=2)
        if self.H_f is not None: self.Hk = self.H_f(self.xk).squeeze(axis=2)

        return norm_df

    def optimize(self,tol:float=1e-6,n_iterations:int=1000) -> Tuple[np.ndarray,float]:
        k = 0
        while True:
            self.step()
            if np.linalg.norm(self.gk) <= tol or k >= n_iterations:
                return self.xk,self.fk
            k += 1