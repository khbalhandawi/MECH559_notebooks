import numpy as np
from abc import ABC, abstractmethod
from scipy.linalg import pinv
from numpy.linalg import norm

def scaling(X, lb, ub, type):
    if type == 1:
        return (X - lb) / (ub - lb)
    elif type == 2:
        return X * (ub - lb) + lb

class AbstractSurrogate(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def train(self):
        pass

    @abstractmethod
    def predict(self, Z):
        pass

class LS(AbstractSurrogate):
    def __init__(self, X, Y, lb, ub, r, d, scale):
        self.X = X
        self.Y = Y
        self.lb = lb
        self.ub = ub
        self.r = r
        self.d = d
        self.scale = scale
        self.n_centers, self.dim_i = X.shape
        self.dim_o = Y.shape[1]
        self.isempty = True
        self.B = np.zeros((self.n_centers, self.dim_i*self.d+1))
        self.W = np.zeros((self.dim_i*self.d+1, self.dim_o))
        self.J = np.eye(self.dim_i*self.d+1) * r
        self.J[0,0] = 0

    def Polynomial(self, d, ζ):
        n = len(ζ)
        b = np.ones((d*n+1))
        for i in range(0, d):
            b[1+i*n:(i+1)*n+1] = ζ**(i+1)
        return b

    def basis(self, d, Z):
        return np.array([self.Polynomial(d, z) for z in Z])

    def train(self):
        self.B = self.basis(self.d, self.X)
        self.W = pinv(self.B.T @ self.B + self.J) @ self.B.T @ self.Y
        self.isempty = False

    def predict(self, Z):
        if Z.shape[1] != self.X.shape[1]:
            raise ValueError("Dimensions of prediction and training sites are not equal")

        if self.isempty:
            raise ValueError("Model is not trained!")

        if self.scale:
            Z = scaling(Z, self.lb, self.ub, 1)

        return self.basis(self.d, Z) @ self.W

class RBF(AbstractSurrogate):
    def __init__(self, X, Y, lb, ub, r, λ, kernel, scale):
        self.X = X
        self.Y = Y
        self.lb = lb
        self.ub = ub
        self.scale = scale
        self.r = r
        self.λ = λ
        self.kernel = kernel
        self.n_centers, self.dim_i = X.shape
        self.dim_o = Y.shape[1]
        self.isempty = True
        self.B = np.zeros((self.n_centers, self.n_centers))
        self.W = np.zeros((self.n_centers, self.dim_o))
        self.J = np.eye(self.n_centers) * r
        self.J[0,0] = 0

    def Gaussian(self, λ, ζ, X):
        return np.exp(-λ * norm(ζ - X,axis=1))[:,None]

    def basis(self, λ, Z):
        return np.hstack([self.Gaussian(λ, z, self.X) for z in Z])

    def train(self):
        self.B = self.basis(self.λ, self.X)
        self.W = pinv(self.B.T @ self.B + self.J) @ self.B.T @ self.Y
        self.isempty = False

    def predict(self, Z):
        if Z.shape[1] != self.X.shape[1]:
            raise ValueError("Dimensions of prediction and training sites are not equal")

        if self.isempty:
            raise ValueError("Model is not trained!")

        if self.scale:
            Z = scaling(Z, self.lb, self.ub, 1)

        return self.basis(self.λ, Z).T @ self.W

class Surrogate:
    def __init__(self, X, Y, type, lb=None, ub=None, r=1e-6, λ=15.0, d=1, kernel=None, scale=True, name=""):
        if X.shape[0] != Y.shape[0]:
            raise ValueError("Number of training inputs does not match number of training outputs")

        if lb is None:
            lb = np.min(X, axis=0)

        if ub is None:
            ub = np.max(X, axis=0)

        if scale:
            X = scaling(X, lb, ub, 1)

        self.name = name
        self.type = type

        if type == "RBF":
            self.model = RBF(X, Y, lb, ub, r, λ, kernel, scale)
        elif type == "LS":
            self.model = LS(X, Y, lb, ub, r, d, scale)

    def train(self):
        self.model.train()

    def predict(self, Z):
        return self.model.predict(Z)
