import numpy as np

def simplex(A:np.ndarray,b:np.ndarray,c:np.ndarray):

    n,m = A.shape

    # optimal solution
    x_opt = np.zeros(m)
    for j in range(m): # columns
        if c[j] == 0:
            x_opt[j] = max(A[:,j].T @ b[:,None],0)

    q = np.argmin(c); # pivot column
    j = c[q]
    if j < 0:
        A_row = b / A[:,q]
        A_row[A_row<0] = np.inf
        p = np.argmin(A_row) # pivot row
        i = A_row[p]
        if i != np.inf:
            An = A.copy(); bn = b.copy(); cn = c.copy()
            for i in range(n): # rows
                if i == p:
                    for j in range(m): # columns
                        A[i,j] = An[i,j]/An[i,q]
                        b[i] = bn[i]/An[i,q]
                else:
                    for j in range(m): # columns
                        A[i,j] = An[i,j] - (An[i,q]/An[p,q])*An[p,j]
                        c[j] = cn[j] - (cn[q]/An[p,q])*An[p,j]
                    
                    b[i] = bn[i] - (An[i,q]/An[p,q])*bn[p]

            stop = False
        else:
            print("Problem is unbounded")
            x_opt = np.inf*np.ones(m)
            stop = True
    else:
        print("solution terminated")
        stop = True

    return x_opt,stop

if __name__ == "__main__":

    A = np.array([
        [2,3,-2,-3,1,0,0,0,0],
        [-5,-2,5,2,0,1,0,0,0],
        [-2,7,2,-7,0,0,1,0,0],
        [0,1,0,-1,0,0,0,-1,0],
        [1,0,-1,0,0,0,0,0,1]], dtype="float64")
    b = np.array([10,-2,8,-5,5], dtype="float64")
    c = np.array([-1,1,1,-1,0,0,0,0,0], dtype="float64")

    cf = c.copy()
    k = 0
    stop = False

    k = 0
    cf = c
    print("==================")
    print("Simplex iterations")
    while not stop:
        x_opt,stop = simplex(A,b,c)

        if not stop:
            k+=1
            print("iteration = %i" %k)
            print("A =", A)
            print("c =",c)
            print("b =",b)

        print("x_opt =",x_opt)

        f_opt = cf[None,:] @ x_opt[:,None]
        z_opt = x_opt[0:3] - x_opt[2:5]

        print("==================")
        # input("hit ENTER to continue")