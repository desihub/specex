
import numpy,scipy,scipy.linalg

def cholesky_solve_and_invert(A,B,overwrite=False,lower=False) :
    a=scipy.linalg.asarray(A)
    b=scipy.linalg.asarray(B)
    posv, = scipy.linalg.get_lapack_funcs(('posv',), (a,))
    L,X,status=posv(a,b,lower=lower,overwrite_a=overwrite)
    
    if status :
        raise Error("cholesky_solve_and_invert error dposv status=",status)
    
    l=scipy.linalg.asarray(L)
    potri, = scipy.linalg.get_lapack_funcs(('potri',), (l,))
    inv,status=potri(l,lower=(not lower)) # this is not a mistake, there is a BUG in scipy!!!!   
        
    if status :
        raise Error("cholesky_solve_and_invert error dpotri status=",status)

    #symmetrize Ai
    if True :
        for i in range(inv.shape[0]) :
            for j in range(i) :
                if  not lower :
                    inv[i,j]=inv[j,i]
                else :
                    inv[j,i]=inv[i,j]
    return X,inv

