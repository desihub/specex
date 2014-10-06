
import numpy,scipy,scipy.linalg

def cholesky_solve(A,B,overwrite=False,lower=False) :
    posv, = scipy.linalg.get_lapack_funcs(('posv',), (A,))
    L,X,status=posv(A,B,lower=lower,overwrite_a=overwrite)
    
    if status :
        print "error cholesky_solve_and_invert error dposv status=",status
        raise Exception("cholesky_solve_and_invert error dposv status=%d"%status)
    
    return X

def cholesky_solve_and_invert(A,B,overwrite=False,lower=False) :
    posv, = scipy.linalg.get_lapack_funcs(('posv',), (A,))
    L,X,status=posv(A,B,lower=lower,overwrite_a=overwrite)
    
    if status :
        print "error cholesky_solve_and_invert error dposv status=",status
        raise Exception("cholesky_solve_and_invert error dposv status=%d"%status)
    
    potri, = scipy.linalg.get_lapack_funcs(('potri',), (L,))
    inv,status=potri(L,lower=(not lower)) # 'not lower' is not a mistake, there is a BUG in scipy!!!!   
        
    if status :
        print "error cholesky_solve_and_invert error dpotri status=",status
        raise Exception("cholesky_solve_and_invert error dpotri status=%d"%status)

    #symmetrize Ai
    if True :
        for i in range(inv.shape[0]) :
            for j in range(i) :
                if  not lower :
                    inv[i,j]=inv[j,i]
                else :
                    inv[j,i]=inv[i,j]
    return X,inv

