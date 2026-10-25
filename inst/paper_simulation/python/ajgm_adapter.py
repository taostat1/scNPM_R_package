
import os, sys, json, numpy as np
from AJGM import AJGM
from convergence_check import checked_ajgm
import prox_tv
# Compatibility only: AJGM supplies lists; prox_tv 3.3 requires ndarray
# for its weight comparison. Values and optimization are unchanged.
_tv1w = prox_tv.tv1w_1d
def _tv1w_array(x,w,*args,**kwargs):
    return _tv1w(np.asarray(x,dtype=float),np.asarray(w,dtype=float),*args,**kwargs)
prox_tv.tv1w_1d = _tv1w_array
root, seed, steps, edge_cut = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
for d in ["T","C"]:
    np.random.seed(seed)
    X = np.log1p(np.loadtxt(os.path.join(root,d+".csv"),delimiter=","))
    scale = X.std(axis=0,ddof=1)
    scale[scale==0] = 1.
    X = X/scale  # Do NOT center: observed zeros must remain zero for imputation.
    print("AJGM condition "+d, flush=True)
    res = checked_ajgm(AJGM, os.path.join(root,d+"_convergence.json"), X=X.copy(),Y=X.copy(),L=3,lambda1=.01,lambda2=.01,lambda3=.1,
               miu=.1,tre=.1,impute_x=True,impute_y="x",
               EMthreshold=.1,MAX_iter=steps)
    P = np.asarray(res["covinv"])
    labels = np.asarray(res["membership"]).astype(int)+1
    if P.shape != (X.shape[1],4*X.shape[1]) or not np.isfinite(P).all():
        raise RuntimeError("Invalid AJGM precision output")
    np.savetxt(os.path.join(root,d+"_labels.csv"),labels,fmt="%d",delimiter=",")
    for k in range(3):
        a = P[:,k*X.shape[1]:(k+1)*X.shape[1]]
        if np.linalg.eigvalsh((a+a.T)/2).min()<=0:
            raise RuntimeError("AJGM returned non-positive-definite precision")
        np.savetxt(os.path.join(root,d+"_P"+str(k+1)+".csv"),a,delimiter=",")
print("AJGM_DONE",flush=True)
