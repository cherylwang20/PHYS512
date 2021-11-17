# seed the pseudorandom number generator
# generate triples number through python

import sys
print(sys.maxsize)
print(2**31)


from random import seed
import numpy as np
import random
# seed random number generator

seed(np.pi)
print(random.randint(0, 2**31))



#@nb.njit
def get_rands_nb(vals):
    n=len(vals)
    for i in range(n):
        vals[i]= random.randint(0,2**31)
    return vals

def get_rands(n):
    vec=np.empty(n,dtype='int64')
    get_rands_nb(vec)
    return vec


n=200000000
vec=get_rands(n*3)
#vv=vec&(2**16-1)

vv=np.reshape(vec,[n,3])
vmax=np.max(vv,axis=1)

maxval=1e8
vv2=vv[vmax<maxval,:]

f=open('rand_points_py.txt','w')
for i in range(vv2.shape[0]):
    myline=repr(vv2[i,0])+' '+repr(vv2[i,1])+' '+ repr(vv2[i,2])+'\n'
    f.write(myline)
f.close()
