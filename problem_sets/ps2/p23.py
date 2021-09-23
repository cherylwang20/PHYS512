import numpy as np


x = np.linspace(0.5,1,100)
y = np.log2(x)



n = np.polynomial.chebyshev.chebfit(x,y,6)

print(n)