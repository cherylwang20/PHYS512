import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0.5,1,100)
y = np.log2(x)


n = np.polynomial.chebyshev.chebfit(x,y,10)
che = np.polynomial.chebyshev.chebval(x,n)

print(np.std(abs(y-che)))


# print(che)