import numpy as np

def inv_guess(x):
    input = np.linspace(0.5,1,64)
    ind = min(range(len(input)), key=lambda i: abs(input[i] - x))
    guess = np.linspace(1,2,64)
    guess = guess[::-1]
    return guess[ind]

def inv(x):
    err = 10**(-8)
    xi = inv_guess(x)
    diff = 1
    count = 0
    while diff > err:
        diff = np.abs(1/x - xi)
        delta = x * xi - 1
        xi = xi*(1-delta)
        count += 1
        print(diff)
    return xi

print(inv(0.555))
print(1/0.555)