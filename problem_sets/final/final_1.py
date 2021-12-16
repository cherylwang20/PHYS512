import numpy as np
import time

# creating the lookup table for our guess value
def inv_guess(x):
    input = np.linspace(0.5,1,64)
    ind = min(range(len(input)), key=lambda i: abs(input[i] - x)) # finding the index of the closest number without using division
    guess = np.linspace(1,2,64)
    guess = guess[::-1]
    return guess[ind] # output the approximate value from the lookup table

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
    return xi,count

def inv2(x):
    err = 10 ** (-8)
    xi = inv_guess(x)
    diff = 1
    count = 0
    while diff > err:
        diff = np.abs(1 / x - xi)
        delta = x * xi - 1
        xi = xi * (1 - delta + delta**2)
        count += 1
    return xi, count


test = np.linspace(0.5, 1, 1000)
count_test = [0]*len(test)
times = [0]*len(test)
for i in range(len(test)):
    start = time.time()
    count_test[i] = inv(test[i])[1]
    end = time.time()
    times[i] = end - start

print(f'The average time is:{np.mean(times)}')
print(np.mean(count_test))

print(inv2(0.666))