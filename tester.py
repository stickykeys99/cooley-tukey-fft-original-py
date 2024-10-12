from fft import fft
from numpy.fft import fft as fft_base
from random import uniform

LB = -10**5
UB = 10**5
A_SIZE = 360

A = [(uniform(LB,UB)+uniform(LB,UB)*1j) for _ in range(A_SIZE)]

result = fft(A)
result2 = fft_base(A)

EPS = 1e-7 # tolerance value
for a,b in zip(result,result2):
    d = abs(a-b)
    if d >= EPS:
        print('Tolerance value check failed')
        exit()

print(f'The two results are within {EPS}')