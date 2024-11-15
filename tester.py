from fft import fft
from numpy.fft import fft as fft_base
from random import uniform

LB = -10**5
UB = 10**5
A_SIZE = 120_960

EPS = 5e-7 # tolerance value

print(f'Testing implemention with tolerance {EPS}')

for _ in range(30):
    A = [(uniform(LB,UB)+uniform(LB,UB)*1j) for _ in range(A_SIZE)]

    result = fft(A)
    result2 = fft_base(A)

    max_d = -float('inf')

    for a,b in zip(result,result2):
        d = abs(a-b)
        if d >= EPS:
            print(f'Failed with d: {d}')
            exit()
        max_d = max(d, max_d)

    print(f'Passed, highest d computed: {d}')