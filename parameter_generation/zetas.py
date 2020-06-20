
def gen_powers(N, q, zeta):
    
    N_ = N
    powers = [0]*N_    
    powers[0] = 1
    powers[1] = N_//2
    i = 1
    while 2**i < N_:
        for j in range(2**i, 2**(i+1), 2):
            powers[j] = powers[j//2]//2
            powers[j+1] = (powers[j//2]+N_)//2
        i = i + 1
    
    return powers

def gen_twiddles(N2, q, zeta, powers):
    twiddle_cnt = len(powers)
    twiddles = [0]*twiddle_cnt
    inv_twiddles = [0]*twiddle_cnt 
    tmp = [0]*(N2)
    for i in range(N2):
        tmp[i] = pow(zeta, i, q)
    for i in range(twiddle_cnt):
        twiddles[i] = tmp[powers[i]]
        inv_twiddles[i] = -(tmp[powers[twiddle_cnt-1-i]])%q    
    return twiddles, inv_twiddles            



q = 132120577
n = 1024
zeta = 73993  # psi(root)

powers = gen_powers(n, q, zeta)
zetas, inv_zetas = gen_twiddles(n, q, zeta, powers)

print(zetas)
print("--------------------")
print(inv_zetas)
