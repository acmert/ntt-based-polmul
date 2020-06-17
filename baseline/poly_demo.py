from math import log
from random import randint

from generate_prime import *
from helper import *
from ntt import *
from poly import *

# Parameter generation

# Determine n and bit-size of q, then find a q satisfying
# the condition: q = 1 (mod 2n) or q = 1 (mod n)
#
# Based on n and q, polynomial multiplication parameters

# Parameters
mod     = 2 # if 1 --> q = 1 (mod n), if 2 --> q = 1 (mod 2n)
n       = 256
q_bit   = 13

q       = 0
w       = 0
w_inv   = 0
psi     = 0
psi_inv = 0

# Generate parameters
wfound = False
while(not(wfound)):
    q = generate_large_prime(q_bit)

    # check q = 1 (mod n or 2n)
    while (not ((q % (mod*n)) == 1)):
        q = generate_large_prime(q_bit)

    # generate NTT parameters
    for i in range(2,q-1):
        wfound = isrootofunity(i,mod*n,q)
        if wfound:
            if mod == 1:
                psi    = 0
                psi_inv= 0
                w      = i
                w_inv  = modinv(w,q)
            else:
                psi    = i
                psi_inv= modinv(psi,q)
                w      = pow(psi,2,q)
                w_inv  = modinv(w,q)
            break

# Print parameters
print("Parameters (NWC)")
print("n      : {}".format(n))
print("q      : {}".format(q))
print("w      : {}".format(w))
print("w_inv  : {}".format(w_inv))
print("psi    : {}".format(psi))
print("psi_inv: {}".format(psi_inv))
print("")

# Parameters (NTRU)
m       = 3*n
mq_bit  = 14

mq      = 0
mw      = 0
mw_inv  = 0

# Generate parameters
wfound = False
while(not(wfound)):
    mq = generate_large_prime(mq_bit)

    # check q = 1 (mod n or 2n)
    while (not ((mq % m) == 1)):
        mq = generate_large_prime(mq_bit)

    # generate NTT parameters
    for i in range(2,mq-1):
        wfound = isrootofunity(i,m,mq)
        if wfound:
            mw      = i
            mw_inv  = modinv(mw,mq)
            break

# m,mq,mw,mw_inv = 192,769,4,577

# Powers of twiddle factors for NTRU (forward and inverse transform)
# Generating necessary powers of twiddle factors for NTRU on-the-fly is really hard.
# Therefore, we create table for powers of twiddle factors prior any operation
nf = [0]*(m//3) # forward

nf[0] = 0
nf[1] = m//6
nf[2] = nf[1]//2
nf[3] = (5*nf[1])//2

i = 2
while (2**i) < (m//3):
    for j in range(2**i, 2**(i+1), 2):
        nf[j]   =  nf[j//2]//2
        nf[j+1] = (nf[j//2]+(m//2))//2
    i = i + 1

ntrupowersf = nf[2:]

ntrupowersi = [] # inverse

idxs, idxe = len(ntrupowersf)-(m//6) ,len(ntrupowersf)
for i in range(int(log(m//6,2))):
    ntrupowersi = ntrupowersi + ntrupowersf[idxs:idxe]
    idxe = idxs
    idxs = idxs - ((m//12)>>i)

ntrupowersb = [0]*(m//3) # basemul

for i in range(m//6):
    ntrupowersb[2*i+0] = ntrupowersi[i]
    ntrupowersb[2*i+1] = ntrupowersi[i] + (m//2)

# print(ntrupowersf)
# print(ntrupowersb)
# print(ntrupowersi)

print("Parameters (NTRU)")
print("m      : {}".format(m))
print("mq     : {}".format(mq))
print("mw     : {}".format(mw))
print("mw_inv : {}".format(mw_inv))
print("")

#NOTE: Comment Out Reference Method for Large Parameters

# Demo
# Random A,B
A = [randint(0,q-1) for _ in range(n)]
B = [randint(0,q-1) for _ in range(n)]

# Random A,B (for ntru)
A_ntru = [randint(0,mq-1) for _ in range(m)]
B_ntru = [randint(0,mq-1) for _ in range(m)]

# Evaluator
Evaluator = Poly()

# reduce functions
pwc  = [-1]+[0]*(n-1)+[1]
nwc  =  [1]+[0]*(n-1)+[1]
ntru =  [1]+[0]*(int(m/2)-1)+[-1]+[0]*(int(m/2)-1)+[1]

print("-------- Sanity check for polynomial multiplication operations --------")
print("")

# Check reference implementations
D0 = Evaluator.SchoolbookPolMul(A,B,q)
D1 = Evaluator.SchoolbookPolMul(A_ntru,B_ntru,mq)
DR0= Evaluator.PolRed(D0,pwc,q) # reduce with x^n-1
DR1= Evaluator.PolRed(D0,nwc,q) # reduce with x^n+1
DR2= Evaluator.PolRed(D1,ntru,mq)# reduce with x^n-x^(n/2)+1
C0 = Evaluator.SchoolbookModPolMul_PWC(A,B,q)
C1 = Evaluator.SchoolbookModPolMul_NWC(A,B,q)
C2 = Evaluator.SchoolbookModPolMul_NTRU(A_ntru,B_ntru,mq)

print("SchoolbookModPolMul_PWC  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(DR0,C0)]) == 0) else "Wrong"))
print("SchoolbookModPolMul_NWC  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(DR1,C1)]) == 0) else "Wrong"))
print("SchoolbookModPolMul_NTRU --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(DR2,C2)]) == 0) else "Wrong"))
print("")

# Check NTT-based polynomial multiplication methods
N0 = Evaluator.NTTBasedPolMul(A,B,psi,psi_inv,q)
N1 = Evaluator.NTTBasedModPolMul_PWC(A,B,w,w_inv,q)
N2 = Evaluator.NTTBasedModPolMul_NWC_v1(A,B,w,w_inv,psi,psi_inv,q)
N3 = Evaluator.NTTBasedModPolMul_NWC_v2(A,B,psi,psi_inv,q)

print("NTTBasedPolMul           --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N0,D0)]) == 0) else "Wrong"))
print("NTTBasedModPolMul_PWC    --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N1,C0)]) == 0) else "Wrong"))
print("NTTBasedModPolMul_NWC_v1 --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N2,C1)]) == 0) else "Wrong"))
print("NTTBasedModPolMul_NWC_v2 --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N3,C1)]) == 0) else "Wrong"))
print("")

# Check CRT-based polynomial multiplication methods
T0 = Evaluator.CRTBasedModPolMul_PWC(A,B,w,w_inv,q)
T1 = Evaluator.CRTBasedModPolMul_NWC_FD1(A,B,psi,psi_inv,q)
T2 = Evaluator.CRTBasedModPolMul_NWC_FD2(A,B,w,w_inv,q)
T3 = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,psi,psi_inv,q,findeg=1)
T4 = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,w,w_inv,q,findeg=2)
T5 = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,w**2 % q,w_inv**2 % q,q,findeg=4)
T6 = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,w**4 % q,w_inv**4 % q,q,findeg=8)
T7 = Evaluator.CRTBasedModPolMul_NTRU_FD3(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq)
T8 = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=3)
T9 = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=6)
T10= Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=12)
T11= Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=24)
# T9 = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw**2 % q,mw_inv**2 % q,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=6)
# T10= Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw**4 % q,mw_inv**4 % q,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=12)
# T11= Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw**8 % q,mw_inv**8 % q,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=24)

print("CRTBasedModPolMul_PWC                  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T0,C0)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NWC_FD1              --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T1,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NWC_FD2              --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T2,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NWC_FDV  (findeg=1)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T3,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NWC_FDV  (findeg=2)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T4,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NWC_FDV  (findeg=4)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T5,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NWC_FDV  (findeg=8)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T6,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NTRU_FD3             --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T7,C2)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NTRU_FDV (findeg=3)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T8,C2)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NTRU_FDV (findeg=6)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T9,C2)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NTRU_FDV (findeg=12) --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T10,C2)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_NTRU_FDV (findeg=24) --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T11,C2)]) == 0) else "Wrong"))
print("")

# Check CRT-based unified polynomial multiplication methods
# ring: 0 --> NWC  (x^n + 1)
# -- findeg: 1
# -- findeg: 2
# -- findeg: 4
# --   ...
# ring: 1 --> NTRU (x^n - x^(n/2) + 1)
# -- findeg: 3
# -- findeg: 6
# -- findeg: 12
# --
ring,findeg = 0,1
R0 = Evaluator.CRTBasedModPolMul_Unified(A,B,psi,psi_inv,q,ring,findeg) # NWC - findeg=2
ring,findeg = 0,2
R1 = Evaluator.CRTBasedModPolMul_Unified(A,B,w,w_inv,q,ring,findeg) # NWC - findeg=2
ring,findeg = 0,4
R2 = Evaluator.CRTBasedModPolMul_Unified(A,B,w**2 % q,w_inv**2 % q,q,ring,findeg) # NWC - findeg=4
ring,findeg = 1,3
R3 = Evaluator.CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=3
ring,findeg = 1,6
R4 = Evaluator.CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=6
ring,findeg = 1,12
R5 = Evaluator.CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=12

print("CRTBasedModPolMul_Unified (NWC  - findeg=1)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R0,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_Unified (NWC  - findeg=2)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R1,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_Unified (NWC  - findeg=4)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R2,C1)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_Unified (NTRU - findeg=3)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R3,C2)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_Unified (NTRU - findeg=6)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R4,C2)]) == 0) else "Wrong"))
print("CRTBasedModPolMul_Unified (NTRU - findeg=12) --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R5,C2)]) == 0) else "Wrong"))

#
# NOTE: We can have extra optimization by combining N_inv with last stage of INTT
# NOTE: Later I can add PWC (x^n - 1)
# NOTE: Later I can add pure NTT/INTT operations
#
