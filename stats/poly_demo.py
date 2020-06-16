from math import log
from random import randint

from generate_prime import *
from helper import *
from ntt import *
from poly import *

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

print("-------- Sanity check for polynomial multiplication operations --------")
print("")

# Check reference implementations
D0 ,D0Mul ,D0Add ,D0Sub ,D0Btf  = Evaluator.SchoolbookPolMul(A,B,q)
D1 ,D1Mul ,D1Add ,D1Sub ,D1Btf  = Evaluator.SchoolbookPolMul(A_ntru,B_ntru,mq)
DR0,DR0Mul,DR0Add,DR0Sub,DR0Btf = Evaluator.PolRed(D0,pwc,q) # reduce with x^n-1
DR1,DR1Mul,DR1Add,DR1Sub,DR1Btf = Evaluator.PolRed(D0,nwc,q) # reduce with x^n+1
DR2,DR2Mul,DR2Add,DR2Sub,DR2Btf = Evaluator.PolRed(D1,ntru,mq)# reduce with x^n-x^(n/2)+1
C0 ,C0Mul ,C0Add ,C0Sub ,C0Btf  = Evaluator.SchoolbookModPolMul_PWC(A,B,q)
C1 ,C1Mul ,C1Add ,C1Sub ,C1Btf  = Evaluator.SchoolbookModPolMul_NWC(A,B,q)
C2 ,C2Mul ,C2Add ,C2Sub ,C2Btf  = Evaluator.SchoolbookModPolMul_NTRU(A_ntru,B_ntru,mq)

print("SchoolbookModPolMul_PWC  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(DR0,C0)]) == 0) else "Wrong"))
print("-- Mul:{}".format(C0Mul))
# print("-- Add:{}".format(C0Add))
# print("-- Sub:{}".format(C0Sub))
# print("-- Btf:{}".format(C0Btf))
print("SchoolbookModPolMul_NWC  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(DR1,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(C1Mul))
# print("-- Add:{}".format(C1Add))
# print("-- Sub:{}".format(C1Sub))
# print("-- Btf:{}".format(C1Btf))
print("SchoolbookModPolMul_NTRU --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(DR2,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(C2Mul))
# print("-- Add:{}".format(C2Add))
# print("-- Sub:{}".format(C2Sub))
# print("-- Btf:{}".format(C2Btf))
print("")

# ------------------------------------------------------------------------------

# Check NTT-based polynomial multiplication methods
N0,N0Mul,N0Add,N0Sub,N0Btf = Evaluator.NTTBasedPolMul(A,B,psi,psi_inv,q)
N1,N1Mul,N1Add,N1Sub,N1Btf = Evaluator.NTTBasedModPolMul_PWC(A,B,w,w_inv,q)
N2,N2Mul,N2Add,N2Sub,N2Btf = Evaluator.NTTBasedModPolMul_NWC_v1(A,B,w,w_inv,psi,psi_inv,q)
N3,N3Mul,N3Add,N3Sub,N3Btf = Evaluator.NTTBasedModPolMul_NWC_v2(A,B,psi,psi_inv,q)

print("NTTBasedPolMul           --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N0,D0)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N0Mul))
# print("-- Add:{}".format(N0Add))
# print("-- Sub:{}".format(N0Sub))
# print("-- Btf:{}".format(N0Btf))
print("NTTBasedModPolMul_PWC    --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N1,C0)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N1Mul))
# print("-- Add:{}".format(N1Add))
# print("-- Sub:{}".format(N1Sub))
# print("-- Btf:{}".format(N1Btf))
print("NTTBasedModPolMul_NWC_v1 --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N2,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N2Mul))
# print("-- Add:{}".format(N2Add))
# print("-- Sub:{}".format(N2Sub))
# print("-- Btf:{}".format(N2Btf))
print("NTTBasedModPolMul_NWC_v2 --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(N3,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N3Mul))
# print("-- Add:{}".format(N3Add))
# print("-- Sub:{}".format(N3Sub))
# print("-- Btf:{}".format(N3Btf))
print("")

# ------------------------------------------------------------------------------

# Check CRT-based polynomial multiplication methods
T0 ,T0Mul ,T0Add ,T0Sub ,T0Btf  = Evaluator.CRTBasedModPolMul_PWC(A,B,w,w_inv,q)
T1 ,T1Mul ,T1Add ,T1Sub ,T1Btf  = Evaluator.CRTBasedModPolMul_NWC_FD1(A,B,psi,psi_inv,q)
T2 ,T2Mul ,T2Add ,T2Sub ,T2Btf  = Evaluator.CRTBasedModPolMul_NWC_FD2(A,B,w,w_inv,q)
T3 ,T3Mul ,T3Add ,T3Sub ,T3Btf  = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,psi,psi_inv,q,findeg=1)
T4 ,T4Mul ,T4Add ,T4Sub ,T4Btf  = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,w,w_inv,q,findeg=2)
T5 ,T5Mul ,T5Add ,T5Sub ,T5Btf  = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,w**2 % q,w_inv**2 % q,q,findeg=4)
T6 ,T6Mul ,T6Add ,T6Sub ,T6Btf  = Evaluator.CRTBasedModPolMul_NWC_FDV(A,B,w**4 % q,w_inv**4 % q,q,findeg=8)
T7 ,T7Mul ,T7Add ,T7Sub ,T7Btf  = Evaluator.CRTBasedModPolMul_NTRU_FD3(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq)
T8 ,T8Mul ,T8Add ,T8Sub ,T8Btf  = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=3)
T9 ,T9Mul ,T9Add ,T9Sub ,T9Btf  = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=6)
T10,T10Mul,T10Add,T10Sub,T10Btf = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=12)
T11,T11Mul,T11Add,T11Sub,T11Btf = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw,mw_inv,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=24)
# T9 = Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw**2 % q,mw_inv**2 % q,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=6)
# T10= Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw**4 % q,mw_inv**4 % q,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=12)
# T11= Evaluator.CRTBasedModPolMul_NTRU_FDV(A_ntru,B_ntru,mw**8 % q,mw_inv**8 % q,ntrupowersf,ntrupowersb,ntrupowersi,mq,findeg=24)

print("CRTBasedModPolMul_PWC                  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T0,C0)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T0Mul))
# print("-- Add:{}".format(T0Add))
# print("-- Sub:{}".format(T0Sub))
# print("-- Btf:{}".format(T0Btf))
print("CRTBasedModPolMul_NWC_FD1              --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T1,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T1Mul))
# print("-- Add:{}".format(T1Add))
# print("-- Sub:{}".format(T1Sub))
# print("-- Btf:{}".format(T1Btf))
print("CRTBasedModPolMul_NWC_FD2              --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T2,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T2Mul))
# print("-- Add:{}".format(T2Add))
# print("-- Sub:{}".format(T2Sub))
# print("-- Btf:{}".format(T2Btf))
print("CRTBasedModPolMul_NWC_FDV  (findeg=1)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T3,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T3Mul))
# print("-- Add:{}".format(T3Add))
# print("-- Sub:{}".format(T3Sub))
# print("-- Btf:{}".format(T3Btf))
print("CRTBasedModPolMul_NWC_FDV  (findeg=2)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T4,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T4Mul))
# print("-- Add:{}".format(T4Add))
# print("-- Sub:{}".format(T4Sub))
# print("-- Btf:{}".format(T4Btf))
print("CRTBasedModPolMul_NWC_FDV  (findeg=4)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T5,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T5Mul))
# print("-- Add:{}".format(T5Add))
# print("-- Sub:{}".format(T5Sub))
# print("-- Btf:{}".format(T5Btf))
print("CRTBasedModPolMul_NWC_FDV  (findeg=8)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T6,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T6Mul))
# print("-- Add:{}".format(T6Add))
# print("-- Sub:{}".format(T6Sub))
# print("-- Btf:{}".format(T6Btf))
print("CRTBasedModPolMul_NTRU_FD3             --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T7,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T7Mul))
# print("-- Add:{}".format(T7Add))
# print("-- Sub:{}".format(T7Sub))
# print("-- Btf:{}".format(T7Btf))
print("CRTBasedModPolMul_NTRU_FDV (findeg=3)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T8,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T8Mul))
# print("-- Add:{}".format(T8Add))
# print("-- Sub:{}".format(T8Sub))
# print("-- Btf:{}".format(T8Btf))
print("CRTBasedModPolMul_NTRU_FDV (findeg=6)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T9,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T9Mul))
# print("-- Add:{}".format(T9Add))
# print("-- Sub:{}".format(T9Sub))
# print("-- Btf:{}".format(T9Btf))
print("CRTBasedModPolMul_NTRU_FDV (findeg=12) --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T10,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T10Mul))
# print("-- Add:{}".format(T10Add))
# print("-- Sub:{}".format(T10Sub))
# print("-- Btf:{}".format(T10Btf))
print("CRTBasedModPolMul_NTRU_FDV (findeg=24) --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(T11,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(T11Mul))
# print("-- Add:{}".format(T11Add))
# print("-- Sub:{}".format(T11Sub))
# print("-- Btf:{}".format(T11Btf))
print("")

# ------------------------------------------------------------------------------
# Parallelism?
# NTRU is best utilized with 3*power-of-2 PEs (can be parallelised up to m//2)
# NWC is best utilized with power-of-2 PEs (can be parallelised up to n//2)
# -- An unified architecture is best utilized with power-of-2 PEs
# All systems needs 2*PE BRAMs for in-place computation

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
R0 ,R0Mul ,R0Add ,R0Sub ,R0Btf = Evaluator.CRTBasedModPolMul_Unified(A,B,psi,psi_inv,q,ring,findeg) # NWC - findeg=2
ring,findeg = 0,2
R1 ,R1Mul ,R1Add ,R1Sub ,R1Btf = Evaluator.CRTBasedModPolMul_Unified(A,B,w,w_inv,q,ring,findeg) # NWC - findeg=2
ring,findeg = 0,4
R2 ,R2Mul ,R2Add ,R2Sub ,R2Btf = Evaluator.CRTBasedModPolMul_Unified(A,B,w**2 % q,w_inv**2 % q,q,ring,findeg) # NWC - findeg=4
ring,findeg = 1,3
R3 ,R3Mul ,R3Add ,R3Sub ,R3Btf = Evaluator.CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=3
ring,findeg = 1,6
R4 ,R4Mul ,R4Add ,R4Sub ,R4Btf = Evaluator.CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=6
ring,findeg = 1,12
R5 ,R5Mul ,R5Add ,R5Sub ,R5Btf = Evaluator.CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=12

print("CRTBasedModPolMul_Unified (NWC  - findeg=1)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R0,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R0Mul))
# print("-- Add:{}".format(R0Add))
# print("-- Sub:{}".format(R0Sub))
# print("-- Btf:{}".format(R0Btf))
print("CRTBasedModPolMul_Unified (NWC  - findeg=2)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R1,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R1Mul))
# print("-- Add:{}".format(R1Add))
# print("-- Sub:{}".format(R1Sub))
# print("-- Btf:{}".format(R1Btf))
print("CRTBasedModPolMul_Unified (NWC  - findeg=4)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R2,C1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R2Mul))
# print("-- Add:{}".format(R2Add))
# print("-- Sub:{}".format(R2Sub))
# print("-- Btf:{}".format(R2Btf))
print("CRTBasedModPolMul_Unified (NTRU - findeg=3)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R3,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R3Mul))
# print("-- Add:{}".format(R3Add))
# print("-- Sub:{}".format(R3Sub))
# print("-- Btf:{}".format(R3Btf))
print("CRTBasedModPolMul_Unified (NTRU - findeg=6)  --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R4,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R4Mul))
# print("-- Add:{}".format(R4Add))
# print("-- Sub:{}".format(R4Sub))
# print("-- Btf:{}".format(R4Btf))
print("CRTBasedModPolMul_Unified (NTRU - findeg=12) --> " + ("Correct" if(sum([abs(x-y) for x,y in zip(R5,C2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R5Mul))
# print("-- Add:{}".format(R5Add))
# print("-- Sub:{}".format(R5Sub))
# print("-- Btf:{}".format(R5Btf))
print("")

# ------------------------------------------------------------------------------

#
# NOTE: We can have extra optimization by combining N_inv with last stage of INTT
# NOTE: Later I can add PWC (x^n - 1)
# NOTE: Later I can add pure NTT/INTT operations
#
