from math import log,ceil
from random import randint

from generate_prime import *
from helper import *

# ----------------------------------------------------------

# make onle 1 of options below
DEBUG_FNTT  = 0
DEBUG_INTT  = 0
DEBUG_NWC1  = 1
DEBUG_NWC2  = 0
DEBUG_NTRU3 = 0

# Display memory patter on console
DEBUG_DISP_FNTT  = DEBUG_FNTT # Forward NTT
DEBUG_DISP_INTT  = DEBUG_INTT # Inverse NTT
DEBUG_DISP_NWC1  = DEBUG_NWC1  # Polynomial Multiplication (NWC, findeg = 1)
DEBUG_DISP_NWC2  = DEBUG_NWC2  # Polynomial Multiplication (NWC, findeg = 2)
DEBUG_DISP_NTRU3 = DEBUG_NTRU3 # Polynomial Multiplication (NTRU, findeg = 3)

# Write Memory Pattern (data and twiddle) into TXT
DEBUG_MEMP_FNTT  = DEBUG_FNTT # Forward NTT
DEBUG_MEMP_INTT  = DEBUG_INTT  # Inverse NTT
DEBUG_MEMP_NWC1  = DEBUG_NWC1  # Polynomial Multiplication (NWC, findeg = 1)
DEBUG_MEMP_NWC2  = DEBUG_NWC2  # Polynomial Multiplication (NWC, findeg = 2)
DEBUG_MEMP_NTRU3 = DEBUG_NTRU3 # Polynomial Multiplication (NTRU, findeg = 3)

# Test case generator
DEBUG_TEST_FNTT  = DEBUG_FNTT # Forward NTT
DEBUG_TEST_INTT  = DEBUG_INTT  # Inverse NTT
DEBUG_TEST_NWC1  = DEBUG_NWC1  # Polynomial Multiplication (NWC, findeg = 1)
DEBUG_TEST_NWC2  = DEBUG_NWC2  # Polynomial Multiplication (NWC, findeg = 2)
DEBUG_TEST_NTRU3 = DEBUG_NTRU3 # Polynomial Multiplication (NTRU, findeg = 3)

WLMONT           = 1 # 0: regular number representation is used / 1: word-level montgomery representation is used

if DEBUG_TEST_FNTT:
    PRM_TXT  = open("FNTT_PARAM.txt","w")
    # Parameters
    # --- n  (ring size)
    # --- q  (modulus)
    # --- K  (bit size of modulus)
    # --- PE (number of PE)
    DIN_TXT  = open("FNTT_DIN.txt","w")
    # Input polynomial
    # --- Stored in order: 0, 1, 2, ..., n-1
    DOUT_TXT = open("FNTT_DOUT.txt","w")
    # Output polynomial
    # --- Stored in bit-reversed order: br(0), br(1), br(2), ..., br(n-1)
    W_TXT    = open("FNTT_W.txt","w")
    # Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
elif DEBUG_TEST_INTT:
    PRM_TXT  = open("INTT_PARAM.txt","w")
    # Parameters
    # --- n  (ring size)
    # --- q  (modulus)
    # --- K  (bit size of modulus)
    # --- PE (number of PE)
    # --- n_inv (n^-1 mod q)
    DIN_TXT  = open("INTT_DIN.txt","w")
    # Input polynomial
    # --- Stored in order: 0, 1, 2, ..., n-1
    DOUT_TXT = open("INTT_DOUT.txt","w")
    # Output polynomial
    # --- Stored in bit-reversed order: br(0), br(1), br(2), ..., br(n-1)
    WINV_TXT = open("INTT_WINV.txt","w")
    # Inverse Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
elif DEBUG_TEST_NWC1:
    PRM_TXT  = open("NWC1_PARAM.txt","w")
    # Parameters
    # --- n  (ring size)
    # --- q  (modulus)
    # --- K  (bit size of modulus)
    # --- PE (number of PE)
    # --- n_inv (n^-1 mod q)
    DIN0_TXT = open("NWC1_DIN0.txt","w")
    DIN1_TXT = open("NWC1_DIN1.txt","w")
    # Input polynomials
    # --- Stored in order: 0, 1, 2, ..., n-1
    DIN0_MFNTT_TXT = open("NWC1_DIN0_MFNTT.txt","w")
    DIN1_MFNTT_TXT = open("NWC1_DIN1_MFNTT.txt","w")
    # Output polynomials (after merged FNTT)
    # --- Stored in bit-reversed order
    DOUT_MINTT_TXT = open("NWC1_DOUT_MINTT.txt","w")
    # Input polynomial (after coefficient-wise multiplication)
    # --- Stored in order: 0, 1, 2, ..., n-1
    DOUT_TXT = open("NWC1_DOUT.txt","w")
    # Output polynomial
    # --- Stored in order: 0, 1, 2, ..., n-1
    W_TXT    = open("NWC1_W.txt","w")     # in this case, it is actually psi
    # Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
    WINV_TXT = open("NWC1_WINV.txt","w")  # in this case, it is actually psi_inv
    # Inverse Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
elif DEBUG_TEST_NWC2:
    PRM_TXT  = open("NWC2_PARAM.txt","w")
    # Parameters
    # --- n  (ring size)
    # --- q  (modulus)
    # --- K  (bit size of modulus)
    # --- PE (number of PE)
    # --- n_inv (n^-1 mod q)
    DIN0_TXT = open("NWC2_DIN0.txt","w")
    DIN1_TXT = open("NWC2_DIN1.txt","w")
    # Input polynomials
    # --- Stored in order: 0, 1, 2, ..., n-1
    DOUT_TXT = open("NWC2_DOUT.txt","w")
    # Output polynomial
    # --- Stored in order: 0, 1, 2, ..., n-1
    W_TXT    = open("NWC2_W.txt","w")
    # Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
    WP_TXT   = open("NWC2_WP.txt","w")
    # Twiddle factors
    # --- For degree-2 polynomial multiplications after NTT operations
    WINV_TXT = open("NWC2_WINV.txt","w")
    # Inverse Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
elif DEBUG_TEST_NTRU3:
    PRM_TXT  = open("NTRU3_PARAM.txt","w")
    # Parameters
    # --- n  (ring size)
    # --- q  (modulus)
    # --- K  (bit size of modulus)
    # --- PE (number of PE)
    # --- W_pre  (power of W used in initial reduction)
    # --- W_post (power of W used in final reconstruction)
    # --- n_inv (n^-1 mod q)
    DIN0_TXT = open("NTRU3_DIN0.txt","w")
    DIN1_TXT = open("NTRU3_DIN1.txt","w")
    # Input polynomials
    # --- Stored in order: 0, 1, 2, ..., n-1
    DOUT_TXT = open("NTRU3_DOUT.txt","w")
    # Output polynomial
    # --- Stored in order: 0, 1, 2, ..., n-1
    W_TXT    = open("NTRU3_W.txt","w")
    # Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt
    WP_TXT   = open("NTRU3_WP.txt","w")
    # Twiddle factors
    # --- For degree-3 polynomial multiplications after NTT operations
    WINV_TXT = open("NTRU3_WINV.txt","w")
    # Inverse Twiddle factors
    # --- Stored as shown in FNTT/INTT_tw_N_PE.txt

# ----------------------------------------------------------

# Parallelism
PE_NUMBER = 16
PE = 2*PE_NUMBER

# ----------------------------------------------------------

# Parameter generation

# Determine n and bit-size of q, then find a q satisfying
# the condition: q = 1 (mod 2n) or q = 1 (mod n)
#
# Based on n and q, polynomial multiplication parameters

# Parameters
mod     = 2 # if 1 --> q = 1 (mod n), if 2 --> q = 1 (mod 2n)
n       = 64
q_bit   = 12

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

fd1_ws = int(log(n,2))+1 # word-size (for WL Montgomery)
fd1_L  = int(ceil(q_bit/fd1_ws)) # iteration (for WL Montgomery)
fd1_R  = 2**(fd1_ws*fd1_L)

fd2_ws = int(log(n,2))   # word-size (for WL Montgomery)
fd2_L  = int(ceil(q_bit/fd2_ws)) # iteration (for WL Montgomery)
fd2_R  = 2**(fd2_ws*fd2_L)

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
mq_bit  = 10

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

fd3_ws = int(log(n,2))   # word-size (for WL Montgomery)
fd3_L  = int(ceil(mq_bit/fd3_ws)) # iteration (for WL Montgomery)
fd3_R  = 2**(fd3_ws*fd3_L)

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

print("Parameters (HW implementation)")
print("n         : {}".format(n))
print("K=log(q)  : {}".format(q_bit))
print("m         : {}".format(m))
print("K=log(mq) : {}".format(mq_bit))
print("PE number : {}".format(PE_NUMBER))
print("")

# ----------------------------------------------------------

# From paper: NTTU: An Area-Efficient Low-POwer NTT-Uncoupled Architecture for NTT-Based Multiplication
# Iterative Radix-2 Decimation-in-Time (DIT) (CT) NTT - NR
# A: input polynomial (standard order)
# W: twiddle factor
# q: modulus
# B: output polynomial (bit-reversed order)
def Radix2_DIT_Iterative_NTT_NR(A,W,q):
    N = len(A)
    B = [_ for _ in A]

    # ---------------------------------
    v = int(math.log(N, 2))
    m = N//PE

    BRAM = []
    BRTW = []

    for i in range(PE):
        BRAM.append([])
        for j in range(v):
            BRAM[i].append([])
            for k in range(m):
                BRAM[i][j].append([])

    for i in range(PE//2):
        BRTW.append([])
        for j in range(v):
            BRTW[i].append([])
            for k in range(m):
                BRTW[i][j].append([])

    bram_counter = 0
    # ---------------------------------

    for s in range(int(log(N,2)),0,-1):
        m = 2**s
        for k in range(int(N/m)):
            TW = pow(W,intReverse(k,int(log(N,2))-s)*int(m/2),q)
            for j in range(int(m/2)):
                u = B[k*m+j]
                t = (TW*B[k*m+j+int(m/2)]) % q

                B[k*m+j]          = (u+t) % q
                B[k*m+j+int(m/2)] = (u-t) % q

                if DEBUG_DISP_FNTT: print("W: "+str(intReverse(k,int(log(N,2))-s)*int(m/2)).ljust(5)+" A0: "+str(k*m+j).ljust(5)+" A1: "+str(k*m+j+int(m/2)).ljust(5))


                # ---------------------------------
                BRAM[(2*(bram_counter >> 0) & (PE-1))+0][int(log(N,2))-s][(bram_counter & ((N//2)-1)) // (PE//2)] = k*m+j
                BRAM[(2*(bram_counter >> 0) & (PE-1))+1][int(log(N,2))-s][(bram_counter & ((N//2)-1)) // (PE//2)] = k*m+j+int(m/2)

                BRTW[bram_counter & ((PE//2)-1)][int(log(N,2))-s][(bram_counter & ((N//2)-1)) // (PE//2)] = intReverse(k,int(log(N,2))-s)*int(m/2)

                bram_counter = bram_counter + 1
                # ---------------------------------

    return B,BRAM,BRTW

# Iterative Radix-2 Decimation-in-Frequency (DIF) (GS) NTT - RN
# A: input polynomial (reversed order)
# W: twiddle factor
# q: modulus
# B: output polynomial (bit-standard order)
def Radix2_DIF_Iterative_NTT_RN(A,W,q):
    N = len(A)
    B = [_ for _ in A]

    # ---------------------------------
    v = int(math.log(N, 2))
    m = N//PE

    BRAM = []
    BRTW = []

    for i in range(PE):
        BRAM.append([])
        for j in range(v):
            BRAM[i].append([])
            for k in range(m):
                BRAM[i][j].append([])

    for i in range(PE//2):
        BRTW.append([])
        for j in range(v):
            BRTW[i].append([])
            for k in range(m):
                BRTW[i][j].append([])

    bram_counter = 0
    # ---------------------------------

    m = 1
    v = N
    d = 1

    while v>1:
        for jf in range(m):
            j = jf
            jt = 0
            while j<(N-1):
                # bit-reversing jt
                TW = pow(W,intReverse(jt,int(log(N>>1,2))),q)

                temp = B[j]

                B[j]   = (temp + B[j+d]) % q
                B[j+d] = (temp - B[j+d])*TW % q

                if DEBUG_DISP_INTT: print("W: "+str(intReverse(jt,int(log(N>>1,2)))).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+d).ljust(5))

                # ---------------------------------
                BRAM[(2*(bram_counter >> 0) & (PE-1))+0][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = j
                BRAM[(2*(bram_counter >> 0) & (PE-1))+1][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = j+d

                BRTW[bram_counter & ((PE//2)-1)][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = intReverse(jt,int(log(N>>1,2)))

                bram_counter = bram_counter + 1
                # ---------------------------------

                jt = jt+1
                j = j + 2*m

        m = 2*m
        v = int(v/2)
        d = 2*d

    return B,BRAM,BRTW

def Radix2_DIF_Iterative_INTT_RN(A,W_inv,q):
    N_inv = modinv(len(A),q)
    B,BRAM,BRTW = Radix2_DIF_Iterative_NTT_RN(A,W_inv,q)
    B = [(x*N_inv) % q for x in B]
    return B,BRAM,BRTW

# Multiplies two "deg" degree polynomial in x^"deg"-w^k where k is some power
# A,B: input polynomials
# wk: w^k
# deg: degree
# q: coefficient modulus
# C: output polynomial
def PolWiseMult(A,B,wk,deg,q):
    C = [0] * (2 * deg)
    D = [0] * deg

    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % q

    for i in range(len(A)):
        D[i] = (C[i] + C[i + len(A)]*wk) % q
    return D

# CRT-based unified structure for NWC and NTRU
# -- ring  : 0->NWC, 1->NTRU
# -- findeg: 1,2,4,... for NWC and 3,6,12,... for NTRU
# A: input polynomial (standard order)
# W: twiddle factor
# q: modulus
# B: output polynomial (bit-reversed order)
def CRT_Iterative_Unified_NR(A,w,q,ring,findeg,powers):
    N = len(A)
    B = [_ for _ in A]

    # ---------------------------------
    if ring == 0:
        v = int(math.log(N//findeg, 2))
        m = N//PE

        BRAM = []
        BRTW = []

        for i in range(PE):
            BRAM.append([])
            for j in range(v):
                BRAM[i].append([])
                for k in range(m):
                    BRAM[i][j].append([])

        for i in range(PE//2):
            BRTW.append([])
            for j in range(v):
                BRTW[i].append([])
                for k in range(m):
                    BRTW[i][j].append([])

        bram_counter = 0
    else:
        v = int(math.log(N//findeg, 2))-1
        m = N//PE

        BRAM = []
        BRTW = []

        for i in range(PE):
            BRAM.append([])
            for j in range(v):
                BRAM[i].append([])
                for k in range(m):
                    BRAM[i][j].append([])

        for i in range(PE//2):
            BRTW.append([])
            for j in range(v):
                BRTW[i].append([])
                for k in range(m):
                    BRTW[i][j].append([])

        bram_counter = 0
        ntru_counter = 0
    # ---------------------------------

    if ring == 0:
        # NWC
        k = 1
        lena = (N//2)
        v = int(log(N//findeg,2))
    else:
        # NTRU
        k = 0
        lena = (N//4)

    while lena >= findeg:
        start = 0
        while start < N:
            if ring == 0:
                # NWC
                W_pow = intReverse(k,v)
            else:
                # NTRU
                W_pow = powers[k]
                # W_pow = (powers[k] // (findeg//3))

            W = pow(w,W_pow,q)
            k = k+1
            j = start
            while(j < (start + lena)):
                t = (W * B[j+lena]) % q

                B[j+lena] = (B[j] - t) % q
                B[j     ] = (B[j] + t) % q

                if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1): print("W: "+str(W_pow).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+lena).ljust(5))
                if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2): print("W: "+str(W_pow).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+lena).ljust(5))
                if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3): print("W: "+str(W_pow).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+lena).ljust(5))

                # ---------------------------------
                if ring == 0:
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+0][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = j
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+1][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = j+lena

                    BRTW[bram_counter & ((PE//2)-1)][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = W_pow

                    bram_counter = bram_counter + 1
                else:
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+0][bram_counter // (N//6)][3*((bram_counter & ((N//6)-1)) // (PE//2))+ntru_counter] = j
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+1][bram_counter // (N//6)][3*((bram_counter & ((N//6)-1)) // (PE//2))+ntru_counter] = j+lena

                    BRTW[bram_counter & ((PE//2)-1)][bram_counter // (N//6)][3*((bram_counter & ((N//6)-1)) // (PE//2))+ntru_counter] = W_pow

                    bram_counter = (bram_counter + 1) if (ntru_counter == 2) else bram_counter
                    ntru_counter = 0 if (ntru_counter == 2) else (ntru_counter + 1)
                # ---------------------------------

                j = j+1
            start = j + lena
        lena = (lena//2)

    return B,BRAM,BRTW

# ICRT-based unified structure for NWC and NTRU
# -- ring  : 0->NWC, 1->NTRU
# -- findeg: 1,2,4,... for NWC and 3,6,12,... for NTRU
# A: input polynomial (bit-reversed order)
# W: twiddle factor
# q: modulus
# B: output polynomial (standard order)
def ICRT_Iterative_Unified_RN(A,w,q,ring,findeg,powers):
    N = len(A)
    B = [_ for _ in A]

    # ---------------------------------
    if ring == 0:
        v = int(math.log(N//findeg, 2))
        m = N//PE

        BRAM = []
        BRTW = []

        for i in range(PE):
            BRAM.append([])
            for j in range(v):
                BRAM[i].append([])
                for k in range(m):
                    BRAM[i][j].append([])

        for i in range(PE//2):
            BRTW.append([])
            for j in range(v):
                BRTW[i].append([])
                for k in range(m):
                    BRTW[i][j].append([])

        bram_counter = 0
    else:
        v = int(math.log(N//findeg, 2))-1
        m = N//PE

        BRAM = []
        BRTW = []

        for i in range(PE):
            BRAM.append([])
            for j in range(v):
                BRAM[i].append([])
                for k in range(m):
                    BRAM[i][j].append([])

        for i in range(PE//2):
            BRTW.append([])
            for j in range(v):
                BRTW[i].append([])
                for k in range(m):
                    BRTW[i][j].append([])

        bram_counter = 0
        ntru_counter = 0
    # ---------------------------------

    k = 0
    lena = findeg

    if ring == 0:
        # NWC
        v = int(log(N//findeg,2))
        lena_limit = (N//2)
    else:
        # NTRU
        powers_new = [_ for _ in powers]
        i = findeg
        r = 1
        while(i >= 6):
            powers_new = powers_new[N//(6*r):]
            i = (i//2)
            r = 2*r
        lena_limit = (N//4)

    while lena <= lena_limit:
        start = 0
        while start < N:
            if ring == 0:
                # NWC
                """
                W_pow = intReverse(k,v)+1
                TW = (-pow(w,W_pow,q)) % q # here, "-" means an extra w^(n/2)
                """
                W_pow = intReverse(k,v)+1
            else:
                # NTRU
                W_pow = powers_new[k]

            TW = pow(w,W_pow,q)
            k = k+1
            j = start
            while(j < (start + lena)):
                t = B[j]

                B[j       ] = (t + B[j + lena]) % q
                B[j + lena] = (t - B[j + lena]) % q
                B[j + lena] = B[j + lena]*TW % q

                if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1):  print("W: "+str(W_pow).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+lena).ljust(5))
                if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2):  print("W: "+str(W_pow).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+lena).ljust(5))
                if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3):  print("W: "+str(W_pow).ljust(5)+" A0: "+str(j).ljust(5)+" A1: "+str(j+lena).ljust(5))

                # ---------------------------------
                if ring == 0:
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+0][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = j
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+1][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = j+lena

                    BRTW[bram_counter & ((PE//2)-1)][bram_counter // (N//2)][(bram_counter & ((N//2)-1)) // (PE//2)] = W_pow

                    bram_counter = bram_counter + 1
                else:
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+0][bram_counter // (N//6)][3*((bram_counter & ((N//6)-1)) // (PE//2))+ntru_counter] = j
                    BRAM[(2*(bram_counter >> 0) & (PE-1))+1][bram_counter // (N//6)][3*((bram_counter & ((N//6)-1)) // (PE//2))+ntru_counter] = j+lena

                    BRTW[bram_counter & ((PE//2)-1)][bram_counter // (N//6)][3*((bram_counter & ((N//6)-1)) // (PE//2))+ntru_counter] = W_pow

                    bram_counter = (bram_counter + 1) if (ntru_counter == 2) else bram_counter
                    ntru_counter = 0 if (ntru_counter == 2) else (ntru_counter + 1)
                # ---------------------------------

                j = j+1
            start = j + lena
        lena = 2*lena

    N_inv = modinv(N//findeg,q)
    for i in range(N):
        B[i] = (B[i] * N_inv) % q

    return B,BRAM,BRTW

# A unified polynomial multiplication algorithm for all methods mentioned above
# ring: 0 --> NWC  (x^n + 1)
# -- findeg: 1
# -- findeg: 2
# -- findeg: 4
# --   ...
# ring: 1 --> NTRU (x^n - x^(n/2) + 1)
# -- findeg: 3
# -- findeg: 6
# -- findeg: 12
# --   ...
# NOTE: Later I can add PWC (x^n - 1)
# NOTE: Later I can add pure NTT/INTT operations
def CRTBasedModPolMul_Unified(A,B,w,w_inv,q,ring,findeg,ntrupowersf=[],ntrupowersb=[],ntrupowersi=[]):
    # --------------------------------------------- Initial step
    if ring == 0:
        """
        NWC requires no initial reduction operation
        """
        A_r = [_ for _ in A]
        B_r = [_ for _ in B]
    else:
        """
        NTRU requires initial reduction
        """
        A_r = [_ for _ in A]
        B_r = [_ for _ in B]

        if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3): print("---- Pre-processing:")

        wk = pow(w,len(A)//6,q)
        for i in range(len(A)//2):
            t1 = (wk*A[i+len(A)//2]) % q
            t2 = (wk*B[i+len(B)//2]) % q
            A_r[i+len(A)//2] = (A[i]+A[i+len(A)//2]-t1)%q
            A_r[i]           = (A[i]               +t1)%q
            B_r[i+len(B)//2] = (B[i]+B[i+len(B)//2]-t2)%q
            B_r[i]           = (B[i]               +t2)%q

            if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3):  print("W: "+str(len(A)//6).ljust(5)+" A0: "+str(i).ljust(5)+" A1: "+str(i+len(A)//2).ljust(5))

    # --------------------------------------------- NTT
    if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1):  print("---- NTT(A)")
    if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2):  print("---- NTT(A)")
    if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3):  print("---- NTT(A)")
    A_ntt,ABR,ATW = CRT_Iterative_Unified_NR(A_r,w,q,ring,findeg,ntrupowersf)
    if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1):  print("---- NTT(B)")
    if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2):  print("---- NTT(B)")
    if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3):  print("---- NTT(B)")
    B_ntt,BBR,BTW = CRT_Iterative_Unified_NR(B_r,w,q,ring,findeg,ntrupowersf)

    if DEBUG_TEST_NWC1:
        for an,bn in zip(A_ntt,B_ntt):
            DIN0_MFNTT_TXT.write(hex(an).replace("L","")[2:]+"\n")
            DIN1_MFNTT_TXT.write(hex(bn).replace("L","")[2:]+"\n")

    # --------------------------------------------- Degree-findeg modular polynomial multiplications
    C_ntt = [0 for _ in range(len(A))]
    if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1):  print("---- Coefficient-wise multiplication:")
    if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2):  print("---- Coefficient-wise multiplication:")
    if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3):  print("---- Coefficient-wise multiplication:")

    for i in range(len(A)//findeg):
        if ring == 0:
            # NWC
            w_pow = 2*intReverse(i,int(log(len(A)//findeg,2)))+1
        else:
            # NTRU
            w_pow = ((findeg//3)*ntrupowersb[i*(findeg//3)])
            # w_pow = ((findeg//3)*ntrupowersb[i*(findeg//3)])//(findeg//3)

        if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1):  print("A: {}".format(i))
        if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2):  print("W: "+str(w_pow).ljust(5)+" A: {}".format(range(findeg*i,findeg*i+findeg)))
        if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3):  print("W: "+str(w_pow).ljust(5)+" A: {}".format(range(findeg*i,findeg*i+findeg)))

        wk    = pow(w,w_pow,q)
        C_ntt[findeg*i:findeg*i+findeg] = PolWiseMult(A_ntt[findeg*i:findeg*i+findeg],B_ntt[findeg*i:findeg*i+findeg],wk,findeg,q)

    # --------------------------------------------- INTT
    if DEBUG_TEST_NWC1:
        for cn in C_ntt:
            DOUT_MINTT_TXT.write(hex(cn).replace("L","")[2:]+"\n")

    if DEBUG_DISP_NWC1  and (ring == 0) and (findeg == 1): print("---- INTT(C)")
    if DEBUG_DISP_NWC2  and (ring == 0) and (findeg == 2): print("---- INTT(C)")
    if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3): print("---- INTT(C)")
    C,CBR,CTW = ICRT_Iterative_Unified_RN(C_ntt,w_inv,q,ring,findeg,ntrupowersi)

    # --------------------------------------------- Final step
    if ring == 0:
        """
        NWC requires no final reconstruction step
        """
        return C,ABR,ATW,CBR,CTW
    else:
        if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3): print("---- Post-processing:")
        """
        NTRU requires final reconstruction step
        """
        wk = modinv((2*pow(w,(len(C)//6),q)-1)%q,q)
        # wk = modinv((2*pow(w,(len(C)//6)//(findeg//3),q)-1)%q,q)

        for i in range(len(C)//2):
            t = ((C[i]-C[i+len(C)//2])*wk)%q   # t = f[i+N//2]
            C[i          ] = (C[i]+C[i+len(C)//2]-t)%q
            C[i+len(C)//2] = (2*t)%q

            if DEBUG_DISP_NTRU3 and (ring == 1) and (findeg == 3): print("W: "+str("P").ljust(5)+" A0: "+str(i).ljust(5)+" A1: "+str(i+len(C)//2).ljust(5))

        return C,ABR,ATW,CBR,CTW

# ----------------------------------------------------------

# Demo
# Random A,B
A = [randint(0,q-1) for _ in range(n)]
B = [randint(0,q-1) for _ in range(n)]

# Random A,B (for ntru)
A_ntru = [randint(0,mq-1) for _ in range(m)]
B_ntru = [randint(0,mq-1) for _ in range(m)]

# reduce functions
pwc  = [-1]+[0]*(n-1)+[1]
nwc  =  [1]+[0]*(n-1)+[1]
ntru =  [1]+[0]*(int(m/2)-1)+[-1]+[0]*(int(m/2)-1)+[1]

# NTT
if DEBUG_DISP_FNTT: print("\n-------- Addressing for NTT --------")
if DEBUG_FNTT     : N0,N0BR,N0TW = Radix2_DIT_Iterative_NTT_NR(A,w,q)
if DEBUG_DISP_INTT: print("\n-------- Addressing for INTT --------")
if DEBUG_INTT     : N1,N1BR,N1TW = Radix2_DIF_Iterative_INTT_RN(N0,w_inv,q)

# POLMUL
if DEBUG_DISP_NWC1: print("\n-------- Addressing for NWC - findeg=1 --------")
if DEBUG_NWC1     : R0,R0BRF,R0TWF,R0BRI,R0TWI = CRTBasedModPolMul_Unified(A,B,psi,psi_inv,q,ring=0,findeg=1) # NWC - findeg=1
if DEBUG_DISP_NWC2: print("\n-------- Addressing for NWC - findeg=2 --------")
if DEBUG_NWC2     : R1,R1BRF,R1TWF,R1BRI,R1TWI = CRTBasedModPolMul_Unified(A,B,w,w_inv,q,ring=0,findeg=2) # NWC - findeg=2
if DEBUG_DISP_NTRU3: print("\n-------- Addressing for NTRU - findeg=3 --------")
ring, findeg = 1,3
if DEBUG_NTRU3    : R2,R2BRF,R2TWF,R2BRI,R2TWI = CRTBasedModPolMul_Unified(A_ntru,B_ntru,mw,mw_inv,mq,ring,findeg,ntrupowersf,ntrupowersb,ntrupowersi) # NTRU - findeg=3

# Print memory structure
def PrintBRAM(BRAM,ring=0,findeg=1):
    if ring == 0:
        v = int(math.log(n//findeg, 2))
        m = n//PE
    else:
        v = int(math.log((3*n)//findeg, 2))-1
        m = (3*n)//PE
    BS = ""
    for j in range(v):
        BS = BS+"*************************************************** stage="+str(j)+"\n"
        BS = BS+"BRAM:"

        for i in range(PE//2):
            BS = BS+"\t|"+str(2*i).ljust(5)+str(2*i+1).ljust(4)+"|"
        BS = BS+"\n"
        BS = BS+"     "
        for i in range(PE//2):
            BS = BS+"\t----------"
        BS = BS+"\n"

        for k in range(m):
            BS = BS + "AD"+str(k)+":"
            for i in range(PE//2):
                BS = BS+"\t|"+str(BRAM[2*i][j][k]).ljust(5)+str(BRAM[2*i+1][j][k]).ljust(4)+"|"
            BS = BS+"\n"

    return BS

def PrintBRTW(BRTW,ring=0,findeg=1):
    if ring == 0:
        v = int(math.log(n//findeg, 2))
        m = n//PE
    else:
        v = int(math.log((3*n)//findeg, 2))-1
        m = (3*n)//PE
    TS = ""
    for j in range(v):
        TS = TS+"*************************************************** stage="+str(j)+"\n"
        TS = TS+"TWID:"

        for i in range(PE//2):
            TS = TS+"\t|"+str(i).ljust(5)+"|"
        TS = TS+"\n"
        TS = TS+"     "
        for i in range(PE//2):
            TS = TS+"\t------"
        TS = TS+"\n"

        for k in range(m):
            TS = TS + "AD"+str(k)+":"
            for i in range(PE//2):
                TS = TS+"\t|"+str(BRTW[i][j][k]).ljust(5)+"|"
            TS = TS+"\n"

    return TS

# Write to txt
if DEBUG_MEMP_FNTT or DEBUG_MEMP_INTT or DEBUG_MEMP_NWC1 or DEBUG_MEMP_NWC2 or DEBUG_MEMP_NTRU3:
    print("")
    print("-------- Generated:")

if DEBUG_MEMP_FNTT:
    FNTT_BR   = PrintBRAM(N0BR)
    FNTT_TW   = PrintBRTW(N0TW)

    # Data
    FNTT_BR_TXT = open("FNTT_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    FNTT_BR_TXT.write(FNTT_BR)
    FNTT_BR_TXT.close()
    # Twiddle
    FNTT_TW_TXT = open("FNTT_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    FNTT_TW_TXT.write(FNTT_TW)
    FNTT_TW_TXT.close()

    print("* FNTT_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")
    print("* FNTT_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")

if DEBUG_MEMP_INTT:
    INTT_BR   = PrintBRAM(N1BR)
    INTT_TW   = PrintBRTW(N1TW)

    # Date
    INTT_BR_TXT = open("INTT_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    INTT_BR_TXT.write(INTT_BR)
    INTT_BR_TXT.close()
    # Twiddle
    INTT_TW_TXT = open("INTT_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    INTT_TW_TXT.write(INTT_TW)
    INTT_TW_TXT.close()

    print("* INTT_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")
    print("* INTT_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")

if DEBUG_MEMP_NWC1:
    NWC1F_BR  = PrintBRAM(R0BRF)
    NWC1I_BR  = PrintBRAM(R0BRI)
    NWC1F_TW  = PrintBRTW(R0TWF)
    NWC1I_TW  = PrintBRTW(R0TWI)

    # Data
    NWC1_BR_TXT = open("NWC1_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    NWC1_BR_TXT.write("---------------------------------------------------------------------- Forward NTT (x2)\n")
    NWC1_BR_TXT.write(NWC1F_BR)
    NWC1_BR_TXT.write("---------------------------------------------------------------------- Coefficient-wise multiplication\n")
    NWC1_BR_TXT.write("---------------------------------------------------------------------- Inverse NTT\n")
    NWC1_BR_TXT.write(NWC1I_BR)
    NWC1_BR_TXT.close()
    # Twiddle
    NWC1_TW_TXT = open("NWC1_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    NWC1_TW_TXT.write("---------------------------------------------------------------------- Forward NTT (x2)\n")
    NWC1_TW_TXT.write(NWC1F_TW)
    NWC1_TW_TXT.write("---------------------------------------------------------------------- Coefficient-wise multiplication\n")
    NWC1_TW_TXT.write("---------------------------------------------------------------------- Inverse NTT\n")
    NWC1_TW_TXT.write(NWC1I_TW)
    NWC1_TW_TXT.close()

    print("* NWC1_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")
    print("* NWC1_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")

if DEBUG_MEMP_NWC2:
    NWC2F_BR  = PrintBRAM(R1BRF,0,2)
    NWC2I_BR  = PrintBRAM(R1BRI,0,2)
    NWC2F_TW  = PrintBRTW(R1TWF,0,2)
    NWC2I_TW  = PrintBRTW(R1TWI,0,2)

    # Data
    NWC2_BR_TXT = open("NWC2_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    NWC2_BR_TXT.write("---------------------------------------------------------------------- Forward NTT (x2)\n")
    NWC2_BR_TXT.write(NWC2F_BR)
    NWC2_BR_TXT.write("---------------------------------------------------------------------- Degree-2 polynomial-wise multiplication\n")
    NWC2_BR_TXT.write("---------------------------------------------------------------------- Inverse NTT\n")
    NWC2_BR_TXT.write(NWC2I_BR)
    NWC2_BR_TXT.close()
    # Twiddle
    NWC2_TW_TXT = open("NWC2_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt","w")
    NWC2_TW_TXT.write("---------------------------------------------------------------------- Forward NTT (x2)\n")
    NWC2_TW_TXT.write(NWC2F_TW)
    NWC2_TW_TXT.write("---------------------------------------------------------------------- Degree-2 polynomial-wise multiplication\n")
    NWC2_TW_TXT.write("---------------------------------------------------------------------- Inverse NTT\n")
    NWC2_TW_TXT.write(NWC2I_TW)
    NWC2_TW_TXT.close()

    print("* NWC2_mem_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")
    print("* NWC2_tw_N"+str(n)+"_PE"+str(PE_NUMBER)+".txt")

if DEBUG_MEMP_NTRU3:
    NTRU3F_BR = PrintBRAM(R2BRF,1,3)
    NTRU3I_BR = PrintBRAM(R2BRI,1,3)
    NTRU3F_TW = PrintBRTW(R2TWF,1,3)
    NTRU3I_TW = PrintBRTW(R2TWI,1,3)

    # Data
    NTRU3_BR_TXT = open("NTRU3_mem_N"+str(m)+"_PE"+str(PE_NUMBER)+".txt","w")
    NTRU3_BR_TXT.write("---------------------------------------------------------------------- First Reduction\n")
    NTRU3_BR_TXT.write("---------------------------------------------------------------------- Forward NTT (x2)\n")
    NTRU3_BR_TXT.write(NTRU3F_BR)
    NTRU3_BR_TXT.write("---------------------------------------------------------------------- Degree-2 polynomial-wise multiplication\n")
    NTRU3_BR_TXT.write("---------------------------------------------------------------------- Inverse NTT\n")
    NTRU3_BR_TXT.write(NTRU3I_BR)
    NTRU3_BR_TXT.write("---------------------------------------------------------------------- Final Reconstruction\n")
    NTRU3_BR_TXT.close()
    # Twiddle
    NTRU3_TW_TXT = open("NTRU3_tw_N"+str(m)+"_PE"+str(PE_NUMBER)+".txt","w")
    NTRU3_TW_TXT.write("---------------------------------------------------------------------- First Reduction\n")
    NTRU3_TW_TXT.write("---------------------------------------------------------------------- Forward NTT (x2)\n")
    NTRU3_TW_TXT.write(NTRU3F_TW)
    NTRU3_TW_TXT.write("---------------------------------------------------------------------- Degree-2 polynomial-wise multiplication\n")
    NTRU3_TW_TXT.write("---------------------------------------------------------------------- Inverse NTT\n")
    NTRU3_TW_TXT.write(NTRU3I_TW)
    NTRU3_TW_TXT.write("---------------------------------------------------------------------- Final Reconstruction\n")
    NTRU3_TW_TXT.close()

    print("* NTRU3_mem_N"+str(m)+"_PE"+str(PE_NUMBER)+".txt")
    print("* NTRU3_tw_N"+str(m)+"_PE"+str(PE_NUMBER)+".txt")

# Generate test vectors
if DEBUG_TEST_FNTT:
    # Parameters
    PRM_TXT.write(hex(n        ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(q        ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(q_bit    ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(PE_NUMBER).replace("L","")[2:]+"\n")
    # Input/Output
    for fnn,inn in zip(A,N0):
        DIN_TXT.write(hex(fnn).replace("L","")[2:]+"\n")
        DOUT_TXT.write(hex(inn).replace("L","")[2:]+"\n")
    # Twiddle factors
    v = int(math.log(n, 2))
    m = n//PE
    for j in range(v):
        for k in range(0,m,max(1,m>>j)):
            for i in range(PE//2):
                if WLMONT:
                    W_TXT.write(hex(pow(w,N0TW[i][j][k],q)*fd1_R % q).replace("L","")[2:]+"\t  ")
                else:
                    W_TXT.write(hex(pow(w,N0TW[i][j][k],q)).replace("L","")[2:]+"\t  ")
            W_TXT.write("\n")
if DEBUG_TEST_INTT:
    # Parameters
    PRM_TXT.write(hex(n          ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(q          ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(q_bit      ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(PE_NUMBER  ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(modinv(n,q)).replace("L","")[2:]+"\n")
    # Input/Output
    for fnn,inn in zip(N0,A):
        DIN_TXT.write(hex(fnn).replace("L","")[2:]+"\n")
        DOUT_TXT.write(hex(inn).replace("L","")[2:]+"\n")
    # Twiddle factors
    v = int(math.log(n, 2))
    m = n//PE
    for j in range(v):
        for k in range(0,m,min(m,1<<j)):
            for i in range(PE//2):
                if WLMONT:
                    WINV_TXT.write(hex(pow(w,N1TW[i][j][k],q)*fd1_R % q).replace("L","")[2:]+"\t  ")
                else:
                    WINV_TXT.write(hex(pow(w,N1TW[i][j][k],q)).replace("L","")[2:]+"\t  ")
            WINV_TXT.write("\n")
if DEBUG_TEST_NWC1:
    # Parameters
    PRM_TXT.write(hex(n          ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(q          ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(q_bit      ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(PE_NUMBER  ).replace("L","")[2:]+"\n")
    PRM_TXT.write(hex(modinv(n,q)).replace("L","")[2:]+"\n")
    # Input/Output
    for mi0,mi1,mo0 in zip(A,B,R0):
        DIN0_TXT.write(hex(mi0).replace("L","")[2:]+"\n")
        DIN1_TXT.write(hex(mi1).replace("L","")[2:]+"\n")
        DOUT_TXT.write(hex(mo0).replace("L","")[2:]+"\n")
    # Twiddle factor
    v = int(math.log(n, 2))
    m = n//PE
    for j in range(v):
        for k in range(0,m,max(1,m>>j)):
            for i in range(PE//2):
                if WLMONT:
                    W_TXT.write(hex(pow(w,R0TWF[i][j][k],q)*fd1_R % q).replace("L","")[2:]+"\t  ")
                else:
                    W_TXT.write(hex(pow(w,R0TWF[i][j][k],q)).replace("L","")[2:]+"\t  ")
            W_TXT.write("\n")
    v = int(math.log(n, 2))
    m = n//PE
    for j in range(v):
        for k in range(0,m,max(1,m>>(v-j-1))):
            for i in range(PE//2):
                if WLMONT:
                    WINV_TXT.write(hex(pow(w,R0TWI[i][j][k],q)*fd1_R % q).replace("L","")[2:]+"\t  ")
                else:
                    WINV_TXT.write(hex(pow(w,R0TWI[i][j][k],q)).replace("L","")[2:]+"\t  ")
            WINV_TXT.write("\n")
    # R1TWF
if DEBUG_TEST_NWC2:
    pass
if DEBUG_TEST_NTRU3:
    pass

#
