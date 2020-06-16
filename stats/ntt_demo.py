from math import log
from random import randint

from generate_prime import *
from helper import *
from ntt import *

# Parameter generation

# Determine n and bit-size of q, then find a q satisfying
# the condition: q = 1 (mod 2n) or q = 1 (mod n)
#
# Based on n and q, generate NTT parameters

mod     = 2 # if 1 --> q = 1 (mod n), if 2 --> q = 1 (mod 2n)
n       = 64
size    = [8,8]
q_bit   = 10

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
print("Parameters")
print("n      : {}".format(n))
print("q      : {}".format(q))
print("w      : {}".format(w))
print("w_inv  : {}".format(w_inv))
print("psi    : {}".format(psi))
print("psi_inv: {}".format(psi_inv))
print("")

#NOTE: Comment Out Naive Method for Large Parameters

# Demo
# Random A
A = [randint(0,q-1) for x in range(n)]

# Reversed A
A_rev = indexReverse(A,int(log(n,2)))

# NTT operation
Evaluator = NTT()

# Reference NTT operation
REF,REFMul,REFAdd,REFSub,REFBtf = Evaluator.NaiveNTT_NN(A,w,q)

# Reversed N0
REF_rev = indexReverse(REF,int(log(n,2)))

# NTT operations
N0 ,N0Mul ,N0Add ,N0Sub ,N0Btf  = Evaluator.NaiveNTT_NR(A,w,q)
N1 ,N1Mul ,N1Add ,N1Sub ,N1Btf  = Evaluator.Radix2_DIT_Recursive_NTT(A,w,q)
N2 ,N2Mul ,N2Add ,N2Sub ,N2Btf  = Evaluator.Radix2_DIF_Recursive_NTT(A,w,q)
N3 ,N3Mul ,N3Add ,N3Sub ,N3Btf  = Evaluator.Radix2_DIF_Iterative_NTT_NR(A,w,q)
N4 ,N4Mul ,N4Add ,N4Sub ,N4Btf  = Evaluator.Radix2_DIF_Iterative_NTT_RN(A_rev,w,q)
N5 ,N5Mul ,N5Add ,N5Sub ,N5Btf  = Evaluator.Radix2_DIF_Iterative_NTT_NN(A,w,q)
N6 ,N6Mul ,N6Add ,N6Sub ,N6Btf  = Evaluator.Radix2_DIT_Iterative_NTT_NR(A,w,q)
N7 ,N7Mul ,N7Add ,N7Sub ,N7Btf  = Evaluator.Radix2_DIT_Iterative_NTT_RN(A_rev,w,q)
N8 ,N8Mul ,N8Add ,N8Sub ,N8Btf  = Evaluator.Radix2_DIT_Iterative_NTT_NN(A,w,q)
N9 ,N9Mul ,N9Add ,N9Sub ,N9Btf  = Evaluator.CRT_Recursive_NTT(A,w,q)
N10,N10Mul,N10Add,N10Sub,N10Btf = Evaluator.CRT_Full_NTT(A,w,q)
N11,N11Mul,N11Add,N11Sub,N11Btf = Evaluator.CooleyTukeyNTT(A,w,q)
N12,N12Mul,N12Add,N12Sub,N12Btf = Evaluator.IterativeNTT(A,w,q)
N13,N13Mul,N13Add,N13Sub,N13Btf = Evaluator.FourStepNTT(A,w,q,size)
N14,N14Mul,N14Add,N14Sub,N14Btf = Evaluator.FourStepNTTv2(A,w_inv,q,size)
N15,N15Mul,N15Add,N15Sub,N15Btf = Evaluator.CTBased_ConstantGeometryNTT(A_rev,w,q)

# Check NTT
print("-------- Sanity check for NTT operations --------")
# print("A         : {}".format(A))
# print("br(A)     : {}".format(A_rev))
# print("NTT(A)    : {}".format(REF))
# print("br(NTT(A)): {}".format(REF_rev))
print("")
print("NaiveNTT_NR                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N0)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N0Mul))
# print("-- Add:{}".format(N0Add))
# print("-- Sub:{}".format(N0Sub))
# print("-- Btf:{}".format(N0Btf))
print("Radix2_DIT_Recursive_NTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N1Mul))
# print("-- Add:{}".format(N1Add))
# print("-- Sub:{}".format(N1Sub))
# print("-- Btf:{}".format(N1Btf))
print("Radix2_DIF_Recursive_NTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N2Mul))
# print("-- Add:{}".format(N2Add))
# print("-- Sub:{}".format(N2Sub))
# print("-- Btf:{}".format(N2Btf))
print("Radix2_DIF_Iterative_NTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N3)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N3Mul))
# print("-- Add:{}".format(N3Add))
# print("-- Sub:{}".format(N3Sub))
# print("-- Btf:{}".format(N3Btf))
print("Radix2_DIF_Iterative_NTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N4)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N4Mul))
# print("-- Add:{}".format(N4Add))
# print("-- Sub:{}".format(N4Sub))
# print("-- Btf:{}".format(N4Btf))
print("Radix2_DIF_Iterative_NTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N5)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N5Mul))
# print("-- Add:{}".format(N5Add))
# print("-- Sub:{}".format(N5Sub))
# print("-- Btf:{}".format(N5Btf))
print("Radix2_DIT_Iterative_NTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N6)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N6Mul))
# print("-- Add:{}".format(N6Add))
# print("-- Sub:{}".format(N6Sub))
# print("-- Btf:{}".format(N6Btf))
print("Radix2_DIT_Iterative_NTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N7)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N7Mul))
# print("-- Add:{}".format(N7Add))
# print("-- Sub:{}".format(N7Sub))
# print("-- Btf:{}".format(N7Btf))
print("Radix2_DIT_Iterative_NTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N8)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N8Mul))
# print("-- Add:{}".format(N8Add))
# print("-- Sub:{}".format(N8Sub))
# print("-- Btf:{}".format(N8Btf))
print("CRT_Recursive_NTT              -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N9)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N9Mul))
# print("-- Add:{}".format(N9Add))
# print("-- Sub:{}".format(N9Sub))
# print("-- Btf:{}".format(N9Btf))
print("CRT_Full_NTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N10)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N10Mul))
# print("-- Add:{}".format(N10Add))
# print("-- Sub:{}".format(N10Sub))
# print("-- Btf:{}".format(N10Btf))
print("CooleyTukeyNTT                 -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N11)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N11Mul))
# print("-- Add:{}".format(N11Add))
# print("-- Sub:{}".format(N11Sub))
# print("-- Btf:{}".format(N11Btf))
print("IterativeNTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N12)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N12Mul))
# print("-- Add:{}".format(N12Add))
# print("-- Sub:{}".format(N12Sub))
# print("-- Btf:{}".format(N12Btf))
print("FourStepNTT                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N13)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N13Mul))
# print("-- Add:{}".format(N13Add))
# print("-- Sub:{}".format(N13Sub))
# print("-- Btf:{}".format(N13Btf))
print("FourStepNTTv2                  -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N14)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N14Mul))
# print("-- Add:{}".format(N14Add))
# print("-- Sub:{}".format(N14Sub))
# print("-- Btf:{}".format(N14Btf))
print("CTBased_ConstantGeometryNTT    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N15)]) == 0) else "Wrong"))
print("-- Mul:{}".format(N15Mul))
# print("-- Add:{}".format(N15Add))
# print("-- Sub:{}".format(N15Sub))
# print("-- Btf:{}".format(N15Btf))
print("")

# INTT operations
R0 ,R0Mul ,R0Add ,R0Sub ,R0Btf  = Evaluator.NaiveINTT_NR(REF,w_inv,q)
R1 ,R1Mul ,R1Add ,R1Sub ,R1Btf  = Evaluator.Radix2_DIT_Recursive_INTT(REF,w_inv,q)
R2 ,R2Mul ,R2Add ,R2Sub ,R2Btf  = Evaluator.Radix2_DIF_Recursive_INTT(REF,w_inv,q)
R3 ,R3Mul ,R3Add ,R3Sub ,R3Btf  = Evaluator.Radix2_DIF_Iterative_INTT_NR(REF,w_inv,q)
R4 ,R4Mul ,R4Add ,R4Sub ,R4Btf  = Evaluator.Radix2_DIF_Iterative_INTT_RN(REF_rev,w_inv,q)
R5 ,R5Mul ,R5Add ,R5Sub ,R5Btf  = Evaluator.Radix2_DIF_Iterative_INTT_NN(REF,w_inv,q)
R6 ,R6Mul ,R6Add ,R6Sub ,R6Btf  = Evaluator.Radix2_DIT_Iterative_INTT_NR(REF,w_inv,q)
R7 ,R7Mul ,R7Add ,R7Sub ,R7Btf  = Evaluator.Radix2_DIT_Iterative_INTT_RN(REF_rev,w_inv,q)
R8 ,R8Mul ,R8Add ,R8Sub ,R8Btf  = Evaluator.Radix2_DIT_Iterative_INTT_NN(REF,w_inv,q)
R9 ,R9Mul ,R9Add ,R9Sub ,R9Btf  = Evaluator.CRT_Recursive_INTT(REF,w_inv,q)
R10,R10Mul,R10Add,R10Sub,R10Btf = Evaluator.CRT_Full_INTT(REF,w_inv,q)
R11,R11Mul,R11Add,R11Sub,R11Btf = Evaluator.CooleyTukeyINTT(REF,w_inv,q)
R12,R12Mul,R12Add,R12Sub,R12Btf = Evaluator.IterativeINTT(REF,w_inv,q)
R13,R13Mul,R13Add,R13Sub,R13Btf = Evaluator.FourStepINTT(REF,w_inv,q,size)
R14,R14Mul,R14Add,R14Sub,R14Btf = Evaluator.FourStepINTTv2(REF,w,q,size)
R15,R15Mul,R15Add,R15Sub,R15Btf = Evaluator.CTBased_ConstantGeometryINTT(REF_rev,w_inv,q)

# Check INTT
print("-------- Sanity check for INTT operations --------")
# print("NTT(A)    : {}".format(REF))
# print("br(NTT(A)): {}".format(REF_rev))
# print("A         : {}".format(A))
# print("br(A)     : {}".format(A_rev))
print("")
print("NaiveINTT_NR                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R0)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R0Mul))
# print("-- Add:{}".format(R0Add))
# print("-- Sub:{}".format(R0Sub))
# print("-- Btf:{}".format(R0Btf))
print("Radix2_DIT_Recursive_INTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R1)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R1Mul))
# print("-- Add:{}".format(R1Add))
# print("-- Sub:{}".format(R1Sub))
# print("-- Btf:{}".format(R1Btf))
print("Radix2_DIF_Recursive_INTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R2)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R2Mul))
# print("-- Add:{}".format(R2Add))
# print("-- Sub:{}".format(R2Sub))
# print("-- Btf:{}".format(R2Btf))
print("Radix2_DIF_Iterative_INTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R3)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R3Mul))
# print("-- Add:{}".format(R3Add))
# print("-- Sub:{}".format(R3Sub))
# print("-- Btf:{}".format(R3Btf))
print("Radix2_DIF_Iterative_INTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R4)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R4Mul))
# print("-- Add:{}".format(R4Add))
# print("-- Sub:{}".format(R4Sub))
# print("-- Btf:{}".format(R4Btf))
print("Radix2_DIF_Iterative_INTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R5)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R5Mul))
# print("-- Add:{}".format(R5Add))
# print("-- Sub:{}".format(R5Sub))
# print("-- Btf:{}".format(R5Btf))
print("Radix2_DIT_Iterative_INTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R6)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R6Mul))
# print("-- Add:{}".format(R6Add))
# print("-- Sub:{}".format(R6Sub))
# print("-- Btf:{}".format(R6Btf))
print("Radix2_DIT_Iterative_INTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R7)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R7Mul))
# print("-- Add:{}".format(R7Add))
# print("-- Sub:{}".format(R7Sub))
# print("-- Btf:{}".format(R7Btf))
print("Radix2_DIT_Iterative_INTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R8)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R8Mul))
# print("-- Add:{}".format(R8Add))
# print("-- Sub:{}".format(R8Sub))
# print("-- Btf:{}".format(R8Btf))
print("CRT_Recursive_INTT              -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R9)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R9Mul))
# print("-- Add:{}".format(R9Add))
# print("-- Sub:{}".format(R9Sub))
# print("-- Btf:{}".format(R9Btf))
print("CRT_Full_INTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R10)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R10Mul))
# print("-- Add:{}".format(R10Add))
# print("-- Sub:{}".format(R10Sub))
# print("-- Btf:{}".format(R10Btf))
print("CooleyTukeyINTT                 -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R11)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R11Mul))
# print("-- Add:{}".format(R11Add))
# print("-- Sub:{}".format(R11Sub))
# print("-- Btf:{}".format(R11Btf))
print("IterativeINTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R12)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R12Mul))
# print("-- Add:{}".format(R12Add))
# print("-- Sub:{}".format(R12Sub))
# print("-- Btf:{}".format(R12Btf))
print("FourStepINTT                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R13)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R13Mul))
# print("-- Add:{}".format(R13Add))
# print("-- Sub:{}".format(R13Sub))
# print("-- Btf:{}".format(R13Btf))
print("FourStepINTTv2                  -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R14)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R14Mul))
# print("-- Add:{}".format(R14Add))
# print("-- Sub:{}".format(R14Sub))
# print("-- Btf:{}".format(R14Btf))
print("CTBased_ConstantGeometryINTT    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R15)]) == 0) else "Wrong"))
print("-- Mul:{}".format(R15Mul))
# print("-- Add:{}".format(R15Add))
# print("-- Sub:{}".format(R15Sub))
# print("-- Btf:{}".format(R15Btf))
print("")

#
