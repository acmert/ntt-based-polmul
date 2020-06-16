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
REF = Evaluator.NaiveNTT_NN(A,w,q)

# Reversed N0
REF_rev = indexReverse(REF,int(log(n,2)))

# NTT operations
N0 = Evaluator.NaiveNTT_NR(A,w,q)
N1 = Evaluator.Radix2_DIT_Recursive_NTT(A,w,q)
N2 = Evaluator.Radix2_DIF_Recursive_NTT(A,w,q)
N3 = Evaluator.Radix2_DIF_Iterative_NTT_NR(A,w,q)
N4 = Evaluator.Radix2_DIF_Iterative_NTT_RN(A_rev,w,q)
N5 = Evaluator.Radix2_DIF_Iterative_NTT_NN(A,w,q)
N6 = Evaluator.Radix2_DIT_Iterative_NTT_NR(A,w,q)
N7 = Evaluator.Radix2_DIT_Iterative_NTT_RN(A_rev,w,q)
N8 = Evaluator.Radix2_DIT_Iterative_NTT_NN(A,w,q)
N9 = Evaluator.CRT_Recursive_NTT(A,w,q)
N10= Evaluator.CRT_Full_NTT(A,w,q)
N11= Evaluator.CooleyTukeyNTT(A,w,q)
N12= Evaluator.IterativeNTT(A,w,q)
N13= Evaluator.FourStepNTT(A,w,q,size)
N14= Evaluator.FourStepNTTv2(A,w_inv,q,size)
N15= Evaluator.CTBased_ConstantGeometryNTT(A_rev,w,q)

# Check NTT
print("-------- Sanity check for NTT operations --------")
# print("A         : {}".format(A))
# print("br(A)     : {}".format(A_rev))
# print("NTT(A)    : {}".format(REF))
# print("br(NTT(A)): {}".format(REF_rev))
print("")
print("NaiveNTT_NR                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N0)]) == 0) else "Wrong"))
print("Radix2_DIT_Recursive_NTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N1)]) == 0) else "Wrong"))
print("Radix2_DIF_Recursive_NTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N2)]) == 0) else "Wrong"))
print("Radix2_DIF_Iterative_NTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N3)]) == 0) else "Wrong"))
print("Radix2_DIF_Iterative_NTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N4)]) == 0) else "Wrong"))
print("Radix2_DIF_Iterative_NTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N5)]) == 0) else "Wrong"))
print("Radix2_DIT_Iterative_NTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N6)]) == 0) else "Wrong"))
print("Radix2_DIT_Iterative_NTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N7)]) == 0) else "Wrong"))
print("Radix2_DIT_Iterative_NTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N8)]) == 0) else "Wrong"))
print("CRT_Recursive_NTT              -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N9)]) == 0) else "Wrong"))
print("CRT_Full_NTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N10)]) == 0) else "Wrong"))
print("CooleyTukeyNTT                 -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N11)]) == 0) else "Wrong"))
print("IterativeNTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF_rev,N12)]) == 0) else "Wrong"))
print("FourStepNTT                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N13)]) == 0) else "Wrong"))
print("FourStepNTTv2                  -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N14)]) == 0) else "Wrong"))
print("CTBased_ConstantGeometryNTT    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(REF,N15)]) == 0) else "Wrong"))
print("")

# INTT operations
R0 = Evaluator.NaiveINTT_NR(REF,w_inv,q)
R1 = Evaluator.Radix2_DIT_Recursive_INTT(REF,w_inv,q)
R2 = Evaluator.Radix2_DIF_Recursive_INTT(REF,w_inv,q)
R3 = Evaluator.Radix2_DIF_Iterative_INTT_NR(REF,w_inv,q)
R4 = Evaluator.Radix2_DIF_Iterative_INTT_RN(REF_rev,w_inv,q)
R5 = Evaluator.Radix2_DIF_Iterative_INTT_NN(REF,w_inv,q)
R6 = Evaluator.Radix2_DIT_Iterative_INTT_NR(REF,w_inv,q)
R7 = Evaluator.Radix2_DIT_Iterative_INTT_RN(REF_rev,w_inv,q)
R8 = Evaluator.Radix2_DIT_Iterative_INTT_NN(REF,w_inv,q)
R9 = Evaluator.CRT_Recursive_INTT(REF,w_inv,q)
R10= Evaluator.CRT_Full_INTT(REF,w_inv,q)
R11= Evaluator.CooleyTukeyINTT(REF,w_inv,q)
R12= Evaluator.IterativeINTT(REF,w_inv,q)
R13= Evaluator.FourStepINTT(REF,w_inv,q,size)
R14= Evaluator.FourStepINTTv2(REF,w,q,size)
R15= Evaluator.CTBased_ConstantGeometryINTT(REF_rev,w_inv,q)

# Check INTT
print("-------- Sanity check for INTT operations --------")
# print("NTT(A)    : {}".format(REF))
# print("br(NTT(A)): {}".format(REF_rev))
# print("A         : {}".format(A))
# print("br(A)     : {}".format(A_rev))
print("")
print("NaiveINTT_NR                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R0)]) == 0) else "Wrong"))
print("Radix2_DIT_Recursive_INTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R1)]) == 0) else "Wrong"))
print("Radix2_DIF_Recursive_INTT       -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R2)]) == 0) else "Wrong"))
print("Radix2_DIF_Iterative_INTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R3)]) == 0) else "Wrong"))
print("Radix2_DIF_Iterative_INTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R4)]) == 0) else "Wrong"))
print("Radix2_DIF_Iterative_INTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R5)]) == 0) else "Wrong"))
print("Radix2_DIT_Iterative_INTT_NR    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R6)]) == 0) else "Wrong"))
print("Radix2_DIT_Iterative_INTT_RN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R7)]) == 0) else "Wrong"))
print("Radix2_DIT_Iterative_INTT_NN    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R8)]) == 0) else "Wrong"))
print("CRT_Recursive_INTT              -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R9)]) == 0) else "Wrong"))
print("CRT_Full_INTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R10)]) == 0) else "Wrong"))
print("CooleyTukeyINTT                 -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R11)]) == 0) else "Wrong"))
print("IterativeINTT                   -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A_rev,R12)]) == 0) else "Wrong"))
print("FourStepINTT                    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R13)]) == 0) else "Wrong"))
print("FourStepINTTv2                  -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R14)]) == 0) else "Wrong"))
print("CTBased_ConstantGeometryINTT    -->" + ("Correct" if(sum([abs(x-y) for x,y in zip(A,R15)]) == 0) else "Wrong"))
print("")

#
