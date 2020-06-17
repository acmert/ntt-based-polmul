from helper import *
from ntt import *

class Poly:
    """
    * These are baseline (not optimized) implementations *

    Reference Implementations
    --Reference Polynomial Multiplication (School-Book)
    --Reference Modular Polynomial Multiplication (School-Book)
    ---- Reduction polynomial: x^n - 1           (Positive wrapped convolution - PWC)
    ---- Reduction polynomial: x^n + 1           (Negative wrapped convolution - NWC)
    ---- Reduction polynomial: x^n - x^(n/2) + 1 (NTRU)
    --Reference Polynomial Reduction (with any reduction polynomial)

    NTT-based Polynomial Multiplication (no reduction)
    -- A(x),B(x): n-1 degree polynomials
    -- C(x)     : 2n-2 degree polynomial
    -- A(x),B(x) should be zero-padded to 2n before operation
    -- C(x)=A(x)*B(X) --> C=INTT_2n(NTT_2n(zero-padded A) . NTT_2n(zero-padded B))
    -- If there is a ring Z_q[x]/f(x), polynomial reduction should be applied separately

    NTT-based Modular Polynomial Multiplication with Carrier Modulus
    -- Polynomial reduction operation is performed separately
    -- If n is a power-of-two
    ---- Select a new ntt-friendly q with bit-size > log(n*q^2)
    ---- Perform "NTT-based Polynomial Multiplication" and then apply reduction
    -- If n is NOT a power-of-two
    ---- Zero-pad input polynomials to the closest power-of-two
    ---- Select a new ntt-friendly q with bit-size > log(n*q^2)
    ---- Perform "NTT-based Polynomial Multiplication" and then apply reduction

    NTT-based Modular Polynomial Multiplication with f(x)=x^n-1
    -- A(x),B(x): n-1 degree polynomials
    -- C(x)     : n-1 degree polynomial
    -- C(x)=A(x)*B(X) --> C=INTT_n(NTT_n(A) . NTT_n(B))

    NTT-based Modular Polynomial Multiplication with f(x)=x^n+1
    -- First implementation
    ---- A(x),B(x): n-1 degree polynomials
    ---- C(x)     : n-1 degree polynomial
    ---- C(x)=A(x)*B(X) --> C=PostProc(INTT_n(NTT_n(PreProc(A)) . NTT_n(PreProc(B))))
    -- Second implementation
    ---- A(x),B(x): n-1 degree polynomials
    ---- C(x)     : n-1 degree polynomial
    ---- C(x)=A(x)*B(X) --> C=INTT_n(MergedNTT_n(A) . MergedNTT_n(B))
    -- First method uses separate pre- and post-processing methods
    -- The second method merges pre- and post-processing with NTT and INTT, respectively

    CRT-based Modular Polynomial Multiplication
    -- f(x)=x^n-1
    -- f(x)=x^n+1 (with w)
    -- f(x)=x^n+1 (with psi)
    -- f(x)=x^n-x^(n/2)+1

    CRT-based Unified Modular Polynomial Multiplication
    -- with any final degree and any f(x)

    Other methods (not implemented in this class)
    -- Karatsuba
    -- Toom-Cook
    -- Schonhage-Strassen

    *******************************
    """

    # A,B: same degree polynomials
    # q: coefficient modulus
    # C: output polynomial
    def SchoolbookPolMul(self,A, B, q):
        C = [0] * (2 * len(A))
        for indexA, elemA in enumerate(A):
            for indexB, elemB in enumerate(B):
                C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % q
        return C

    # A: input polynomial
    # F: reduction polynomial
    # q: coefficient modulus
    # D: output polynomial
    # Assuming coefficient of largest degree of F is 1
    def PolRed(self,A,F,q):
        if len(A) < len(F):
            return A
        else:
            D = [_ for _ in A]
            R = [(-x) % q for x in F[0:len(F)-1]]
            for i in range(len(D)-1,len(F)-2,-1):
                for j in range(len(R)):
                    D[i-1-j] = (D[i-1-j] + D[i]*R[len(R)-1-j]) % q
                D[i] = 0
            return D[0:len(F)-1]

    # A,B: input polynomials in x^n-1
    # q: coefficient modulus
    # D: output polynomial in x^n-1
    def SchoolbookModPolMul_PWC(self,A, B, q):
        C = [0] * (2 * len(A))
        D = [0] * (len(A))
        for indexA, elemA in enumerate(A):
            for indexB, elemB in enumerate(B):
                C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % q

        for i in range(len(A)):
            D[i] = (C[i] + C[i + len(A)]) % q
        return D

    # A,B: input polynomials in x^n+1
    # q: coefficient modulus
    # D: output polynomial in x^n+1
    def SchoolbookModPolMul_NWC(self,A, B, q):
        C = [0] * (2 * len(A))
        D = [0] * (len(A))
        for indexA, elemA in enumerate(A):
            for indexB, elemB in enumerate(B):
                C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % q

        for i in range(len(A)):
            D[i] = (C[i] - C[i + len(A)]) % q
        return D

    # A,B: input polynomials in x^n-x^(n/2)+1 where n = 3*2^t
    # q: coefficient modulus
    # D: output polynomial in x^n-x^(n/2)+1 where n = 3*2^t
    def SchoolbookModPolMul_NTRU(self,A, B, q):
        C = [0] * (2 * len(A))
        for indexA, elemA in enumerate(A):
            for indexB, elemB in enumerate(B):
                C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % q
        D = [_ for _ in C]
        for i in range(2*len(A)-1,len(A)-1,-1):
            D[i-len(A)+int(len(A)/2)] = (D[i-len(A)+int(len(A)/2)] + D[i]) % q
            D[i-len(A)              ] = (D[i-len(A)              ] - D[i]) % q
            D[i] = 0
        return D[0:len(A)]

    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: 2n-2 degree polynomial
    # C(x)=A(x)*B(X) --> C=INTT_2n(NTT_2n(zero-padded A) . NTT_2n(zero-padded B))
    def NTTBasedPolMul(self,A,B,w,w_inv,q):
        N = len(A)
        A_padded = A + [0]*N
        B_padded = B + [0]*N

        Evaluator = NTT()

        A_ntt = Evaluator.Radix2_DIT_Iterative_NTT_NR(A_padded,w,q)
        B_ntt = Evaluator.Radix2_DIT_Iterative_NTT_NR(B_padded,w,q)

        C_ntt = [(x*y) % q for x,y in zip(A_ntt,B_ntt)]

        C = Evaluator.Radix2_DIF_Iterative_INTT_RN(C_ntt,w_inv,q)
        return C

    # NTT-Based polynomial multiplication with carrier modulus
    # This can be implemented with other techniques we used so far
    # for k-bit modulus q, we need a log(n)+2*k-bit carrier modulus Q
    # After multiplication operation, separate polynomial reduction is required
    def NTTBasedPolMulwithCM(self,A,B,qw,qw_inv,q,Qw,Qw_inv,Q):
        """
        -- Polynomial reduction operation is performed separately
        -- If n is a power-of-two
        ---- Select a new ntt-friendly q with bit-size > log(n*q^2)
        ---- Perform "NTT-based Polynomial Multiplication" and then apply reduction
        -- If n is NOT a power-of-two
        ---- Zero-pad input polynomials to the closest power-of-two
        ---- Select a new ntt-friendly q with bit-size > log(n*q^2)
        ---- Perform "NTT-based Polynomial Multiplication" and then apply reduction
        """
        pass

    # NTT-Based Modular Polynomial Multiplication with f(x)=x^n-1 (Positive Wrapped Convolution)
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    # C(x)=A(x)*B(X) --> C=INTT_n(NTT_n(A) . NTT_n(B))
    def NTTBasedModPolMul_PWC(self,A,B,w,w_inv,q):
        Evaluator = NTT()

        A_ntt = Evaluator.Radix2_DIT_Iterative_NTT_NR(A,w,q)
        B_ntt = Evaluator.Radix2_DIT_Iterative_NTT_NR(B,w,q)

        C_ntt = [(x*y) % q for x,y in zip(A_ntt,B_ntt)]

        C = Evaluator.Radix2_DIF_Iterative_INTT_RN(C_ntt,w_inv,q)
        return C

    # NTT-Based Modular Polynomial Multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- with separate pre-processing and post-processing
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    # C(x)=A(x)*B(X) --> C=PostProc(INTT_n(NTT_n(PreProc(A)) . NTT_n(PreProc(B))))
    def NTTBasedModPolMul_NWC_v1(self,A,B,w,w_inv,psi,psi_inv,q):
        Evaluator = NTT()

        A_p = [(x*pow(psi,i,q)) % q for i,x in enumerate(A)]
        B_p = [(x*pow(psi,i,q)) % q for i,x in enumerate(B)]

        A_ntt = Evaluator.Radix2_DIT_Iterative_NTT_NR(A_p,w,q)
        B_ntt = Evaluator.Radix2_DIT_Iterative_NTT_NR(B_p,w,q)

        C_ntt = [(x*y) % q for x,y in zip(A_ntt,B_ntt)]

        C_p = Evaluator.Radix2_DIF_Iterative_INTT_RN(C_ntt,w_inv,q)

        C = [(x*pow(psi_inv,i,q)) % q for i,x in enumerate(C_p)]
        return C

    # NTT-Based Modular Polynomial Multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- with merged pre-processing and post-processing
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    # C(x)=A(x)*B(X) --> C=INTT_n(MergedNTT_n(A) . MergedNTT_n(B))
    def NTTBasedModPolMul_NWC_v2(self,A,B,psi,psi_inv,q):

        A_ntt = self.CTBasedMergedNTT_NR(A,psi,q)
        B_ntt = self.CTBasedMergedNTT_NR(B,psi,q)

        C_ntt = [(x*y) % q for x,y in zip(A_ntt,B_ntt)]

        C = self.GSBasedMergedINTT_RN(C_ntt,psi_inv,q)
        return C

    # CRT-based modular polynomial multiplication with f(x)=x^n-1 (Positive Wrapped Convolution)
    # It is using CRT-based NTT instead of regular NTT
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    # C(x)=A(x)*B(X) --> C=INTT_n(NTT_n(A) . NTT_n(B))
    def CRTBasedModPolMul_PWC(self,A,B,w,w_inv,q):
        # Note: if you use CRT_Recursive_NTT, output of NTT operation will be in bit-reversed order
        Evaluator = NTT()

        A_ntt = Evaluator.CRT_Full_NTT(A,w,q)
        B_ntt = Evaluator.CRT_Full_NTT(B,w,q)

        C_ntt = [(x*y) % q for x,y in zip(A_ntt,B_ntt)]

        C = Evaluator.CRT_Full_INTT(C_ntt,w_inv,q)
        return C

    # CRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing CRT instead of NTT
    # -- it is using psi (q = 1 mod 2n) and final degree of CRT reduction is 1
    # -- It is same as MergedNTT method (as shown in CRTBasedModPolMul_PWC)
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    def CRTBasedModPolMul_NWC_FD1(self,A,B,psi,psi_inv,q):
        A_ntt = self.CRT_Iterative_NWC_FD1(A,psi,q)
        B_ntt = self.CRT_Iterative_NWC_FD1(B,psi,q)

        C_ntt = [(x*y) % q for x,y in zip(A_ntt,B_ntt)]

        C = self.ICRT_Iterative_NWC_FD1(C_ntt,psi_inv,q)
        return C

    # CRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing CRT instead of NTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is 2
    # -- it is using Iterative version of reduction function
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    def CRTBasedModPolMul_NWC_FD2(self,A,B,w,w_inv,q):
        A_ntt = self.CRT_Iterative_NWC_FD2_NR(A,w,q)
        B_ntt = self.CRT_Iterative_NWC_FD2_NR(B,w,q)

        C_ntt = [0 for _ in range(len(A))]

        # Degree-2 modular polynomial multiplications
        for i in range(len(A)//2):
            w_pow = 2*intReverse(i,int(log(len(A)//2,2)))+1
            wk    = pow(w,w_pow,q)
            C_ntt[2*i:2*i+2] = self.PolWiseMult(A_ntt[2*i:2*i+2],B_ntt[2*i:2*i+2],wk,2,q)

        # NOTE: it is using w. (We cen convert it into w_inv by modification)
        C = self.ICRT_Iterative_NWC_FD2_RN(C_ntt,w_inv,q)

        return C

    # CRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing CRT instead of NTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is a variable (power-of-two, 2, 4, ...)
    # -- it is using Iterative version of reduction function
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    def CRTBasedModPolMul_NWC_FDV(self,A,B,w,w_inv,q,findeg):
        A_ntt = self.CRT_Iterative_NWC_FDV_NR(A,w,q,findeg)
        B_ntt = self.CRT_Iterative_NWC_FDV_NR(B,w,q,findeg)

        C_ntt = [0 for _ in range(len(A))]

        # Degree-findeg modular polynomial multiplications
        for i in range(len(A)//findeg):
            w_pow = 2*intReverse(i,int(log(len(A)//findeg,2)))+1
            wk    = pow(w,w_pow,q)
            C_ntt[findeg*i:findeg*i+findeg] = self.PolWiseMult(A_ntt[findeg*i:findeg*i+findeg],B_ntt[findeg*i:findeg*i+findeg],wk,findeg,q)

        C = self.ICRT_Iterative_NWC_FDV_RN(C_ntt,w_inv,q,findeg)

        return C

    # CRT-based modular polynomial multiplication with f(x)=x^n-x^(n/2)+1 (Negative Wrapped Convolution)
    # -- utilizing CRT instead of NTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is 3
    # -- it is using Iterative version of reduction function
    # A,B: n-1 degree polynomials
    # w, w_inv: twiddle factors
    # q: coefficient modulus
    # C: n-1 degree polynomial
    def CRTBasedModPolMul_NTRU_FD3(self,A,B,w,w_inv,ntrupowersf,ntrupowersb,ntrupowersi,q):
        # Initial reduction
        A_r = [_ for _ in A]
        B_r = [_ for _ in B]

        wk = pow(w,len(A)//6,q)
        for i in range(len(A)//2):
            t1 = (wk*A[i+len(A)//2]) % q
            t2 = (wk*B[i+len(B)//2]) % q
            A_r[i+len(A)//2] = (A[i]+A[i+len(A)//2]-t1)%q
            A_r[i]           = (A[i]               +t1)%q
            B_r[i+len(B)//2] = (B[i]+B[i+len(B)//2]-t2)%q
            B_r[i]           = (B[i]               +t2)%q

        # NTT
        A_ntt = self.CRT_Iterative_NTRU_FD3_NR(A_r,w,ntrupowersf,q)
        B_ntt = self.CRT_Iterative_NTRU_FD3_NR(B_r,w,ntrupowersf,q)

        C_ntt = [0 for _ in range(len(A))]

        # Degree-findeg modular polynomial multiplications
        for i in range(len(A)//3):
            w_pow = ntrupowersb[i]
            wk    = pow(w,w_pow,q)
            C_ntt[3*i:3*i+3] = self.PolWiseMult(A_ntt[3*i:3*i+3],B_ntt[3*i:3*i+3],wk,3,q)

        # INTT
        C = self.ICRT_Iterative_NTRU_FD3_RN(C_ntt,w_inv,ntrupowersi,q)

        # Final reconstruction
        C_r = [_ for _ in C]

        wk = modinv((2*pow(w,len(C)//6,q)-1)%q,q)

        for i in range(len(C)//2):
            t = ((C[i]-C[i+len(C)//2])*wk)%q   # t = f[i+N//2]
            C_r[i          ] = (C[i]+C[i+len(C)//2]-t)%q
            C_r[i+len(C)//2] = (2*t)%q

        return C_r


    def CRTBasedModPolMul_NTRU_FDV(self,A,B,w,w_inv,ntrupowersf,ntrupowersb,ntrupowersi,q,findeg):
        # Initial reduction
        A_r = [_ for _ in A]
        B_r = [_ for _ in B]

        wk = pow(w,(len(A)//6),q)
        # wk = pow(w,(len(A)//6)//(findeg//3),q)

        for i in range(len(A)//2):
            t1 = (wk*A[i+len(A)//2]) % q
            t2 = (wk*B[i+len(B)//2]) % q
            A_r[i+len(A)//2] = (A[i]+A[i+len(A)//2]-t1)%q
            A_r[i]           = (A[i]               +t1)%q
            B_r[i+len(B)//2] = (B[i]+B[i+len(B)//2]-t2)%q
            B_r[i]           = (B[i]               +t2)%q

        # NTT
        A_ntt = self.CRT_Iterative_NTRU_FDV_NR(A_r,w,ntrupowersf,q,findeg)
        B_ntt = self.CRT_Iterative_NTRU_FDV_NR(B_r,w,ntrupowersf,q,findeg)

        C_ntt = [0 for _ in range(len(A))]

        # Degree-findeg modular polynomial multiplications
        for i in range(len(A)//findeg):
            w_pow = ((findeg//3)*ntrupowersb[i*(findeg//3)])
            # w_pow = ((findeg//3)*ntrupowersb[i*(findeg//3)])//(findeg//3)

            wk    = pow(w,w_pow,q)
            C_ntt[findeg*i:findeg*i+findeg] = self.PolWiseMult(A_ntt[findeg*i:findeg*i+findeg],B_ntt[findeg*i:findeg*i+findeg],wk,findeg,q)

        # INTT
        C = self.ICRT_Iterative_NTRU_FDV_RN(C_ntt,w_inv,ntrupowersi,q,findeg)

        # Final reconstruction
        C_r = [_ for _ in C]

        wk = modinv((2*pow(w,(len(C)//6),q)-1)%q,q)
        # wk = modinv((2*pow(w,(len(C)//6)//(findeg//3),q)-1)%q,q)

        for i in range(len(C)//2):
            t = ((C[i]-C[i+len(C)//2])*wk)%q   # t = f[i+N//2]
            C_r[i          ] = (C[i]+C[i+len(C)//2]-t)%q
            C_r[i+len(C)//2] = (2*t)%q

        return C_r

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
    def CRTBasedModPolMul_Unified(self,A,B,w,w_inv,q,ring,findeg,ntrupowersf=[],ntrupowersb=[],ntrupowersi=[]):
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

            wk = pow(w,len(A)//6,q)
            for i in range(len(A)//2):
                t1 = (wk*A[i+len(A)//2]) % q
                t2 = (wk*B[i+len(B)//2]) % q
                A_r[i+len(A)//2] = (A[i]+A[i+len(A)//2]-t1)%q
                A_r[i]           = (A[i]               +t1)%q
                B_r[i+len(B)//2] = (B[i]+B[i+len(B)//2]-t2)%q
                B_r[i]           = (B[i]               +t2)%q

        # --------------------------------------------- NTT
        A_ntt = self.CRT_Iterative_Unified_NR(A_r,w,q,ring,findeg,ntrupowersf)
        B_ntt = self.CRT_Iterative_Unified_NR(B_r,w,q,ring,findeg,ntrupowersf)

        # --------------------------------------------- Degree-findeg modular polynomial multiplications
        C_ntt = [0 for _ in range(len(A))]
        for i in range(len(A)//findeg):
            if ring == 0:
                # NWC
                w_pow = 2*intReverse(i,int(log(len(A)//findeg,2)))+1
            else:
                # NTRU
                w_pow = ((findeg//3)*ntrupowersb[i*(findeg//3)])
                # w_pow = ((findeg//3)*ntrupowersb[i*(findeg//3)])//(findeg//3)

            wk    = pow(w,w_pow,q)
            C_ntt[findeg*i:findeg*i+findeg] = self.PolWiseMult(A_ntt[findeg*i:findeg*i+findeg],B_ntt[findeg*i:findeg*i+findeg],wk,findeg,q)

        # --------------------------------------------- INTT
        C = self.ICRT_Iterative_Unified_RN(C_ntt,w_inv,q,ring,findeg,ntrupowersi)

        # --------------------------------------------- Final step
        if ring == 0:
            """
            NWC requires no final reconstruction step
            """
            return C
        else:
            """
            NTRU requires final reconstruction step
            """
            wk = modinv((2*pow(w,(len(C)//6),q)-1)%q,q)
            # wk = modinv((2*pow(w,(len(C)//6)//(findeg//3),q)-1)%q,q)

            for i in range(len(C)//2):
                t = ((C[i]-C[i+len(C)//2])*wk)%q   # t = f[i+N//2]
                C[i          ] = (C[i]+C[i+len(C)//2]-t)%q
                C[i+len(C)//2] = (2*t)%q

            return C

    ############################################################################################ (Helper Function for Pol Mul Operations)

    # Multiplies two "deg" degree polynomial in x^"deg"-w^k where k is some power
    # A,B: input polynomials
    # wk: w^k
    # deg: degree
    # q: coefficient modulus
    # C: output polynomial
    def PolWiseMult(self,A,B,wk,deg,q):
        C = [0] * ((2 * deg)-1)
        # D = [0] * ((2 * deg)-1)

        if deg == 1:
            # if final degree is 1
            D = [(x*y)%q for x,y in zip(A,B)]
            return D[0:deg]
        else:
            # if final degree is larger than 1
            for indexA, elemA in enumerate(A):
                for indexB, elemB in enumerate(B):
                    C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % q

            D = [_ for _ in C]
            for i in range(len(A)-1):
                D[i] = (C[i] + C[i + len(A)]*wk) % q

        return D[0:deg]

    # -------------------------------------------------------------------------- Iterative

    # Merged NTT with pre-processing (optimized) (iterative)
    # This is not NTT, this is pre-processing + NTT
    # (see: https://eprint.iacr.org/2016/504.pdf)
    # A: input polynomial (standard order)
    # Psi: 2n-th root of unity
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CTBasedMergedNTT_NR(self,A,Psi,q):
        N = len(A)
        B = [_ for _ in A]

        l = int(log(N,2))

        t = N
        m = 1
        while(m<N):
            t = int(t/2)
            for i in range(m):
                j1 = 2*i*t
                j2 = j1 + t - 1
                Psi_pow = intReverse(m+i,l)
                S = pow(Psi,Psi_pow,q)
                for j in range(j1,j2+1):
                    U = B[j]
                    V = (B[j+t]*S) % q

                    B[j]   = (U+V) % q
                    B[j+t] = (U-V) % q
            m = 2*m

        return B

    # Merged INTT with post-processing (optimized) (iterative)
    # This is not NTT, this is pre-processing + NTT
    # (see: https://eprint.iacr.org/2016/504.pdf)
    # A: input polynomial (Bit-reversed order)
    # Psi: 2n-th root of unity
    # q: modulus
    # B: output polynomial (standard order)
    def GSBasedMergedINTT_RN(self,A,Psi,q):
        N = len(A)
        B = [_ for _ in A]

        l = int(log(N,2))

        t = 1
        m = N
        while(m>1):
            j1 = 0
            h = int(m/2)
            for i in range(h):
                j2 = j1 + t - 1
                Psi_pow = intReverse(h+i,l)
                S = pow(Psi,Psi_pow,q)
                for j in range(j1,j2+1):
                    U = B[j]
                    V = B[j+t]

                    B[j]   = (U+V) % q
                    B[j+t] = (U-V)*S % q
                j1 = j1 + 2*t
            t = 2*t
            m = int(m/2)

        N_inv = modinv(N, q)
        for i in range(N):
            B[i] = (B[i] * N_inv) % q

        return B

    # CRT-based Merged NTT for f(x)=x^n+1 (Negative Wrapped Convolution) (iterative)
    # Actually, same as CTBasedMergedNTT_NR() function
    # This is not NTT, this is pre-processing + NTT
    # A: input polynomial (standard order)
    # Psi: 2n-th root of unity
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CRT_Iterative_NWC_FD1(self,A,psi,q):
        return self.CTBasedMergedNTT_NR(A,psi,q)

    # ICRT-based Merged INTT for f(x)=x^n+1 (Negative Wrapped Convolution) (iterative)
    # Actually, same as GSBasedMergedINTT_RN() function
    # This is not NTT, this is pre-processing + NTT
    # A: input polynomial (bit-reversed order)
    # Psi: 2n-th root of unity
    # q: modulus
    # B: output polynomial (standard order)
    def ICRT_Iterative_NWC_FD1(self,A,psi,q):
        return self.GSBasedMergedINTT_RN(A,psi,q)

    # CRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing CRT instead of NTT (Iterative Version)
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is 2
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    # NOTE: 1 iteration less executed version of "CTBasedMergedNTT_NR"
    def CRT_Iterative_NWC_FD2_NR(self,A,w,q):
        N = len(A)
        B = [_ for _ in A]

        k=1
        lena = (N//2)

        v = int(log(lena,2))

        while lena >= 2:
            start = 0
            while start < N:
                W_pow = intReverse(k,v)
                W = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = (W * B[j+lena]) % q

                    B[j+lena] = (B[j] - t) % q
                    B[j     ] = (B[j] + t) % q

                    j = j+1
                start = j + lena
            lena = (lena//2)

        return B

    # ICRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing ICRT instead of INTT (Iterative Version)
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is 2
    # A: input polynomial (bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    # NOTE: 1 iteration less executed version of "GSBasedMergedINTT_NR"
    def ICRT_Iterative_NWC_FD2_RN(self,A,w,q):
        N = len(A)
        B = [_ for _ in A]

        k = 0
        lena = 2

        v = int(log(N//2,2))

        while lena <= (N//2):
            start = 0
            while start < N:
                W_pow = intReverse(k,v)+1
                TW = pow(w,W_pow,q)
                """
                W_pow and TW below use "w" instead of "w_inv"

                W_pow = (N//2) - 1 - intReverse(k,v)
                TW = (-pow(w,W_pow,q)) % q # here, "-" means an extra w^(n/2)
                """
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = B[j]

                    B[j       ] = (t + B[j + lena]) % q
                    B[j + lena] = (t - B[j + lena]) % q
                    B[j + lena] = B[j + lena]*TW % q

                    j = j+1
                start = j + lena
            lena = 2*lena

        N_inv = modinv(N//2,q)
        for i in range(N):
            B[i] = (B[i] * N_inv) % q

        return B

    # CRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing CRT instead of NTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is a variable (power-of-two, 2, 4, ...)
    # -- it is using Iterative version of reduction function
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CRT_Iterative_NWC_FDV_NR(self,A,w,q,findeg):
        N = len(A)
        B = [_ for _ in A]

        k=1
        lena = (N//2)

        v = int(log(N//findeg,2))

        while lena >= findeg:
            start = 0
            while start < N:
                W_pow = intReverse(k,v)
                W = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = (W * B[j+lena]) % q

                    B[j+lena] = (B[j] - t) % q
                    B[j     ] = (B[j] + t) % q

                    j = j+1
                start = j + lena
            lena = (lena//2)

        return B

    # ICRT-based modular polynomial multiplication with f(x)=x^n+1 (Negative Wrapped Convolution)
    # -- utilizing ICRT instead of INTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is a variable (power-of-two, 2, 4, ...)
    # -- it is using Iterative version of reduction function
    # A: input polynomial (bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def ICRT_Iterative_NWC_FDV_RN(self,A,w,q,findeg):
        N = len(A)
        B = [_ for _ in A]

        k = 0
        lena = findeg

        v = int(log(N//findeg,2))

        while lena <= (N//2):
            start = 0
            while start < N:
                """
                W_pow = intReverse(k,v)+1
                TW = (-pow(w,W_pow,q)) % q # here, "-" means an extra w^(n/2)
                """
                W_pow = intReverse(k,v)+1
                TW = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = B[j]

                    B[j       ] = (t + B[j + lena]) % q
                    B[j + lena] = (t - B[j + lena]) % q
                    B[j + lena] = B[j + lena]*TW % q

                    j = j+1
                start = j + lena
            lena = 2*lena

        N_inv = modinv(N//findeg,q)
        for i in range(N):
            B[i] = (B[i] * N_inv) % q

        return B

    # CRT-based modular polynomial multiplication with f(x)=x^n-x^n/2+1 (NTRU)
    # -- utilizing CRT instead of NTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is 3
    # -- it is using Iterative version of reduction function
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CRT_Iterative_NTRU_FD3_NR(self,A,w,powers,q):
        N = len(A)
        B = [_ for _ in A]

        k=0
        lena = (N//4)

        while lena >= 3:
            start = 0
            while start < N:
                W_pow = powers[k]
                W = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = (W * B[j+lena]) % q

                    B[j+lena] = (B[j] - t) % q
                    B[j     ] = (B[j] + t) % q

                    j = j+1
                start = j + lena
            lena = (lena//2)

        return B

    # ICRT-based modular polynomial multiplication with f(x)=x^n-x^n/2+1 (NTRU)
    # -- utilizing ICRT instead of INTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is 3
    # -- it is using Iterative version of reduction function
    # A: input polynomial (bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def ICRT_Iterative_NTRU_FD3_RN(self,A,w,powers,q):
        N = len(A)
        B = [_ for _ in A]

        k = 0
        lena = 3

        v = int(log(N//3,2))

        while lena <= (N//4):
            start = 0
            while start < N:
                W_pow = powers[k]
                TW = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = B[j]

                    B[j       ] = (t + B[j + lena]) % q
                    B[j + lena] = (t - B[j + lena]) % q
                    B[j + lena] = B[j + lena]*TW % q

                    j = j+1
                start = j + lena
            lena = 2*lena

        N_inv = modinv(N//3,q)
        for i in range(N):
            B[i] = (B[i] * N_inv) % q

        return B

    # CRT-based modular polynomial multiplication with f(x)=x^n-x^n/2+1 (NTRU)
    # -- utilizing CRT instead of NTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is a variable (3,6,12,...)
    # -- it is using Iterative version of reduction function
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CRT_Iterative_NTRU_FDV_NR(self,A,w,powers,q,findeg):
        N = len(A)
        B = [_ for _ in A]

        k=0
        lena = (N//4)

        while lena >= findeg:
            start = 0
            while start < N:
                W_pow = powers[k]
                # W_pow = (powers[k] // (findeg//3))

                W = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = (W * B[j+lena]) % q

                    B[j+lena] = (B[j] - t) % q
                    B[j     ] = (B[j] + t) % q

                    j = j+1
                start = j + lena
            lena = (lena//2)

        return B

    # ICRT-based modular polynomial multiplication with f(x)=x^n-x^n/2+1 (NTRU)
    # -- utilizing ICRT instead of INTT
    # -- it is using w (q = 1 mod n) and final degree of CRT reduction is a variable (3,6,12,...)
    # -- it is using Iterative version of reduction function
    # A: input polynomial (bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def ICRT_Iterative_NTRU_FDV_RN(self,A,w,powers,q,findeg):
        N = len(A)
        B = [_ for _ in A]

        k = 0
        lena = findeg

        # Powers need to be adjusted accordingly
        powers_new = [_ for _ in powers]
        i = findeg
        r = 1
        while(i >= 6):
            powers_new = powers_new[N//(6*r):]
            i = (i//2)
            r = 2*r

        while lena <= (N//4):
            start = 0
            while start < N:
                W_pow = powers_new[k]
                # W_pow = (powers_new[k] // (findeg//3))

                TW = pow(w,W_pow,q)
                k = k+1
                j = start
                while(j < (start + lena)):
                    t = B[j]

                    B[j       ] = (t + B[j + lena]) % q
                    B[j + lena] = (t - B[j + lena]) % q
                    B[j + lena] = B[j + lena]*TW % q

                    j = j+1
                start = j + lena
            lena = 2*lena

        N_inv = modinv(N//findeg,q)
        for i in range(N):
            B[i] = (B[i] * N_inv) % q

        return B

    # CRT-based unified structure for NWC and NTRU
    # -- ring  : 0->NWC, 1->NTRU
    # -- findeg: 1,2,4,... for NWC and 3,6,12,... for NTRU
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CRT_Iterative_Unified_NR(self,A,w,q,ring,findeg,powers):
        N = len(A)
        B = [_ for _ in A]

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

                    j = j+1
                start = j + lena
            lena = (lena//2)

        return B

    # ICRT-based unified structure for NWC and NTRU
    # -- ring  : 0->NWC, 1->NTRU
    # -- findeg: 1,2,4,... for NWC and 3,6,12,... for NTRU
    # A: input polynomial (bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def ICRT_Iterative_Unified_RN(self,A,w,q,ring,findeg,powers):
        N = len(A)
        B = [_ for _ in A]

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

                    j = j+1
                start = j + lena
            lena = 2*lena

        N_inv = modinv(N//findeg,q)
        for i in range(N):
            B[i] = (B[i] * N_inv) % q

        return B

#
