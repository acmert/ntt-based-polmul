from math import log,floor,ceil
from helper import *

matrix = lambda polynomial, col_length: list(zip(*[polynomial[i:i + col_length] for i in range(0, len(polynomial), col_length)]))

# Cooley-Tukey Butterfly Structure
# A0,A1: input coefficients
# W: twiddle factor
# q: modulus
# B0,B1: output coefficients
def CT_Butterfly(A0,A1,W,q):
    """
    A0 -------\--|+|-- B0
               \/
               /\
    A1 --|x|--/--|-|-- B1
    """
    M = (A1 * W) % q

    B0 = (A0 + M) % q
    B1 = (A0 - M) % q

    return B0,B1

# Gentleman-Sandle Butterfly Structure
# A0,A1: input coefficients
# W: twiddle factor
# q: modulus
# B0,B1: output coefficients
def GS_Butterfly(A0,A1,W,q):
    """
    A0 --\--|+|------- B0
          \/
          /\
    A1 --/--|-|--|x|-- B1
    """
    M0 = (A0 + A1) % q
    M1 = (A0 - A1) % q

    B0 = M0
    B1 = (M1 * W) % q

    return B0,B1

class NTT:
    """
    - Definition of NTT:

    Existence condition: q = 1 (mod n) and w: n-th root of unity

    [a_0, a_1, ..., a_n-1] --> [A_0, A_1, ..., A_n-1]

    Forward NTT: A_i = sum{j from 0 to n-1}(a_j * w^ij mod q) for i from 0 to n-1
    Inverse NTT: a_i = sum{j from 0 to n-1}(A_j * w^-ij mod q) for i from 0 to n-1
    """

    """
    List of NTT Algorithms: (Inside the FFT Black Box, by Chu and George)
    -- Naive NTT (see Wikipedia definition of NTT operation)
    -- Radix-2 Decimation-in-Time (DIT) Recursive NTT (Cooley-Tukey)
    -- Radix-2 Decimation-in-Frequency (DIF) Recursive NTT (Gentleman-Sandle)
    -- Radix-2 Decimation-in-Time (DIT) Iterative NTT
    ---- NR (N: Natural order, R: Reversed Order)
    ---- RN
    ---- NN
    -- Radix-2 Decimation-in-Time (DIF) Iterative NTT
    ---- NR
    ---- RN
    ---- NN
    """

    """
    Note: Any forward NTT function can be used for inverse NTT if you give input
    in proper order and w^-1 instead of w. Besides, INTT requires output
    coefficients to be multiplied with n^-1 mod q.
    """

    """
    - What is standard order? : 0, 1, ..., n-1
    - What is reversed/bit-reversed (scrambled) order? : 0, br(1), ..., br(n-1)
    where br() function bit-revese the integer with log(n) bits
    """

    # Naive NTT
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def NaiveNTT_NN(self,A,W,q):
        """
        Very slow baseline implementation. Do not use for large parameter set.
        """
        N = len(A)
        B = [0]*N

        for i in range(N):
            for j in range(N):
                B[i] = (B[i] + A[j]*(W**(i*j))) % q

        return B

    # Naive NTT
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def NaiveNTT_NR(self,A,W,q):
        """
        Very slow baseline implementation. Do not use for large parameter set.
        """
        N = len(A)
        B = [0]*N

        v = int(log(N,2))

        for i in range(N):
            for j in range(N):
                W_pow = intReverse(i,v)*j
                B[i] = (B[i] + A[j]*(W**W_pow)) % q
        return B

    # Recursive Radix-2 Decimation-in-Time (DIT) (CT) NTT
    # A: input polynomial (standard order --> it becomes reversed after recursions)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def Radix2_DIT_Recursive_NTT(self,A,W,q):
        N = len(A)
        B = [0]*N

        if N == 2:
            # When N is 2, perform butterfly operation with W = 1
            B[0] = (A[0] + A[1]) % q
            B[1] = (A[0] - A[1]) % q

            return B
        else:
            # Divide input into two (even indices, odd indices)
            AE = [A[i] for i in range(0,N,2)]
            AO = [A[i] for i in range(1,N,2)]

            # Reduce twiddle factor for the next recursion
            W_new = pow(W,2,q)

            # Apply NTT operations to the even and odd indices of the input separately
            BE = self.Radix2_DIT_Recursive_NTT(AE,W_new,q)
            BO = self.Radix2_DIT_Recursive_NTT(AO,W_new,q)

            # Outputs of first and second NTT operations go to the first and second
            # half of the array (output array)
            B = BE+BO

            # Perform CT-Butterfly where first and second inputs of butterfly
            # operation are from first and second half of the output respectively
            # First and second outputs of the butterfly operation go to first and
            # second half of the array (input array) respectively
            for i in range(int(N/2)):
                B[i], B[i+int(N/2)] = CT_Butterfly(B[i],B[i+int(N/2)],pow(W,i,q),q)

            return B

    # Recursive Radix-2 Decimation-in-Frequency (DIF) (GS) NTT
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def Radix2_DIF_Recursive_NTT(self,A,W,q):
        N = len(A)
        B = [0]*N

        if N == 2:
            # When N is 2, perform butterfly operation with W = 1
            B[0] = (A[0] + A[1]) % q
            B[1] = (A[0] - A[1]) % q

            return B
        else:
            # Divide input into two (first half, second half)

            # Perform GS-Butterfly where first and second inputs of butterfly
            # operation are from first and second half of the input respectively
            # First and second outputs of the butterfly operation go to first and
            # second half of the array (input array) respectively
            for i in range(int(N/2)):
                B[i], B[i+int(N/2)] = GS_Butterfly(A[i],A[i+int(N/2)],pow(W,i,q),q)

            # Reduce twiddle factor for the next recursion
            W_new = pow(W,2,q)

            # Apply NTT operations to the first and second half of the input separately
            BE = self.Radix2_DIF_Recursive_NTT(B[0:int(N/2)],W_new,q)
            BO = self.Radix2_DIF_Recursive_NTT(B[int(N/2):N],W_new,q)

            # Outputs of first and second NTT operations go to the first and second
            # half of the array (output array)
            B = BE+BO

            return B

    # From paper: NTTU: An Area-Efficient Low-POwer NTT-Uncoupled Architecture for NTT-Based Multiplication
    # Iterative Radix-2 Decimation-in-Time (DIT) (CT) NTT - NR
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def Radix2_DIT_Iterative_NTT_NR(self,A,W,q):
        N = len(A)
        B = [_ for _ in A]

        for s in range(int(log(N,2)),0,-1):
            m = 2**s
            for k in range(int(N/m)):
                TW = pow(W,intReverse(k,int(log(N,2))-s)*int(m/2),q)
                for j in range(int(m/2)):
                    u = B[k*m+j]
                    t = (TW*B[k*m+j+int(m/2)]) % q

                    B[k*m+j]          = (u+t) % q
                    B[k*m+j+int(m/2)] = (u-t) % q

        return B

    # Iterative Radix-2 Decimation-in-Time (DIT) (CT) NTT - RN
    # A: input polynomial (bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def Radix2_DIT_Iterative_NTT_RN(self,A,W,q):
        N = len(A)
        B = [_ for _ in A]

        v = int(N/2)
        m = 1
        d = 1

        while m<N:
            np = 2*m
            lp = np*(v-1)
            for k in range(m):
                j = k
                jl = k + lp
                jt = k*v
                TW = pow(W,jt,q)
                while j < (jl+1):
                    temp = (TW*B[j+d]) % q
                    B[j+d] = (B[j] - temp) % q
                    B[j]   = (B[j] + temp) % q
                    j = j+np
            v = int(v/2)
            m = 2*m
            d = 2*d

        return B

    # Iterative Radix-2 Decimation-in-Time (DIT) (CT) NTT - NN
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def Radix2_DIT_Iterative_NTT_NN(self,A,W,q):
        N = len(A)
        B = [_ for _ in A]
        C = [_ for _ in A]
        # C = [0]*N

        v = int(N/2)
        m = 1
        d = int(N/2)

        if int(log(v))%2 == 0:
            nsi = True
        else:
            nsi = False

        while m<N:
            if nsi:
                l = 0
                for k in range(m):
                    jf = 2*k*v
                    jl = jf + v - 1
                    jt = k*v

                    TW = pow(W,jt,q)

                    for j in range(jf,jl+1):
                        temp = (TW*B[j+d]) % q

                        C[l]          = (B[j] + temp) % q
                        C[l+int(N/2)] = (B[j] - temp) % q

                        l = l+1
                nsi = False
            else:
                l = 0
                for k in range(m):
                    jf = 2*k*v
                    jl = jf + v - 1
                    jt = k*v

                    TW = pow(W,jt,q)

                    for j in range(jf,jl+1):
                        temp = (TW*C[j+d]) % q

                        B[l]          = (C[j] + temp) % q
                        B[l+int(N/2)] = (C[j] - temp) % q

                        l = l+1
                nsi = True
            v = int(v/2)
            m = 2*m
            d = int(d/2)

        return C

    # Iterative Radix-2 Decimation-in-Frequency (DIF) (GS) NTT - NR
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def Radix2_DIF_Iterative_NTT_NR(self,A,W,q):
        N = len(A)
        B = [_ for _ in A]

        m = 1
        v = N

        while v>1:
            s = int(v/2)
            for k in range(m):
                jf = k * v
                jl = jf + s - 1
                jt = 0
                for j in range(jf,jl+1):
                    TW = pow(W,jt,q)

                    temp = B[j]

                    B[j  ] = (temp + B[j+s]) % q
                    B[j+s] = (temp - B[j+s])*TW % q

                    jt = jt + m
            m = 2*m
            v = s

        return B

    # Iterative Radix-2 Decimation-in-Frequency (DIF) (GS) NTT - RN
    # A: input polynomial (reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-standard order)
    def Radix2_DIF_Iterative_NTT_RN(self,A,W,q):
        N = len(A)
        B = [_ for _ in A]

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

                    jt = jt+1
                    j = j + 2*m
            m = 2*m
            v = int(v/2)
            d = 2*d

        return B

    # Iterative Radix-2 Decimation-in-Frequency (DIF) (GS) NTT - NN
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def Radix2_DIF_Iterative_NTT_NN(self,A,W,q):
        N = len(A)
        # requires two memory
        B = [_ for _ in A]
        C = [_ for _ in A]
        # C = [0] * N

        m = 1
        v = N
        d = 1

        if int(log(v//2))%2 == 0:
            nsi = True
        else:
            nsi = False

        while v>1:
            if nsi:
                for jf in range(m):
                    j = jf
                    jt = 0
                    k = jf
                    while j<(N-1):
                        TW = pow(W,jt,q)

                        C[j]   = (B[k] + B[k+int(N/2)]) % q
                        C[j+d] = (B[k] - B[k+int(N/2)])*TW % q

                        jt = jt + m
                        j = j + 2*m
                        k = k + m
                nsi = False
            else:
                for jf in range(m):
                    j = jf
                    jt = 0
                    k = jf
                    while j<(N-1):
                        TW = pow(W,jt,q)

                        B[j]   = (C[k] + C[k+int(N/2)]) % q
                        B[j+d] = (C[k] - C[k+int(N/2)])*TW % q

                        jt = jt + m
                        j = j + 2*m
                        k = k + m
                nsi = True
            m = 2*m
            v = int(v/2)
            d = 2*d

        return C

    ######################################################################## (INTT)
    """
    List of INTT Algorithms: NTT algorithms with extra n^-1 mod q multiplication
    -- Naive NTT (see Wikipedia definition of NTT operation)
    -- Radix-2 Decimation-in-Time (DIT) Recursive NTT (Cooley-Tukey)
    -- Radix-2 Decimation-in-Frequency (DIF) Recursive NTT (Gentleman-Sandle)
    -- Radix-2 Decimation-in-Time (DIT) Iterative NTT
    ---- NR (N: Natural order, R: Reversed Order)
    ---- RN
    ---- NN
    -- Radix-2 Decimation-in-Time (DIF) Iterative NTT
    ---- NR
    ---- RN
    ---- NN
    """

    def NaiveINTT_NN(self,A,W_inv,q):
        """
        Very slow baseline implementation. Do not use for large parameter set.
        """
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.NaiveNTT_NN(A,W_inv,q)]
        return B

    def NaiveINTT_NR(self,A,W_inv,q):
        """
        Very slow baseline implementation. Do not use for large parameter set.
        """
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.NaiveNTT_NR(A,W_inv,q)]
        return B

    def Radix2_DIT_Recursive_INTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIT_Recursive_NTT(A,W_inv,q)]
        return B

    def Radix2_DIF_Recursive_INTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIF_Recursive_NTT(A,W_inv,q)]
        return B

    def Radix2_DIT_Iterative_INTT_NR(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIT_Iterative_NTT_NR(A,W_inv,q)]
        return B

    def Radix2_DIT_Iterative_INTT_RN(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIT_Iterative_NTT_RN(A,W_inv,q)]
        return B

    def Radix2_DIT_Iterative_INTT_NN(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIT_Iterative_NTT_NN(A,W_inv,q)]
        return B

    def Radix2_DIF_Iterative_INTT_NR(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIF_Iterative_NTT_NR(A,W_inv,q)]
        return B

    def Radix2_DIF_Iterative_INTT_RN(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIF_Iterative_NTT_RN(A,W_inv,q)]
        return B

    def Radix2_DIF_Iterative_INTT_NN(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.Radix2_DIF_Iterative_NTT_NN(A,W_inv,q)]
        return B

    """
    CRT-based NTT (it is also used for polynomial multiplication in x^n-1)

    Example for 8-pt NTT (w -> 8th root of unity)
    Input  -> Standard Order
    Output -> Bit-reversed Order (We can make it standard order)

                                      x^8 - 1 --------------------------------------------- Stage #0
                                      /     \
                                     /       \
                                    /         \
                                   /           \
                                  /             \
                                 /               \
                                /                 \
                               /                   \
                              /                     \
                             /                       \
                      x^4 - 1                         x^4 + 1 ----------------------------- Stage #1
                         ||                              ||
                      x^4 - 1                         x^4 - w^4
                     /  \                                   /  \
                    /    \                                 /    \
                   /      \                               /      \
            x^2 - 1        x^2 + 1               x^2 - w^2        x^2 + w^2 --------------- Stage #2
               ||             ||                     ||               ||
            x^2 - 1       x^2 - w^4              x^2 - w^2        x^2 - w^6
           / \               / \                   / \               / \
          /   \             /   \                 /   \             /   \
         /     \           /     \               /     \           /     \
    x - 1     x + 1   x - w^2   x + w^2       x - w   x + w   x - w^3   x + w^3 ----------- Stage #3
      ||        ||       ||       ||            ||      ||       ||       ||
    x - 1    x - w^4  x - w^2   x - w^6      x - w   x - w^5  x - w^3   x - w^7

    -- Recursive
    -- Full
    -- Iterative (converted to an optimized algorithm) --> Already presented above.
    ---- CT
    ---- GS
    """

    # CRT-based NTT (recursive)
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def CRT_Recursive_NTT(self,A,W,q,pw=0,depth=1):
        N = len(A)

        if N == 1:
            return A
        else:
            B = [0]*N

            W_N = pow(W,pw,q)

            # reduce
            B[0:int(N/2)] = [(A[i] + A[i+int(N/2)]*W_N) % q for i in range(int(N/2))]
            B[int(N/2):N] = [(A[i] - A[i+int(N/2)]*W_N) % q for i in range(int(N/2))]

            # recall functions
            B[0:int(N/2)] = self.CRT_Recursive_NTT(B[0:int(N/2)], W,q,int(pw/2)                 ,2*depth)
            B[int(N/2):N] = self.CRT_Recursive_NTT(B[int(N/2):N], W,q,int(pw/2)+int((N/4)*depth),2*depth)

            return B

    # CRT-based NTT (full)
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def CRT_Full_NTT(self,A,W,q):
        N = len(A)
        B = [0]*N

        # If i or j is bit-reversed, output will be in bit-reversed order
        for j in range(N):
            C = [x*pow(W**j,i,q) % q for i,x in enumerate(A)]
            B[j] = sum(C) % q

        return B

    ######################################################################## (INTT)

    """
    CRT-based INTT (it is also used for polynomial multiplication in x^n-1)
    It is NTT algorithms with extra n^-1 mod q multiplication

    -- Recursive
    -- Full
    -- Iterative (converted to an optimized algorithm) --> Already stated algorithms above!
    ---- CT
    ---- GS
    """

    def CRT_Recursive_INTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.CRT_Recursive_NTT(A,W_inv,q)]
        return B

    def CRT_Full_INTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.CRT_Full_NTT(A,W_inv,q)]
        return B

    """
    List of NTT Algorithms: (from literature)
    -- Recursive Cooley-Tukey (CT) NTT (see http://people.scs.carleton.ca/~maheshwa/courses/5703COMP/16Fall/FFT_Report.pdf)
    -- Iterative NTT (see https://eprint.iacr.org/2019/109.pdf)
    -- Constant-Geometry NTT (see https://tches.iacr.org/index.php/TCHES/article/view/8344/7692 or https://eprint.iacr.org/2014/646.pdf)
       (NOTE: There are typos in the Algorithm presented in the papers)
    -- Stockham NTT (see https://ieeexplore.ieee.org/document/8717615)
    -- Four-Step NTT (see https://eprint.iacr.org/2015/818.pdf)
    """

    # Cooley-Tukey NTT
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def CooleyTukeyNTT(self,A,W,q):
        N = len(A)

        if (N == 2):
            B = [0] * N

            B[0] = (A[0] + A[1]) % q
            B[1] = (A[0] - A[1]) % q

            return B
        else:
            B = [0] * N
            w = 1

            A_even = [0] * (N >> 1)
            A_odd  = [0] * (N >> 1)

            for i in range(N >> 1):
                A_even[i] = A[2 * i]
                A_odd[i]  = A[2 * i + 1]

            B_even = self.CooleyTukeyNTT(A_even,(W * W % q),q)
            B_odd  = self.CooleyTukeyNTT(A_odd, (W * W % q),q)

            for i in range(N >> 1):
                B[i]            = (B_even[i] + w * B_odd[i]) % q
                B[i + (N >> 1)] = (B_even[i] - w * B_odd[i]) % q

                w = w * W

        return B

    # Iterative NTT
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (bit-reversed order)
    def IterativeNTT(self,A,W,q):
        N = len(A)
        B = [x for x in A]

        v = int(log(N, 2))

        for i in range(0, v):
            for j in range(0, (2 ** i)):
                for k in range(0, (2 ** (v - i - 1))):
                    s = j * (2 ** (v - i)) + k
                    t = s + (2 ** (v - i - 1))

                    w = (W ** ((2 ** i) * k)) % q

                    as_temq = B[s]
                    at_temq = B[t]

                    B[s] = (as_temq + at_temq) % q
                    B[t] = ((as_temq - at_temq) * w) % q

        return B

    # Four-Step NTT
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # size: input polynomial partition
    # B: output polynomial (standard order)
    def FourStepNTT(self,A,W,q,size):
        """
        This is a unified four-step NTT algorithm for both forward and inverse
        NTT operations. The coefficients of input polynomial should be given in
        standard order. The output is generated in standard order. Forward NTT
        uses twiddle factors and inverse NTT uses modular inverse of twiddle factors.

        This algorithm divides NTT operation into smaller parts. "size" input
        determines the size of these small NTT operations. For details of the
        algorithm, see the paper: https://eprint.iacr.org/2015/818.pdf
        """
        N = len(A)

        poly = [_ for _ in A]

        size0 = size[0]
        size1 = size[1]

        temp0 = 1
        # STEP.1
        if isinstance(size0, list):
            for i in size0:
                temp0 = temp0 * i
            STEP_1 = matrix(poly, N/temp0)
            W_0 = (W ** (N/temp0)) % q
            for i in range(int(N/temp0)):
                STEP_1[i] = self.FourStepNTT(STEP_1[i],W_0,q,size0)
        else:
            temp0 = size0
            STEP_1 = matrix(poly, int(N/temp0))
            W_0 = (W ** int(N/temp0)) % q
            for i in range(int(N/temp0)):
                STEP_1[i] =  self.CooleyTukeyNTT(STEP_1[i],W_0,q)

        # STEP.2 - Transpose
        STEP_2 = [ [row[c] for row in STEP_1 if c < len(row)] for c in range(0, max([len(row) for row in STEP_1])) ]
        # STEP_2 = list(zip(*STEP_1))

        # STEP.3 - Multiply with twiddle factor of N-pt NTT
        STEP_3 = [[0]*int(N/temp0)]*size0
        for i in range(temp0):
            STEP_3[i] = [(STEP_2[i][k] * (W ** (i*k)) % q) for k in range(int(N/temp0))]

        temp1 = 1
        #STEP.4
        if isinstance(size1, list):
            for i in size1:
                temp1 = temp0 * i
            W_1 = (W ** int(N/temp1)) % q
            for i in range(int(N/temp1)):
                STEP_3[i] = self.FourStepNTT(STEP_3[i],W_1,q,size1)
        else:
            temp1 = size1
            W_1 = (W ** int(N/temp1)) % q
            for i in range(int(N/temp1)):
                STEP_3[i] = self.CooleyTukeyNTT(STEP_3[i],W_1,q)

        # Final transpose
        STEP_4 = [ [row[c] for row in STEP_3 if c < len(row)] for c in range(0, max([len(row) for row in STEP_3])) ]
        # STEP_4 = list(zip(*STEP_3))

        # Convert matrix into array
        STEP_4 = [item for sublist in STEP_4 for item in sublist]

        return STEP_4

    # Four-Step NTT v2
    # A: input polynomial (standard order)
    # W: twiddle factor
    # q: modulus
    # size: input polynomial partition
    # B: output polynomial (standard order)
    def FourStepNTTv2(self,A,W,q,size):
        """
        This is a four-step NTT algorithm for both forward and inverse NTT
        operations. The coefficients of tnput polynomial should be given in
        standard order. The output is generated in standard order. Forward NTT
        uses modular inverse of twiddle factors and inverse NTT uses twiddle factors.

        This algorithm divides NTT operation into smaller parts. "size" input
        determines the size of these small NTT operations. For details of the
        algorithm, see the paper: https://eprint.iacr.org/2015/818.pdf
        """
        N = len(A)

        # If this is an inverse transform operation
        N_inv = modinv(N, q)
        # Re-order input
        poly = [A[0]] + list(reversed(A[1:]))

        size0 = size[0]
        size1 = size[1]

        temp0 = 1
        # STEP.1
        if isinstance(size0, list):
            for i in size0:
                temp0 = temp0 * i
            STEP_1 = matrix(poly, N/temp0)
            W_0 = (W ** (N/temp0)) % q
            for i in range(int(N/temp0)):
                STEP_1[i] = self.FourStepNTT(STEP_1[i],W_0,q,size0)
        else:
            temp0 = size0
            STEP_1 = matrix(poly, int(N/temp0))
            W_0 = (W ** int(N/temp0)) % q
            for i in range(int(N/temp0)):
                STEP_1[i] =  self.CooleyTukeyNTT(STEP_1[i],W_0,q)

        # STEP.2 - Transpose
        STEP_2 = [ [row[c] for row in STEP_1 if c < len(row)] for c in range(0, max([len(row) for row in STEP_1])) ]
        # STEP_2 = list(zip(*STEP_1))

        # STEP.3 - Multiply with twiddle factor of N-pt NTT
        STEP_3 = [[0]*int(N/temp0)]*size0
        for i in range(temp0):
            STEP_3[i] = [(STEP_2[i][k] * (W ** (i*k)) % q) for k in range(int(N/temp0))]

        temp1 = 1
        #STEP.4
        if isinstance(size1, list):
            for i in size1:
                temp1 = temp0 * i
            W_1 = (W ** int(N/temp1)) % q
            for i in range(int(N/temp1)):
                STEP_3[i] = self.FourStepNTT(STEP_3[i],W_1,q,size1)
        else:
            temp1 = size1
            W_1 = (W ** int(N/temp1)) % q
            for i in range(int(N/temp1)):
                STEP_3[i] = self.CooleyTukeyNTT(STEP_3[i],W_1,q)

        # Final transpose
        STEP_4 = [ [row[c] for row in STEP_3 if c < len(row)] for c in range(0, max([len(row) for row in STEP_3])) ]
        # STEP_4 = list(zip(*STEP_3))

        # Convert matrix into array
        STEP_4 = [item for sublist in STEP_4 for item in sublist]

        return STEP_4

    # CT-Based Constant-Geometry NTT
    # A: input polynomial (Bit-reversed order)
    # W: twiddle factor
    # q: modulus
    # B: output polynomial (standard order)
    def CTBased_ConstantGeometryNTT(self,A,W,q):
        N = len(A)
        v = int(log(N,2))

        #B = indexReverse(A,v)
        B = [_ for _ in A]
        C = [0 for _ in range(N)]

        for s in range(1,v+1):
            for j in range(int(N/2)):
                k = int(floor(j/(2**(v-s)))*(2**(v-s)))

                TW = pow(W,k,q)

                C[j           ] = (B[2*j] + B[2*j+1]*TW) % q
                C[j + int(N/2)] = (B[2*j] - B[2*j+1]*TW) % q

            if s != v:
                B = [_ for _ in C]

        return C

    ######################################################################## (INTT)
    """
    List of INTT Algorithms: (from literature): NTT algorithms with extra n^-1 mod q multiplication
    -- Recursive Cooley-Tukey (CT) NTT (see http://people.scs.carleton.ca/~maheshwa/courses/5703COMP/16Fall/FFT_Report.pdf)
    -- Iterative NTT (see https://eprint.iacr.org/2019/109.pdf)
    -- Constant-Geometry NTT (see https://tches.iacr.org/index.php/TCHES/article/view/8344/7692 or https://eprint.iacr.org/2014/646.pdf)
    -- Stockham NTT (see https://ieeexplore.ieee.org/document/8717615)
    -- Four-Step NTT (see https://eprint.iacr.org/2015/818.pdf)
    """

    def CooleyTukeyINTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.CooleyTukeyNTT(A,W_inv,q)]
        return B

    def IterativeINTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.IterativeNTT(A,W_inv,q)]
        return B

    def FourStepINTT(self,A,W_inv,q,size):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.FourStepNTT(A,W_inv,q,size)]
        return B

    def FourStepINTTv2(self,A,W_inv,q,size):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.FourStepNTTv2(A,W_inv,q,size)]
        return B

    def CTBased_ConstantGeometryINTT(self,A,W_inv,q):
        N_inv = modinv(len(A),q)
        B = [(x*N_inv) % q for x in self.CTBased_ConstantGeometryNTT(A,W_inv,q)]
        return B
#
