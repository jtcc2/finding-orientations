from sage.all import matrix, diagonal_matrix, vector, is_pseudoprime, QQ, ZZ, GF, lcm, sqrt, log, ceil, floor, denominator, round, random_prime, next_prime, QuaternionAlgebra, kronecker, ceil, DiagonalQuadraticForm
import itertools
from numpy import argmax
from cornacchia import all_cornacchia
from hnf import lower_hnf_basis
from algorithm_5_1 import Fp_to_int
from ideals import random_ideal
from random import randint

def factors_easily(n, B=2**20):
    """
    Given a number n, checks if a n is "Cornacchia Friendly" (= easily factorable)
    """
    n = ZZ(n)
    if n < 0: return False
    if n < 2**160: return True
    l,_ = n.factor(limit=B)[-1]
    return l < 2**160 or is_pseudoprime(l)

def quat_algs(p):
    """
    Generate 3 isomorphic quaternion algebras ramified at p \neq 2 and infity, with abs(i^2) small
    """
    Bs = []
    mod = 4
    if p % 4 == 3:
        Bs.append(QuaternionAlgebra(-1, -p))
    q = 1
    while len(Bs) < 3:
        q = next_prime(q)
        if (-q) % mod == 1 and kronecker(-q, p) == -1:
            Bs.append(QuaternionAlgebra(-q, -p))
    assert all([B.ramified_primes() == [p] for B in Bs ])
    return Bs

def isomorphism_gamma(B_old, B):
    """
    Defines an isomorphism of quaternion algebras, See Lemma 10 [EPSV23]
    """
    if B_old == B:
        return B(1), B(1)
    i_old, j_old, k_old = B_old.gens()
    q_old = -ZZ(i_old**2)
    i, j, k = B.gens()
    q = -ZZ(i**2) 
    p = -ZZ(j**2)
    x, y = DiagonalQuadraticForm(QQ, [1,p]).solve(q_old/q)
    return x + j*y, (x + j_old*y)**(-1)

def eval_isomorphism(alpha, B, gamma):
    """
    Evaluates a quaternion in an isomorphism of quaternion algebras
    """
    i, j, k = B.gens()
    return sum([coeff*b for coeff, b in zip(alpha.coefficient_tuple(), [1, i*gamma, j, k*gamma])])

def find_element_defining_embedding_randomized(O, d, t, filter_func=None):
    """
        Continuously randomizes basis for the order O, until an element of trace t norm d can be found without doing any hard factorizations.
            Maps the element back into the starting order and returns it.
        
        Returns two values:
        - The element of trace t norm d.
        - A "confidence" boolean.
            If no element found, True if we're certain it's not possible an element exists. False if it still might be possible as we may have skipped it.
    """
    def find_element_defining_embedding_with_skips(O, d, t):
        """
        Attempts to find an element in quaternion order 'O' with trace 't' and norm 'd', but may miss solutions.
            In solving x^2+|a|y^2=v with Coracchias, skips if v if it is hard to factor.
            Returns solution, and 'confidence' boolean that is True if no solutions have been skipped.
        """
        # Put basis in HNF
        basis_hnf = lower_hnf_basis(B, O.basis())
        e00, e01, e02, e03 = basis_hnf[0]
        _,   e11, e12, e13 = basis_hnf[1]
        _,   _,   e22, e23 = basis_hnf[2]
        _,   _,   _,   e33 = basis_hnf[3]
        if (e00 == 0) or (e11 == 0) or (e22 == 0) or (e33 == 0):
            return None, True
        # Find alpha_0
        alpha_0 = t / (2 * e00)
        if (alpha_0 not in ZZ) or (d not in ZZ):
            # If none works, we're confident no solution exists in any conjugate order
            return None, True
        # Compute N
        N = lcm([e.denominator() for e in [e00,e01,e02,e03,e11,e12,e13,e22,e23,e33]])
        N2 = N**2
        # Find residues of alpha_1 mod p
        Fp = GF(p)
        sq_mod_p = Fp(d - (alpha_0 * e00)**2) / Fp(q)
        rt1 = sqrt(sq_mod_p)
        if rt1 not in Fp:
            return None, True
        rt2 = -rt1
        residues = [Fp_to_int((rt1 - Fp(alpha_0 * e01)) / Fp(e11)), Fp_to_int((rt2 - Fp(alpha_0 * e01)) / Fp(e11))]
        # compute maximum value of k - for each residue
        temp1 = d - (alpha_0**2)*(e00**2)
        temp1_scaled = N2 * temp1
        temp2 = sqrt(temp1 / q) - alpha_0*e01
        ks = [floor((temp2 - ZZ(r)*e11)/(p*e11)) for r in residues]
        # loop over k decreasing, for each residue
        max_iter = sum([k + 1 for k in ks if k >= 0])
        skipped_v = False
        while max(ks) >= 0:
            k_index = argmax(ks)
            k = ks[k_index]
            r = residues[k_index]
            ks[k_index] = ks[k_index] - 1
            # Compute u and v (v = RHS for Cornacchia)
            alpha_1 = ZZ(r) + k*p
            gamma_1 = alpha_0*e01 + alpha_1*e11
            u = q * N2 * gamma_1**2
            v = ZZ((temp1_scaled - u) / p)
            if factors_easily(v):
                # find all solutions to Cornacchia's
                betas = all_cornacchia(q, v)
                for beta_pair in betas:
                    # Check if this gives a solution with integral alpha_2 and alpha_3
                    alpha_2 = (beta_pair[0] - N*alpha_1*e12 - N*alpha_0*e02) / (N*e22)
                    alpha_3 = (beta_pair[1] - N*alpha_1*e13 - N*alpha_2*e23 - N*alpha_0*e03) / (N*e33)
                    if (alpha_2 in ZZ) and (alpha_3 in ZZ):
                        alpha = alpha_0*basis_hnf[0] + alpha_1*basis_hnf[1] + alpha_2*basis_hnf[2] + alpha_3*basis_hnf[3]
                        valid_sln = True
                        if filter_func != None:
                            valid_sln = filter_func(alpha, k)
                        if valid_sln:
                            return alpha, (not skipped_v)
            else:
                skipped_v = True
        # If we didn't skip any v's we know no solution exists
        return None, (not skipped_v)

    p = -ZZ(O.quaternion_algebra().gens()[1]**2)
    Bs = quat_algs(p)
    for B in Bs:
        # The maximal order with small denominator O_0
        O0 = B.maximal_order()
        # Compute isomorphism between the quat algebras
        gamma, gamma_inv = isomorphism_gamma(O.quaternion_algebra(), B)
        # Transfer the maximal order to new quaternion algebra
        O_in_new_quatalg = B.quaternion_order([eval_isomorphism(alpha, B, gamma) for alpha in O.gens()])
        print(f"\nFinding solution in {O_in_new_quatalg}")
        q, p = [ZZ(abs(l)) for l in B.invariants()]
        # Find connecting ideal
        I = O0 * O_in_new_quatalg
        I = I * denominator(I.norm())
        # Reduced basis to find other small equivalent ideals, which gives suitable isomorphisms of O
        basis_hnf = lower_hnf_basis(B, I.basis())
        M = matrix(QQ, [ai.coefficient_tuple() for ai in basis_hnf])
        S = 2**ceil(log(p*q, 2))
        D = diagonal_matrix(round(S * sqrt(g.reduced_norm())) for g in B.basis())
        reduced_basis = (M * D).LLL() * ~D
        # Define constants for conjugating order
        used = []
        max_size = round(p**(1/1.8),5) + 10
        bound = max(round(log(p,2)/10), 10)
        # Try a bunch of small connecting ideals
        for (a1,a2,a3,a4) in itertools.product(range(0,bound+1), range(-bound,bound+1), range(-bound,bound+1), range(-bound,bound+1)):
            coeffvec = vector(QQ, [a1,a2,a3,a4])
            y = coeffvec * reduced_basis * vector(B.basis())
            Jnorm = y.reduced_norm() / I.norm()
            if y in used or Jnorm > max_size:
                continue
            used.append(y)
            y = y.conjugate() / I.norm()
            J = I * y
            beta, confidence = find_element_defining_embedding_with_skips(J.right_order(), d, t)
            if beta:
                beta_new =  y * beta * y**(-1)
                return eval_isomorphism(beta_new, O.quaternion_algebra(), gamma_inv)
            else:
                if confidence: return None, True
    return None, False