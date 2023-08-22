from sage.all import ZZ, GF, lcm, sqrt, floor
from cornacchia import all_cornacchia
from numpy import argmax
from hnf import lower_hnf_basis
from ideals import small_equivalent_ideal

def Fp_to_int(n):
    """
    Returns element of GF(p) as integer in range (-p/2 ... +p/2]
    """
    if ZZ(n) > n.parent().order() / 2: return ZZ(n) - n.parent().order()
    return ZZ(n)

def find_element_defining_embedding(O, d, t, all_slns=False, filter_func=None):
    """
    Finds an element in quaternion order 'O' with trace 't' and norm 'd'. Set 'all_slns=True' to get all solutions.
        A function can be provided as 'filter_func' which is called when a solution is found to see if it should be counted or not. We use this for filtering primitive solutions.
    """
    slns = [] if all_slns else None
    # Compute the connecting ideal, and find smaller equivalent ideal, to give right order with lower N
    B = O.quaternion_algebra()
    I = B.maximal_order() * O
    J, y = small_equivalent_ideal(I)
    O_new = J.right_order()
    # Put basis in HNF
    basis_hnf = lower_hnf_basis(B, O_new.basis())
    e00, e01, e02, e03 = basis_hnf[0]
    _,   e11, e12, e13 = basis_hnf[1]
    _,   _,   e22, e23 = basis_hnf[2]
    _,   _,   _,   e33 = basis_hnf[3]
    if (e00 == 0) or (e11 == 0) or (e22 == 0) or (e33 == 0):
        return slns
    # Find alpha_0
    alpha_0 = t / (2 * e00)
    if (alpha_0 not in ZZ) or (d not in ZZ):
        return slns
    # Compute a, b, N
    q, p = [ZZ(abs(l)) for l in B.invariants()]
    N = lcm([e.denominator() for e in [e00,e01,e02,e03,e11,e12,e13,e22,e23,e33]])
    N2 = N**2
    # Find residues of alpha_1 mod p
    Fp = GF(p)
    sq_mod_p = Fp(d - (alpha_0 * e00)**2) / Fp(q)
    rt1 = sqrt(sq_mod_p)
    if rt1 not in Fp:
        return slns
    rt2 = -rt1
    residues = [Fp_to_int((rt1 - Fp(alpha_0 * e01)) / Fp(e11)), Fp_to_int((rt2 - Fp(alpha_0 * e01)) / Fp(e11))]
    # compute maximum value of k - for each residue
    temp1 = d - (alpha_0**2)*(e00**2)
    temp1_scaled = N2 * temp1
    temp2 = sqrt(temp1 / q) - alpha_0*e01
    ks = [floor((temp2 - ZZ(r)*e11)/(p*e11)) for r in residues]
    # loop over k decreasing, for each residue
    max_iter = sum([k + 1 for k in ks if k >= 0])
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
        # find all solutions to Cornacchia's
        betas = all_cornacchia(q, v)
        for beta_pair in betas:
            # Check if this gives a solution with integral alpha_2 and alpha_3
            alpha_2 = (beta_pair[0] - N*alpha_1*e12 - N*alpha_0*e02) / (N*e22)
            alpha_3 = (beta_pair[1] - N*alpha_1*e13 - N*alpha_2*e23 - N*alpha_0*e03) / (N*e33)
            if (alpha_2 in ZZ) and (alpha_3 in ZZ):
                alpha = alpha_0*basis_hnf[0] + alpha_1*basis_hnf[1] + alpha_2*basis_hnf[2] + alpha_3*basis_hnf[3]
                alpha_in_O = y * alpha * y**(-1) # map alpha back in to original order
                valid_sln = True
                if filter_func != None:
                    valid_sln = filter_func(alpha_in_O, k)
                if valid_sln:
                    if all_slns: slns.append(alpha_in_O)
                    if not all_slns: return alpha_in_O
    return slns