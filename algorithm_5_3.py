from sage.all import matrix, gcd
from hnf import upper_hnf_matrix, basis_to_matrix
from algorithm_5_1 import find_element_defining_embedding

def check_primitive(alpha, O, to_f_matrix=None):
    """
    Given solution alpha in an order O, returns True if solution is primitive, False otherwise.

    Works by computing the matrix mapping alpha's coefficients in 1,i,j,k to coefficients of upper triangular HNF matrix.
        For efficienty, to avoid repeating this computation you can provide it as 'to_f_matrix'.
    """
    # compute coeffs of alpha in upper triangular HNF basis
    if to_f_matrix == None:
        to_f_matrix = upper_hnf_matrix(basis_to_matrix(O.basis())).inverse().transpose()
    f_coeffs = (matrix(alpha.coefficient_tuple()) * to_f_matrix).list()
    # remove duplicate values
    S = [f_coeffs[1]]
    if f_coeffs[2] not in S: S.append(f_coeffs[2])
    if f_coeffs[3] not in S: S.append(f_coeffs[3])
    S = [m for m in S if m != 0]
    if len(S) == 0: return False
    # check gcd
    return gcd(S) == 1

def find_element_defining_primitive_embedding(O, d, t, all_slns=False):
    """
    Finds element of order 'O' of trace t and norm d that is a primitive solution.
    """
    to_f_matrix = basis_to_matrix(O.basis()).inverse().transpose()

    def check_sln(alpha, k):
        return check_primitive(alpha, O, to_f_matrix)
    
    # We pass the above function into the search as a validity check, so we know the solution we get must be primitive 
    return find_element_defining_embedding(O, d, t, all_slns, filter_func=check_sln)