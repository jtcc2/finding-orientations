# Code for computing Hermite normal form of basis and matrices
#   Note: The use a column style HNF convention which is different to Sage's built in .hermite_form()

from sage.all import matrix, QQ, ZZ, lcm

def basis_to_matrix(basis):
    """
    Given a basis of an ideal as a list, convert it to a column-style basis matrix.
    """
    return matrix(QQ,[[basis[m][r] for r in range(0, len(basis))] for m in range(0, len(basis))]).transpose()

def matrix_to_basis(B, M):
    """
    Given a quaternion algebra B, and a column-style basis matrix M, return the basis as elements of B.
    """
    return [a[0] for a in (M.transpose() * matrix([B(1)] + list(B.gens())).transpose())]

def lower_hnf_matrix(M):
    """
    Reduces matrix to lower triangular form.
    """
    denom = lcm([l.denominator() for l in M.list()])
    M_ZZ = matrix(ZZ,[[M[k][l]*denom for l in range(0, M.dimensions()[1])] for k in range(0, M.dimensions()[0])])
    return (1/denom) * (M_ZZ.transpose().hermite_form().transpose())

def upper_hnf_matrix(M):
    """
    Reduces matrix to upper triangular form.
    """
    return lower_hnf_matrix(M[::-1,:])[::-1,::-1]

def lower_hnf_basis(B, basis):
    """
    Reduces basis to lower triangular form.
    """
    return matrix_to_basis(B, lower_hnf_matrix(basis_to_matrix(basis)))

def upper_hnf_basis(B, basis):
    """
    Reduces basis to upper triangular form.
    """
    return matrix_to_basis(B, upper_hnf_matrix(basis_to_matrix(basis)))