from sage.all import ZZ, ceil, log, gcd
from random import randint

def small_equivalent_ideal(I):
    """
    Returns a left-ideal J of smaller norm in the left ideal class, and return the element y such that I = yJ
    """
    _,mn = I.quadratic_form().__pari__().qfminim(None,None,1)
    el = sum(ZZ(c)*g for c,g in zip(mn, I.basis()))
    y = el.conjugate() / I.norm()
    I *= y
    return I, y

def random_ideal(O):
    """
    FROM DftP: Sample an integral left O-ideal whose class is approximately uniform.
    """
    Q = O.quaternion_algebra()
    i,j,k = Q.gens()
    p = ZZ(Q.discriminant())
    assert p.is_pseudoprime()

    I = O.unit_ideal()
    for it in reversed(range(55+2*ceil(log(p,2)))):   #TODO figure out good bound
        O = I.right_order()
        while True:
            beta = sum(randint(1,10**4)*a for a in O.basis())
            if gcd(beta.reduced_norm(), 4) == 2:
                break
        J = O*2 + O*beta
        assert J.norm() == 2
        I *= J
        if not it or I.norm() >= p**(ZZ(2)/3):
            # find an ideal of smaller norm in the class
            _,mn = I.quadratic_form().__pari__().qfminim(None,None,1)
            el = sum(ZZ(c)*g for c,g in zip(mn, I.basis()))
            I *= el.conjugate() / I.norm()
    return I