from sage.all import factor

def get_isogeny(E, q, i):
    """
    Returns the 'i'th isogeny of degree 'q' from elliptic curve 'E'.
    Note that the base field of 'E' must be large enough for a basis of the q-torsion to exist.
    """
    P, Q = E.torsion_basis(q) # Gives error if field extension too small
    return E.isogeny(P + i*Q, algorithm="factored")

def get_automorphisms_of_order(E, ord):
    """
    Returns automorphism of order 'ord' on curve E. Returns None if such an automorphism doesn't exist.
    """
    for iso in E.isomorphisms(E):
        sf = iso.scaling_factor()
        if (sf != 1 and sf != -1 and sf**ord == -1):
            return iso
    return None

def is_endo_trace_zero(endomorphism):
    """
    Returns true if endomorphism has trace zero, false otherwise.
    """
    return endomorphism == -endomorphism.dual()

def find_endomorphisms_from_degrees(start_curve, remaining_qs, on_found_endo=None, isog_chain=[]):
    """
    Finds endomorphisms by walking in isogeny graph with a sequence of given degrees.

    Given starting curve 'start_curve',
        a list of isogenies representing a walk from the starting curve 'isog_chain' to the current curve,
        and a list of primes q, 'remaining_qs', where it remains to evaluate all possible q-isogenies,
        recursively performs a depth-first search of the isogeny graph to find endomorphisms.
    Input 'on_found_endo' provides a function that is called whenever an endomorphism is found. If it returns True, we stop.

    Example:
        evaluate_path(E, [3,3,7])
        Explore all paths of 3-isogenies from E, then from those codomains all 3-isogenies, and from those codomains all 7-isogenies, recording all endomorphisms.
    """
    # Take the next prime from the 'remaining_qs' list
    if remaining_qs == []:
        return
    q = remaining_qs[0]
    remaining_qs = remaining_qs[1:]

    # We continue walking from the codomain of the last isogeny
    cur_curve = start_curve if len(isog_chain) == 0 else isog_chain[-1].codomain()

    for i in range(0, q+1):
        # Walk to next curve
        phi = get_isogeny(cur_curve, q, i)
        new_curve = phi.codomain()
        if (new_curve.is_isomorphic(start_curve) == True):
            # We have found an endomorphism, add it to the list
            iso = new_curve.isomorphism_to(start_curve)
            endo = phi
            for j in range(0, len(isog_chain)):
                endo = endo.pre_compose(isog_chain[-j-1])
            endo = endo.post_compose(iso)
            # Call function to handle finding an endomorphism
            if on_found_endo(endo):
                return True
        # Perform next step in walk
        new_isog_chain = isog_chain.copy()
        new_isog_chain.append(phi)
        if find_endomorphisms_from_degrees(start_curve, remaining_qs, on_found_endo, new_isog_chain):
            return True

def solve_orienting_problem(E, d, find_all=True, find_suborder_orientations=True, do_print=True):
    """
    Solves the O-Orienting Problem for a curve 'E' and imaginary quadratic order in the form Z[w]=Z[sqrt{-d}], i.e. generator w has trace 0 and norm 'd'.
    Note:
    - E must be given over a large enough field extension so a basis for the q-torsion exists for each q | d.
    - If d has too many factors the algorithm will not terminate.
    """
    # Get special automorphisms - these are automorphisms which can change the trace of the resulting endomorphism
    auts = []
    if E.j_invariant() == 1728:
        auts.append(get_automorphisms_of_order(E, 2))
    if E.j_invariant() == 0:
        aut = get_automorphisms_of_order(E, 3)
        auts.append(aut, aut**2)

    # Construct the array of isogeny degrees we want in the isogeny path. Factors of d could appear multiple times depending on exponent.
    qs = factor(d)
    remaining_qs = []
    for q in qs:
        remaining_qs.extend([q[0]] * q[1])

    # Find all isogeny paths with degrees given by 'remaining_qs' array (or partial paths) which give endomorphisms
    #    The following function is called when an endomorphism is found
    valid_endomoprhisms = []
    endomorphisms_orienting_suborders = []
    def on_found_endo(endo):
        # Check trace is zero, or if is there an automorphism that makes it trace zero
        trace_zero = is_endo_trace_zero(endo)
        if not trace_zero:
            aut_trace_zero = [is_endo_trace_zero(endo.post_compose(aut)) for aut in auts]
            if True in aut_trace_zero:
                endo = endo.post_compose(auts[aut_trace_zero.index(True)])
                trace_zero = True
        if not trace_zero:
            if do_print: print("Found endomorphism of degree: " + str(endo.degree()) + " - not trace zero")
            return False
        # Check the degree is what we want
        if endo.degree() == d:
            if do_print: print("Found trace zero endomorphism of correct degree!")
            valid_endomoprhisms.append(endo)
            if not find_all: return True
        else:
            if do_print: print("Found trace zero endomorphism found of smaller degree: " + str(endo.degree()) + " so a suborder might orient the curve.")
            endomorphisms_orienting_suborders.append(endo)

    find_endomorphisms_from_degrees(E, remaining_qs, on_found_endo)

    if do_print:
        print()
        if len(valid_endomoprhisms) > 0:
            print("Success, found " + str(len(valid_endomoprhisms)) + " trace zero endomorphisms of correct degree.")
        if len(valid_endomoprhisms) == 0 and len(endomorphisms_orienting_suborders) > 0:
            print("Failed to find trace zero endomorphisms of correct degree. But, at least" + str(len(endomorphisms_orienting_suborders)) + " endomorphisms of smaller degree exist, each orienting a suborder.")
        if len(valid_endomoprhisms) == 0 and len(endomorphisms_orienting_suborders) == 0:
            print("Failed to find any orientations.")

    return valid_endomoprhisms, endomorphisms_orienting_suborders