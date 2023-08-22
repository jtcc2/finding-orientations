from sage.all import factor, floor, QQ, ZZ, RR, PolynomialRing, Integers, sqrt, round

def all_cornacchia(d, m):
    """
    Returns all solutions to x^2 + dy^2 = m, including imprimitive solutions.
    """
    if m < 0: return []
    if m == 0: return [(0, 0)]
    sols = []
    # Iterate over g such that g^2 divides m
    #   Writing m = q1^e1 * q2^e2 *..., we do this by storing an array [f1, f2, ...] where g = q1^f1 * q2^f2 * ..., and 0 <= fi <= floor(ei / 2)
    #   And we increase it in each iteration [0, 0, 0, ...] -> [1, 0, 0, ...] -> [2, 0, 0, ...] -> (then when f1 is maximum) [0, 1, 0, ...] -> [1, 1, 0, ...] -> ...
    fs = factor(m)
    g_fac_arr = [-1] + [0]*len(fs) # store the exponents of g in this array
    while True:
        g_fac_arr[0] += 1
        g = 1
        for i in range(0, len(fs)):
            if g_fac_arr[i] == floor(fs[i][1]/QQ(2)) + 1: 
                g_fac_arr[i] = 0
                g_fac_arr[i+1] += 1
            # expand value of g
            g *= fs[i][0]**g_fac_arr[i]
        if g_fac_arr[-1] != 0: break
        tempm = ZZ(m/(g**2))
        # Find primitive solutions to x^2 + dy^2 = m/g^2 using Cornacchias and scale them to solution (gx, gy)
        P = PolynomialRing(Integers(tempm), 'X')
        X = P.gens()[0]
        rs =[ZZ(r) for r in (X**2 + d).roots(multiplicities=False)]
        bound = round(sqrt(tempm),5)
        for r in rs:
            n = tempm
            while r > bound:
                n, r = r, n%r
            s = sqrt((tempm - r**2)/d)
            if s in ZZ:
                sols.extend([(g*r, g*s), (g*r, -g*s), (-g*r, g*s), (-g*r, -g*s)])
                if d == 1:
                    sols.extend([(g*s, g*r), (g*s, -g*r), (-g*s, g*r), (-g*s, -g*r)])
    return list(set(sols))