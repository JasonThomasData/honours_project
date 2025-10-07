def g(s, R0, k):
    return (1 + R0/k *(1 - s))**(-k)

def fixed_point_iteration(g, R0_i, k):
    # Burden, Annette & Burden, Richard & Faires, J.. (2011). Numerical Analysis, 9th ed.

    p0 = 0.5 # this is a safe bet, since probabilities in [0,1]
    tolerance = 0.0001
    N = 1000 # max iterations

    for i in range(1,N):
        p = g(p0, R0_i, k)
        if abs(p - p0) < tolerance:
            return p
            break
        p0 = p

def generate_for_k(g, R0, k):
    p_series = []
    for _, R0_i in enumerate(R0):
        p = fixed_point_iteration(g, R0_i, k)
        p_series.append(p)
    return p_series
