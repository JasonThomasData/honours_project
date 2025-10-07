def to_shape_scale_per_k(k_dispersion_values, R0):
    shape = []
    scale = []
    for k in k_dispersion_values:
        if k == "inf":
            k = 10**7
        shape.append(k)
        scale.append(R0/k)

    results = {}
    for i, k, in enumerate(k_dispersion_values):
        results[str(k)] = {
            "shape": shape[i],
            "scale": scale[i]
        }

    return results

def to_shape_scale_per_alpha(k, m, alpha_values):
    shape = []
    scale = []
    for alpha in alpha_values:
        if k == "inf":
            k = 10**7
        R0 = m*alpha # See thesis for reason
        shape.append(k)
        scale.append(R0/k)

    results = {}
    for i, alpha, in enumerate(alpha_values):
        results[str(alpha)] = {
            "shape": shape[i],
            "scale": scale[i]
        }

    return results

