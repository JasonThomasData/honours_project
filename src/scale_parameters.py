def to_shape_scale(k_dispersion_values, R0):

    shape = []
    scale = []
    for k in k_dispersion_values:
        if k == "inf":
            k = 10**7
        
        # is the index of dispersion, k = (\mu)/(\sigma^2)        
        # use the method of moments:
        shape.append(k)
        scale.append(R0/k)

    results = {}
    for i, k, in enumerate(k_dispersion_values):
        results[str(k)] = {
            "shape": shape[i],
            "scale": scale[i]
        }

    return results
