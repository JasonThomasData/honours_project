def simulate_outbreak(offspring_distribution_generator, target):
    
    outbreak = {
        "history": [[1]], # Generation 0 
        "extinct": False,
        "population": 1
    }
    
    target_reached = False
    i=0
    
    while not target_reached:
        next_generation = []
        for j in range(0,len(outbreak["history"][i])):
            if outbreak["history"][i][j] == 0:
                continue
            for _ in range(0, outbreak["history"][i][j]):
                offspring_for_this_individual = offspring_distribution_generator.get()
                next_generation.append(offspring_for_this_individual)

        target_reached = (sum(next_generation) >= target)

        outbreak["history"].append(next_generation)
        outbreak["population"] += sum(next_generation)
        if sum(next_generation) == 0:
            outbreak["extinct"] = True
            return outbreak
        
        i += 1

    return outbreak
    