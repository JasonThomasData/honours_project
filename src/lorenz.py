from scipy.integrate import simpson

def gini_coefficient(line_pts):
    curve_area = simpson(line_pts, dx=1)
    # In all cases, the x-axis is 1:100, and y is 0:1, so the area under the identity function is (1*100)/2 = 50
    B = 50
    A = abs(curve_area - B) # This way, the curve can be above or below, or both, compared to the identity function
    return A/B

def get_Perkins_curve(data, removal_count=100, remove_zero=False):
    # For aggregation and coaggregation data
    if remove_zero:
        data = [x for x in data if x != 0]
    descending_data = sorted(data)
    descending_data.reverse()
    m = len(descending_data)
    v = sum(descending_data)
    
    h_step = m/removal_count
    
    transmission_potential = []
    for i in list(range(0, m)):
        v_i = descending_data[i]/v
        h_i = (i+1)/m
        transmission_potential.append(v_i**2 * h_i)
    total_transmission_potential = sum(transmission_potential)

    # QUESTION - SHOULD WE SORT THE LIST ONCE THE CALCULATION IS DONE? THE INCREASING h_i TERM MEANS THIS LIST WILL BE UNSORTED
    # I SAY NO, IF YOU HAVE A LOOK AT THE PERKINS CHART, THE CURVE THAT INTERPOLATES THE POINTS CHANGES GRADIENTS OFTEN
    #print(transmission_potential)
    #transmission_potential = sorted(transmission_potential)
    #transmission_potential.reverse()
    #print(transmission_potential)
    
    transmission_potential_percentages = []

    h_percent = list(range(0,removal_count+1))

    for j in h_percent:
        if j==0:
            potential_of_removed_hosts = 0
        else:
            h_j = round(j*h_step) # this many vertebrates in calculation
            potential_of_removed_hosts = sum(transmission_potential[0:h_j])/total_transmission_potential
        transmission_potential_percentages.append(potential_of_removed_hosts)

    return transmission_potential_percentages

def get_JR_curve(bnbl, bl, removal_count=100, remove_zero=False):
    # bnbl corresponds to the list of b_n * b_l in thesis
    # bl corresponds to the list of b_l
    if remove_zero:
        bnbl = [x for x in bnbl if x != 0]
    descending_bnbl = sorted(bnbl)
    descending_bnbl.reverse()
    m = len(descending_bnbl)
    bl_sum = sum(bl)

    h_step = m/removal_count

    h_percent = list(range(0,removal_count+1))

    transmission_potential = []
    for i in list(range(0, m)):
        transmission_potential.append(descending_bnbl[i]/bl_sum)
    total_transmission_potential = sum(transmission_potential)

    transmission_potential_percentages = []
    for j in h_percent:
        if j==0:
            potential_of_removed_hosts = 0
        else:
            h_j = round(j*h_step) # this many vertebrates in calculation
            potential_of_removed_hosts = sum(transmission_potential[0:h_j])/total_transmission_potential
        transmission_potential_percentages.append(potential_of_removed_hosts)

    return transmission_potential_percentages

def get_aggregation_curve(data, removal_count=100, remove_zero=False):
    # For aggregation and coaggregation data
    if remove_zero:
        data = [x for x in data if x != 0]
    descending_data = sorted(data)
    descending_data.reverse()
    m = len(descending_data)
    v = sum(descending_data)
    
    h_step = m/removal_count
    
    h_percent = list(range(0,removal_count+1))
    
    aggregation_cummulative = []
    for i in list(range(0, m)):
        v_i = descending_data[i]/v
        aggregation_cummulative.append(v_i)
    aggregation_total = sum(aggregation_cummulative)
    
    aggregation_percentages = []
    for j in h_percent:
        if j==0:
            aggregation_of_removed_hosts = 0
        else:
            h_j = round(j*h_step) # this many vertebrates in calculation
            aggregation_of_removed_hosts = sum(aggregation_cummulative[0:h_j])/aggregation_total
        aggregation_percentages.append(aggregation_of_removed_hosts)

    return aggregation_percentages
