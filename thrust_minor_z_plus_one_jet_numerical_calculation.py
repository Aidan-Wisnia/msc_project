import numpy as np


# General definitions

N = 1000000
epsilon = 1*10**(-6)
R_prime_list = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
decimal_places = 6

print(f"""Number of iterations: {N}.
Cutoff value: {epsilon}.""")
# defining functions used across the contrbutions

def generate_i_emissions(epsilon_cutoff):
    
    first_zeta_emission = 1
    
    # Creating two lists to store the emissions of zeta and phi:
    
    zeta_emissions_for_event = []
    phi_emissions_for_event = []
    
    
    # first random (not for zeta) generation of phi/zeta:
    
    zeta = first_zeta_emission
    phi = 2 * np.pi* np.random.rand()
    
    # storing outputs
    zeta_emissions_for_event.append(zeta)
    phi_emissions_for_event.append(phi)
    
    
    # while loop for all subsequent emissions until zeta > epsilon

    while zeta > epsilon_cutoff:

        # nth random generation of phi/zeta:

        # define the number of emissions so far

        number_of_emissions = len(phi_emissions_for_event)


        zeta = zeta_emissions_for_event[number_of_emissions - 1] * (np.random.rand())**(1/R_prime)
        phi = 2 * np.pi* np.random.rand()

        # storing outputs
        zeta_emissions_for_event.append(zeta)
        phi_emissions_for_event.append(phi)

    return zeta_emissions_for_event, phi_emissions_for_event



def generate_random_special_emission(epsilon_cutoff):
    
    zeta_special_emission = np.random.rand()
    while zeta_special_emission < epsilon_cutoff:
        zeta_special_emission = np.random.rand()
        
    phi_special_emission = 2 * np.pi* np.random.rand()
    return zeta_special_emission, phi_special_emission



def generate_N_events(N):

    zeta_event_list = []
    phi_event_list = []
    zeta_special_event_list = []
    phi_special_event_list = []



    for j in range(N):

        zeta_i, phi_i = generate_i_emissions(epsilon)
        zeta_special, phi_special = generate_random_special_emission(epsilon)


        zeta_event_list.append(zeta_i)
        phi_event_list.append(phi_i)
        zeta_special_event_list.append(zeta_special)
        phi_special_event_list.append(phi_special)
        
    return zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list



def x_calculation(zeta, phi):
    '''
    Calculates x for a given zeta and phi.
    '''
    
    x = zeta * (np.sin(phi) / np.abs(np.sin(phi)))
    return x



def x_i_sum_calculation(zeta_list, phi_list):
    '''
    Calculates x for each emission i that is not special or the first emission.
    '''
    assert zeta_list[0] != 1 , f"First zeta in {zeta_list} is 1"
    assert len(zeta_list) == len(phi_list), f"The length of {zeta_list} and {phi_list} are not the same"
    
    x_i_sum = 0
    
    for i in range(len(zeta_list)):
        x_i_sum = x_i_sum + x_calculation(zeta_list[i], phi_list[i])
        
    return x_i_sum



def separate_first_and_ith_emissions(zeta_list, phi_list):
    '''
    Creates 2 new variables and 2 new lists, the zeta and phi for the first emissions and the 
    zetas and phis for the rest of the non-first emissions.
    '''
    
    zeta_one = zeta_list[0]
    phi_one = phi_list[0]
    zeta_i = zeta_list[1:]
    phi_i = phi_list[1:]
    
    return zeta_one, phi_one, zeta_i, phi_i



def zeta_i_sum_calculation(zeta_list):
    '''
    Calculates the sum of zetas from a list when the first emission is not = 1
    '''
    assert zeta_list[0] != 1 , f"The first emission of {zeta_list} is 1."
    
    zeta_i_sum = 0
    
    for i in range(len(zeta_list)):
        zeta_i_sum = zeta_i_sum + zeta_list[i]
        
    return zeta_i_sum



def generate_z():
    return np.random.rand()



# F_NLL calculation

F_NLL_result_list = []
F_NLL_error_list = []

for R_prime in R_prime_list:

    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)

    F_NLL_result = 0 
    F_NLL_error = 0

    for j in range(N):

        # separate out the zeta_1, phi_1 from the rest of the soft emissions
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        # calculate the x for the soft emissions
        x_sum = x_i_sum_calculation(zeta_i, phi_i)

        # calculate zeta for the soft emissions
        zeta_sum = zeta_i_sum_calculation(zeta_i)

        # calcualte the x contribution due to the first emission
        x_one = x_calculation(zeta_one, phi_one)

        # calculate the weight for F_NLL 
        F_NLL_weight = (1/2 * (zeta_sum + zeta_one + np.abs(x_sum + x_one)))**(-R_prime)

        F_NLL_result = F_NLL_result + ( F_NLL_weight / N )
        F_NLL_error = F_NLL_error + ( (F_NLL_weight**2) / N )

    F_NLL_error = np.sqrt( ( F_NLL_error - F_NLL_result**2 ) / N )

    # collect all total results
    F_NLL_result_list.append(round(F_NLL_result, decimal_places))
    F_NLL_error_list.append(round(F_NLL_error, decimal_places))

print('F_NLL results:')
print("R' values: ",R_prime_list)
print("Result:  ", F_NLL_result_list)
print("Error: ", F_NLL_error_list)



# F_SC_Calculation


F_SC_result_list = []
F_SC_error_list = []

for R_prime in R_prime_list:
    
    # Generate events    
    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)


    SC_result = 0
    SC_error = 0

    for j in range(N):

        # separate out the zeta_1, phi_1 from the rest of the soft emissions
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        # calculate the x for the soft emissions
        x_sum = x_i_sum_calculation(zeta_i, phi_i)

        # calculate zeta for the soft emissions
        zeta_sum = zeta_i_sum_calculation(zeta_i)

        # Get phi and zeta special emissions
        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]

        # SC MINOR calc
        SC_MINOR_zeta_special = zeta_special

        # get x contribution from first emission
        SC_MINOR_x_one = x_calculation(zeta_one, phi_one)

        # get x contribution from special emission
        SC_MINOR_x_special = x_calculation(SC_MINOR_zeta_special, phi_special)


        # define the "real" part of the weight calculation
        SC_MINOR_real = 1/2 * ( (zeta_sum + zeta_one + SC_MINOR_zeta_special) \
                               + np.abs(x_sum + SC_MINOR_x_one + SC_MINOR_x_special) \
                              )


        # define the "virtual" part of the weight calculation
        SC_MINOR_virtual = max(SC_MINOR_zeta_special, \
                                 (1/2 * ( \
                                        zeta_sum + zeta_one + np.abs(x_sum + SC_MINOR_x_one)
                                        ) \
                                 ) \
                                )

        # do weight calculation
        SC_MINOR_weight =  (1/SC_MINOR_zeta_special) * \
        ( \
        (SC_MINOR_real**(-R_prime) * ((1/R_prime) + np.log(SC_MINOR_real) + np.log(1/SC_MINOR_zeta_special))) - \
        (SC_MINOR_virtual**(-R_prime) * ((1/R_prime) + np.log(SC_MINOR_virtual) + np.log(1/SC_MINOR_zeta_special)))\
        )

        # SC MAJOR calc

        # set the zeta special emission to be 1
        SC_MAJOR_zeta_special = 1

        # get x contribution from the redefined special emission
        SC_MAJOR_x_special = x_calculation(SC_MAJOR_zeta_special, phi_special)

        # define "real" part of the weight calc
        SC_MAJOR_real = 1/2 * ( \
                               (SC_MAJOR_zeta_special + zeta_sum) + np.abs(SC_MAJOR_x_special + x_sum) \
                              )

        # define the "virtual" part of the weight 
        SC_MAJOR_virtual = max(1, \
                               1/2 * ( \
                                      zeta_sum + np.abs(x_sum) \
                                     ) \
                              )

        # do weight calc
        SC_MAJOR_weight = (1/R_prime) * \
        ( \
         ( SC_MAJOR_real**(-R_prime) * ((1/R_prime) + np.log(SC_MAJOR_real)) ) - \
         ( SC_MAJOR_virtual**(-R_prime) * ((1/R_prime) + np.log(SC_MAJOR_virtual)) )
        )

        # get results and errors
        SC_result = SC_result + ((SC_MINOR_weight + SC_MAJOR_weight) / N)
        SC_error = SC_error + (((SC_MINOR_weight + SC_MAJOR_weight)**2) / N)
        
    SC_error = np.sqrt((SC_error - SC_result**2)/N)
    
    # collect total results
    F_SC_result_list.append(round(SC_result, decimal_places))
    F_SC_error_list.append(round(SC_error, decimal_places))


print('F_SC Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_SC_result_list)
print("Error: ", F_SC_error_list)


# F_NNLL_HC

# F_SC_Calculation


F_HC_result_list = []
F_HC_error_list = []

for R_prime in R_prime_list:


    # Generate events
    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)


    HC_result = 0
    HC_error = 0

    for j in range(N):

        # separate out the zeta_1, phi_1 from the rest of the soft emissions
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        # calculate the x for the soft emissions
        x_sum = x_i_sum_calculation(zeta_i, phi_i)

        # calculate zeta for the soft emissions
        zeta_sum = zeta_i_sum_calculation(zeta_i)

        # Get phi and zeta special emissions
        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]

        # HC MINOR calc
        HC_MINOR_zeta_special = zeta_special

        # get x contribution from first emission
        HC_MINOR_x_one = x_calculation(zeta_one, phi_one)

        # get x contribution from special emission
        HC_MINOR_x_special = x_calculation(HC_MINOR_zeta_special, phi_special)


        # define the "real" part of the weight calculation
        HC_MINOR_real = 1/2 * (\
                               (zeta_sum + zeta_one + HC_MINOR_zeta_special) + \
                               np.abs(x_sum + HC_MINOR_x_one + HC_MINOR_x_special)
                               )



        # define the "virtual" part of the weight calculation
        HC_MINOR_virtual = max(HC_MINOR_zeta_special, \
                                 (1/2 * ( \
                                        (zeta_sum + zeta_one) + np.abs(x_sum + HC_MINOR_x_one)
                                        ) \
                                 ) \
                                )

        # do weight calculation
        HC_MINOR_weight =  (1/HC_MINOR_zeta_special) * \
        ( (HC_MINOR_real**(-R_prime)) - (HC_MINOR_virtual**(-R_prime)) )


        # HC MAJOR calc

        # set the zeta special emission to be 1
        HC_MAJOR_zeta_special = 1

        # get x contribution from the redefined special emission
        HC_MAJOR_x_special = x_calculation(HC_MAJOR_zeta_special, phi_special)


        # define "real" part of the weight calc
        HC_MAJOR_real = 1/2 * ( \
                               (HC_MAJOR_zeta_special + zeta_sum) + np.abs(HC_MAJOR_x_special + x_sum) \
                              )

        # define the "virtual" part of the weight 
        HC_MAJOR_virtual = max(1, \
                               1/2 * ( \
                                      zeta_sum + np.abs(x_sum) \
                                     ) \
                              )

        # do weight calc
        HC_MAJOR_weight = (1/R_prime) * \
        ( (HC_MAJOR_real**(-R_prime)) - (HC_MAJOR_virtual**(-R_prime)) )

        # get results and errors
        HC_result = HC_result + ((HC_MINOR_weight + HC_MAJOR_weight) / N)
        HC_error = HC_error + (((HC_MINOR_weight + HC_MAJOR_weight)**2) / N)

    HC_error = np.sqrt((HC_error - HC_result**2)/N)

    # collect total results
    F_HC_result_list.append(round(HC_result, decimal_places))
    F_HC_error_list.append(round(HC_error, decimal_places))


print('F_HC Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_HC_result_list)
print("Error: ", F_HC_error_list)



# delta F_rec correction, with C_F pre-factor

F_REC_result_list = []
F_REC_error_list = []

for R_prime in R_prime_list:


    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)

    recoil_result = 0
    recoil_error = 0

    for j in range(N):
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]


        zeta_sum = zeta_i_sum_calculation(zeta_i)

        x_sum = x_i_sum_calculation(zeta_i, phi_i)
        x_one = x_calculation(zeta_one, phi_one)
        x_special = x_calculation(zeta_special, phi_special)

        z = np.random.rand()
        while z < epsilon:
            z = np.random.rand()

        P_q_splitting_function = (z-2) + 2/z
        
        weight_recoil_part_MINOR = 1/2 * ( \
            zeta_sum + zeta_one + \
            np.abs(x_special + (1 - z) * (x_one + x_sum)) + \
            np.abs(x_special - z * (x_one + x_sum)) \
                                         ) ** (-R_prime)

        weight_soft_part_MINOR = 1/2 * ( \
            zeta_sum + zeta_one + zeta_special + \
            np.abs(x_special + (x_one + x_sum)) \
                                         ) ** (-R_prime)

        recoil_MINOR_weight = 1/zeta_special * P_q_splitting_function * (\
                                 weight_recoil_part_MINOR - weight_soft_part_MINOR
                                                )

        zeta_special = 1
        
        x_special = x_calculation(zeta_special, phi_special)

        weight_recoil_part_MAJOR = 1/2 * ( \
            zeta_sum + \
            np.abs(x_special + (1 - z) * (x_sum)) + \
            np.abs(x_special - z * (x_sum)) \
                                         ) ** (-R_prime)

        weight_soft_part_MAJOR = 1/2 * ( \
            zeta_sum + zeta_special + \
            np.abs(x_special + x_sum) \
                                         ) ** (-R_prime)

        recoil_MAJOR_weight = 1/R_prime * P_q_splitting_function * (\
                                 weight_recoil_part_MAJOR - weight_soft_part_MAJOR
                                                )


        recoil_result = recoil_result + ((recoil_MINOR_weight + recoil_MAJOR_weight) / N)

        recoil_error = recoil_error + (((recoil_MINOR_weight + recoil_MAJOR_weight)**2) / N)

    recoil_error = np.sqrt((recoil_error - recoil_result**2)/N)

    
    # collect total results
    F_REC_result_list.append(round(recoil_result, decimal_places))
    F_REC_error_list.append(round(recoil_error, decimal_places))
    

print('delta F_rec, with C_F pre-factor Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_REC_result_list)
print("Errors: ", F_REC_error_list)



# delta F_rec correction, with C_A pre-factor

F_REC_result_list = []
F_REC_error_list = []

for R_prime in R_prime_list:


    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)

    recoil_result = 0
    recoil_error = 0

    for j in range(N):
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]


        zeta_sum = zeta_i_sum_calculation(zeta_i)

        x_sum = x_i_sum_calculation(zeta_i, phi_i)
        x_one = x_calculation(zeta_one, phi_one)
        x_special = x_calculation(zeta_special, phi_special)

        z = np.random.rand()
        while z < epsilon:
            z = np.random.rand()

            
        P_g_C_A_splitting_function = z*(1-z)-2 + 2/z
        
        
        weight_recoil_part_MINOR = 1/2 * ( \
            zeta_sum + zeta_one + \
            np.abs(x_special + (1 - z) * (x_one + x_sum)) + \
            np.abs(x_special - z * (x_one + x_sum)) \
                                         ) ** (-R_prime)

        weight_soft_part_MINOR = 1/2 * ( \
            zeta_sum + zeta_one + zeta_special + \
            np.abs(x_special + (x_one + x_sum)) \
                                         ) ** (-R_prime)

        recoil_MINOR_weight = 1/zeta_special * P_g_C_A_splitting_function * (\
                                 weight_recoil_part_MINOR - weight_soft_part_MINOR
                                                )

        zeta_special = 1
        
        x_special = x_calculation(zeta_special, phi_special)

        weight_recoil_part_MAJOR = 1/2 * ( \
            zeta_sum + \
            np.abs(x_special + (1 - z) * (x_sum)) + \
            np.abs(x_special - z * (x_sum)) \
                                         ) ** (-R_prime)

        weight_soft_part_MAJOR = 1/2 * ( \
            zeta_sum + zeta_special + \
            np.abs(x_special + x_sum) \
                                         ) ** (-R_prime)

        recoil_MAJOR_weight = 1/R_prime * P_g_C_A_splitting_function * (\
                                 weight_recoil_part_MAJOR - weight_soft_part_MAJOR
                                                )


        recoil_result = recoil_result + ((recoil_MINOR_weight + recoil_MAJOR_weight) / N)

        recoil_error = recoil_error + (((recoil_MINOR_weight + recoil_MAJOR_weight)**2) / N)

    recoil_error = np.sqrt((recoil_error - recoil_result**2)/N)

    
    # collect total results
    F_REC_result_list.append(round(recoil_result, decimal_places))
    F_REC_error_list.append(round(recoil_error, decimal_places))
    

print('delta F_rec, with C_A pre-factor Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_REC_result_list)
print("Error: ", F_REC_error_list)


# delta F_rec correction, with n_f T_R

F_REC_result_list = []
F_REC_error_list = []

for R_prime in R_prime_list:


    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)

    recoil_result = 0
    recoil_error = 0

    for j in range(N):
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]


        zeta_sum = zeta_i_sum_calculation(zeta_i)

        x_sum = x_i_sum_calculation(zeta_i, phi_i)
        x_one = x_calculation(zeta_one, phi_one)
        x_special = x_calculation(zeta_special, phi_special)

        z = np.random.rand()
        while z < epsilon:
            z = np.random.rand()

            
        P_g_T_R_n_f_splitting_function = 1-2*z*(1-z)
            
        weight_recoil_part_MINOR = 1/2 * ( \
            zeta_sum + zeta_one + \
            np.abs(x_special + (1 - z) * (x_one + x_sum)) + \
            np.abs(x_special - z * (x_one + x_sum)) \
                                         ) ** (-R_prime)

        weight_soft_part_MINOR = 1/2 * ( \
            zeta_sum + zeta_one + zeta_special + \
            np.abs(x_special + (x_one + x_sum)) \
                                         ) ** (-R_prime)

        recoil_MINOR_weight = 1/zeta_special * P_g_T_R_n_f_splitting_function * (\
                                 weight_recoil_part_MINOR - weight_soft_part_MINOR
                                                )

        zeta_special = 1
        
        x_special = x_calculation(zeta_special, phi_special)

        weight_recoil_part_MAJOR = 1/2 * ( \
            zeta_sum + \
            np.abs(x_special + (1 - z) * (x_sum)) + \
            np.abs(x_special - z * (x_sum)) \
                                         ) ** (-R_prime)

        weight_soft_part_MAJOR = 1/2 * ( \
            zeta_sum + zeta_special + \
            np.abs(x_special + x_sum) \
                                         ) ** (-R_prime)

        recoil_MAJOR_weight = 1/R_prime * P_g_T_R_n_f_splitting_function * (\
                                 weight_recoil_part_MAJOR - weight_soft_part_MAJOR
                                                )


        recoil_result = recoil_result + ((recoil_MINOR_weight + recoil_MAJOR_weight) / N)

        recoil_error = recoil_error + (((recoil_MINOR_weight + recoil_MAJOR_weight)**2) / N)

    recoil_error = np.sqrt((recoil_error - recoil_result**2)/N)

    
    # collect total results
    F_REC_result_list.append(round(recoil_result, decimal_places))
    F_REC_error_list.append(round(recoil_error, decimal_places))
    

print('delta F_rec, with n_f T_R pre-factor Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_REC_result_list)
print("Error: ", F_REC_error_list)


# Capital Delta F_REC with the p_l pre-factor

F_D_REC_result_list = []
F_D_REC_error_list = []

for R_prime in R_prime_list:

    
    # Generate events    
    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)


    D_REC_result = 0
    D_REC_error = 0

    for j in range(N):

        # separate out the zeta_1, phi_1 from the rest of the soft emissions
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        # calculate the x for the soft emissions
        x_sum = x_i_sum_calculation(zeta_i, phi_i)

        # calculate zeta for the soft emissions
        zeta_sum = zeta_i_sum_calculation(zeta_i)

        # Get phi and zeta special emissions
        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]

        # D_REC MINOR calc
        D_REC_MINOR_zeta_special = zeta_special

        # get x contribution from first emission
        D_REC_MINOR_x_one = x_calculation(zeta_one, phi_one)

        # get x contribution from special emission
        D_REC_MINOR_x_special = x_calculation(D_REC_MINOR_zeta_special, phi_special)

        # generate random z
        z = np.random.rand()
        while z < epsilon:
            z = np.random.rand()
        
        p_l_weight = (4 * z * (1 - z) * (2 * (np.cos(phi_special))**2 - 1))
        
        # define the "real" part of the weight calculation
        D_REC_MINOR_real = 1/2 * ( (zeta_sum + zeta_one) \
                               + np.abs(D_REC_MINOR_x_special + (1-z)*(x_sum + D_REC_MINOR_x_one)) \
                               + np.abs(D_REC_MINOR_x_special - z*(x_sum + D_REC_MINOR_x_one)) \
                              )


        # define the "virtual" part of the weight calculation
        D_REC_MINOR_virtual = max(D_REC_MINOR_zeta_special, \
                                 (1/2 * ( \
                                        zeta_sum + zeta_one + np.abs(x_sum + D_REC_MINOR_x_one) \
                                        ) \
                                 ) \
                                )

        # do weight calculation
        D_REC_MINOR_weight =  (1/D_REC_MINOR_zeta_special) * p_l_weight * \
        ( \
        (D_REC_MINOR_real**(-R_prime)) - \
        (D_REC_MINOR_virtual**(-R_prime)) \
        )

        # D_REC MAJOR calc

        # set the zeta special emission to be 1
        D_REC_MAJOR_zeta_special = 1

        # get x contribution from the redefined special emission
        D_REC_MAJOR_x_special = x_calculation(D_REC_MAJOR_zeta_special, phi_special)

        # define "real" part of the weight calc
        D_REC_MAJOR_real = 1/2 * ( \
                               (zeta_sum) + \
                               (np.abs(D_REC_MAJOR_x_special + (1-z)*x_sum)) + \
                               (np.abs(D_REC_MAJOR_x_special - z*x_sum)) \
                              )

        # define the "virtual" part of the weight 
        D_REC_MAJOR_virtual = max(1, \
                               1/2 * ( \
                                      zeta_sum + np.abs(x_sum) \
                                     ) \
                              )

        # do weight calc
        D_REC_MAJOR_weight = (1/R_prime) * p_l_weight * \
        ( \
         ( D_REC_MAJOR_real**(-R_prime)) - \
         ( D_REC_MAJOR_virtual**(-R_prime)) \
        )

        # get results and errors
        D_REC_result = D_REC_result + ((D_REC_MINOR_weight + D_REC_MAJOR_weight) / N)
        D_REC_error = D_REC_error + (((D_REC_MINOR_weight + D_REC_MAJOR_weight)**2) / N)

    D_REC_error = np.sqrt((D_REC_error - D_REC_result**2)/N)

    # collect all total results
    F_D_REC_result_list.append(round(D_REC_result, decimal_places))
    F_D_REC_error_list.append(round(D_REC_error, decimal_places))


print('Capital Delta F_REC with pre-factor Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_D_REC_result_list)
print("Error: ", F_D_REC_error_list)


# delta F correl with the n_f pre-factor

F_CORREL_result_list = []
F_CORREL_error_list = []

for R_prime in R_prime_list:


    zeta_event_list, phi_event_list, zeta_special_event_list, phi_special_event_list = generate_N_events(N)

    corr_result = 0
    corr_error = 0

    for j in range(N):
        zeta_one, phi_one, zeta_i, phi_i = separate_first_and_ith_emissions(zeta_event_list[j], phi_event_list[j])

        zeta_special = zeta_special_event_list[j]
        phi_special = phi_special_event_list[j]

        x_sum = x_i_sum_calculation(zeta_i, phi_i)
        zeta_sum = zeta_i_sum_calculation(zeta_i)

        x_one = x_calculation(zeta_one, phi_one)
        x_special = x_calculation(zeta_special, phi_special)


        # this z is the splitting of the parent not from the leg
        t = np.random.rand()
        z = (np.cos( np.pi/2 * (1-t)))**2
        changed_mu = np.random.rand()


        # This is the new phi introduced for this section:
        relative_phi = np.random.rand() * 2 * np.pi




        nf_pre_factor = \
                    ( \
                     1 - \
                     ( (z * (1 - z)) / (1 + (changed_mu / (1-changed_mu))) ) * \
                     ( \
                      (2 * np.cos(relative_phi)) + \
                      ( ( (1 - 2 * z) / (np.sqrt(z * (1 - z))) ) * (np.sqrt( changed_mu / (1 - changed_mu) )) ) \
                     )**2 \
                    )

        f_correl = np.sqrt( ( changed_mu / (1-changed_mu) ) * ( z*(1-z) ) ) * \
                   zeta_special/2 * (np.abs(np.sin(relative_phi + phi_special)) / np.abs(np.sin(phi_special)))

        corr_MINOR_weight = \
        (nf_pre_factor/(zeta_special*changed_mu*2)) * (np.pi * np.sqrt( z*(1-z) )) * ( \
                                         ( \
                                         1/2 * ( \
                                                (zeta_one + zeta_sum) + \
                                                np.abs((1-z)*x_special + f_correl) + \
                                                np.abs(z*x_special - f_correl) + \
                                                np.abs(x_one + x_sum + x_special) \
                                               ) \
                                         )**(-R_prime) - \
                                         ( \
                                         1/2 * ( \
                                                (zeta_one + zeta_sum + zeta_special) + \
                                                np.abs(x_one + x_sum + x_special) \
                                               ) \
                                         )**(-R_prime) \
                                        )

        zeta_special = 1
        x_special = x_calculation(zeta_special, phi_special)


        f_correl = np.sqrt( ( changed_mu / (1-changed_mu) ) * ( z*(1-z) ) ) * \
                   zeta_special/2 * (np.abs(np.sin(relative_phi + phi_special)) / np.abs(np.sin(phi_special)))



        corr_MAJOR_weight = \
        (nf_pre_factor/(R_prime*changed_mu*2)) * (np.pi * np.sqrt( z*(1-z) )) * ( \
                                    ( \
                                     1/2 * ( \
                                            (zeta_sum) + \
                                            np.abs((1-z)*x_special + f_correl) + \
                                            np.abs(z*x_special - f_correl) + \
                                            np.abs(x_sum + x_special) \
                                           ) \
                                    )**(-R_prime) - \
                                    ( \
                                     ( \
                                     1/2 * ( \
                                            (zeta_sum + zeta_special) + \
                                            np.abs(x_sum + x_special) \
                                           ) \
                                     ) \
                                    )**(-R_prime) \
                                   )

        corr_result = corr_result + ((corr_MINOR_weight + corr_MAJOR_weight) / N)

        corr_error = corr_error + (((corr_MINOR_weight + corr_MAJOR_weight)**2) / N)

    corr_error = np.sqrt((corr_error - corr_result**2) / N)

    # collect all total results
    F_CORREL_result_list.append(round(corr_result, decimal_places))
    F_CORREL_error_list.append(round(corr_error, decimal_places))


print('delta F correl with the n_f pre-factor Results:')
print("R' values: ", R_prime_list)
print("Result: ", F_CORREL_result_list)
print("Error: ", F_CORREL_error_list)

