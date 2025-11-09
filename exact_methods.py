from imports import *

def Liouvillian(H, Ls, hbar = 1):
    d = len(H) # dimension of the system
    superH = -1j/hbar * ( np.kron(np.eye(d),H)-np.kron(H.T,np.eye(d)) ) # Hamiltonian part
    superL = sum([np.kron(L.conjugate(),L) 
                  - 1/2 * ( np.kron(np.eye(d),L.conjugate().T.dot(L)) +
                            np.kron(L.T.dot(L.conjugate()),np.eye(d)) 
                          ) for L in Ls])
    return superH + superL

def Enlarge_Matrix_site_j(j, N, matrix):
    # I⊗...⊗I⊗M⊗I...⊗I: Convert local operators into global operators.
    # j: site, starts in 0.
    
    M = np.eye(len(matrix))
    if j == 0: M = matrix
    
    for i in range(1,N):
        if i == j: M = np.kron(M, matrix)
        else: M = np.kron(M, np.eye(len(matrix)))        

    return M

def Correlation_Matrix_i_Matrix_j(i,j,N,matrix_i, matrix_j):
    # I⊗...⊗I⊗M⊗I...⊗I⊗M⊗I⊗I...⊗I

    M = np.eye(len(matrix_i))
    
    if j == 0: M = matrix_j
    elif i == 0: M = matrix_i
    
    for k in range(1,N):
        if k == j: M = np.kron(M, matrix_j)
        elif k == i: M = np.kron(M, matrix_i)
        else: M = np.kron(M, np.eye(len(matrix_i)))        

    return M

def S_Term(N, cte_list, SigmaMatrix):
    # I⊗...⊗I⊗ΔSigma⊗Sigma⊗I...⊗I: Sigma can be Sigmax, Sigmay or Sigmaz.
    # cte_list = [cte_L, cte_M, cte_R], can be J or ∆.

    Matrix_Sigma = np.zeros((2**N, 2**N))

    cte_L = cte_list[0]
    cte_M = cte_list[1]
    cte_R = cte_list[2]    
    
    for i in range(0,N-1):
        M = np.eye(len(SigmaMatrix))
        
        if i == 0: M = SigmaMatrix
       
        for j in range(1,N):
            if j == i or j == i + 1: M = np.kron(M, SigmaMatrix)
            else: M = np.kron(M, np.eye(len(SigmaMatrix)))        

        if i < N/2 - 1: cte = cte_L
        elif i > N/2 - 1: cte = cte_R
        else: cte = cte_M

        Matrix_Sigma = Matrix_Sigma + M*cte #cte can be ∆_i or J_i

    return Matrix_Sigma #∑ I⊗...⊗I⊗ΔSigma⊗Sigma⊗I...⊗I

def build_number_op_list(N):
    """
    Builds a list of number operators, where each element corresponds
    to the number operator acting on a specific site.
    """
    return [Enlarge_Matrix_site_j(j, N, numberop) for j in range(N)]






# def perform_exact_diag(gamma_L, F_L, gamma_R, F_R, dt, nt, initial_state,H,N):

#     print("This is N: ", N)
    
#     L_K = []
    
#     L_K.append(np.sqrt(gamma_L*(1-F_L))*Enlarge_Matrix_site_j(0, N, Sigma_minus))  
#     L_K.append(np.sqrt(gamma_L*F_L)*Enlarge_Matrix_site_j(0, N, Sigma_plus)) 
#     L_K.append(np.sqrt(gamma_R *(1-F_R))*Enlarge_Matrix_site_j(N-1, N, Sigma_minus))
#     L_K.append(np.sqrt(gamma_R * F_R)*Enlarge_Matrix_site_j(N-1, N, Sigma_plus))
     
#     Superoperator = Liouvillian(H, L_K)



#     # Create time evolution operator
#     d = len(H)
#     U = scipy.linalg.expm(Superoperator * dt)
#     rho_t = initial_state.reshape(d**2,1)  # Vectorized  state
#     number_ops = build_number_op_list(N)

#     expectation_value_history= [[] for qubit in range(N)]

#     rho_matrix = initial_state / np.trace(initial_state)
#     for site in range(N):
#         expectation_value_history[site].append(np.trace(number_ops[site] @ rho_matrix)/np.trace(initial_state))
#         print("The number op I am using!", number_ops[site])



#     print("Initial expectation value of number operator:", expectation_value_history[0])
#     # Time evolution loop
#     for step in range(1,nt+1):
        
#         rho_t = U @ rho_t
#         rho_matrix = rho_t.reshape(d ,d)
#         rho_matrix = rho_matrix / np.trace(rho_matrix)

#         for site in range(N):
#             expectation_value_history[site].append(np.trace(number_ops[site] @ rho_matrix))

#     return expectation_value_history



def perform_exact_diag(gamma_L, F_L,gamma_R, F_R, dt, nt, initial_state, H, N):

        # Build Lindblad operators
        L_K = [    np.sqrt(gamma_L*F_L)  * Enlarge_Matrix_site_j(0, N, Sigma_minus),
            np.sqrt(gamma_L*(1-F_L)) * Enlarge_Matrix_site_j(0, N, Sigma_plus),
            np.sqrt(gamma_R*(1-F_R)) * Enlarge_Matrix_site_j(N-1, N, Sigma_plus),
            np.sqrt(gamma_R*F_R) * Enlarge_Matrix_site_j(N-1, N, Sigma_minus)]

        # Construct superoperator
        Superoperator = Liouvillian(H, L_K)
        d = len(H)


        # Time evolution operator
        U = scipy.linalg.expm(Superoperator * dt)
        
        # Initialize state (column-major vectorization)
        rho_t = initial_state.reshape(d**2, 1, order='F')
        number_ops = build_number_op_list(N)
        expectation_value_history = [[] for _ in range(N)]
     

        # Initial measurements
        rho_matrix = initial_state.copy()
       
        for site in range(N):
            exp_val = np.real(np.trace(number_ops[site] @ rho_matrix))
            expectation_value_history[site].append(exp_val)
            print("Exact diag initial:" , exp_val)
           

        # # Time evolution loop
        for step in range(1, nt+1):
                print("Exact step ", step, " out of", nt)
                rho_t = U @ rho_t
                rho_matrix = rho_t.reshape(d, d, order='F')
                # Store results
   
                rho_matrix = rho_matrix / np.trace(rho_matrix)

                for site in range(N):
                    exp_val = np.real(np.trace(number_ops[site] @ rho_matrix))
                    expectation_value_history[site].append(exp_val)

        # for step in range(1, nt+1):
        #         rho_t = U @ rho_t
        #         rho_matrix = rho_t.reshape(d, d, order='F')
        #         # Store results
   
          

        #         for site in range(N):
        #             exp_val = np.real(np.trace(number_ops[site] @ rho_matrix))
        #             expectation_value_history[site].append(exp_val)

        return expectation_value_history, Superoperator


def build_exact_diag_hamiltonian(J, epsilon):
    N = len(epsilon)

    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    
    # On-site energy terms (ε a_j^† a_j)
    for j in range(N):
        # a_j^† a_j = (σ_j^- σ_j^+) = (1 - σ_j^z)/2
        Z_j = Enlarge_Matrix_site_j(j, N, Sigma_z)
        H += epsilon[j] * 0.5 * (np.eye(dim) - Z_j)
    
    # Hopping terms 
    for j in range(N-1):
        H += J*Correlation_Matrix_i_Matrix_j(j,j+1,N, Sigma_x, Sigma_x)
        H += J*Correlation_Matrix_i_Matrix_j(j,j+1,N, Sigma_y, Sigma_y) 
    return H


