def build_initial_states(ham_real, layers):


    """
    Builds Initial States for Exact Diagonalization and VQTE simulations.
    
    This function creates:
    1. A parameterized quantum circuit (ansatz) for VQTE
    2. Initial parameter values for the ansatz
    3. The corresponding quantum statevector
    4. A density matrix representation for exact diagonalization
    
    Inputs:
        ham_real : The real Hamiltonian for the VQTE simulation (determines number of qubits)
    
    Returns:
        init_state : Statevector - Initial quantum state for VQTE
        initial_state : np.matrix - Density matrix representation for exact diagonalization
        ansatz : QuantumCircuit - Parameterized circuit used for VQTE
        init_param_values : dict - Dictionary of initial parameter values for the ansatz
    """
    # Create an ansatz circut with reps
    ansatz = EfficientSU2(ham_real.num_qubits, reps = layers)

    #Initialize param dictionary
    init_param_values = {}

    # Set all params to 2π initially
    for i in range(len(ansatz.parameters)):
        init_param_values[ansatz.parameters[i]] = (
        2*np.pi)
      
    # Assign params to the ansatz
    vqte_init_state = Statevector(ansatz.assign_parameters(init_param_values))
    

    # Copy initial state data to a vector
    psi_vector = vqte_init_state.data

    # Reshape to a matrix
    rho_matrix = psi_vector.reshape(2 ,2, order='F')
    exact_diag_initial_state = np.matrix(rho_matrix)

    return vqte_init_state, exact_diag_initial_state, init_param_values, ansatz

def output_results(vqte_results, exact_diag_results, time, nt,time_points):

    # Build graph
    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt+1)
 
    # Plot results
    plt.plot(time_points, exact_diag_results, label='Exact Diag Result', marker='', linestyle='solid')
    plt.plot(time_axis, vqte_results,marker='', linestyle='dashed', label='VQTE Result', color='blue')

    # Titles and finalize plot
    plt.title("Comparison of VQTE and Exact Time Evolution for a Qubit coupled to two thermal reserviors")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()


def verify_density_matrix(rho):

    """
    Verifies that a Density Matrix is Hermatian

    Inputs:
        rho: (density matrix) The density matrix to be checked.
    
    """
        
    # Check Hermitian
    hermitian = np.allclose(rho, rho.conj().T)
    print(f"Is Hermitian: {hermitian}")
    
    # Check trace is 1
    trace = np.trace(rho)
    print(f"Trace: {trace} (should be 1)")
    
    # Check positive semidefinite
    eigenvalues = np.linalg.eigvalsh(rho)
    print(f"Eigenvalues: {eigenvalues}")
    print(f"All eigenvalues ≥ 0: {np.all(eigenvalues >= 0)}")
    
    # Check purity 
    purity = np.trace(rho @ rho)
    print(f"Purity (Tr(ρ²)): {purity} (should be 1 for pure state)")


def calculate_fidelity(vqte_results, exact_results):
    """
    """
    fidelities = []
    
    for i in range(min(len(vqte_results), len(exact_results))):


        rho_vqte = vqte_results[i]
        rho_exact = exact_results[i]
            
        sqrt_rho_vqte = scipy.linalg.sqrtm(rho_vqte)
        product = sqrt_rho_vqte @ rho_exact @ sqrt_rho_vqte 
        sqrt_product = scipy.linalg.sqrtm(product)
        fidelity_value = np.real(np.trace(sqrt_product))**2
        fidelities.append(1-fidelity_value)

    return fidelities

def extract_density_matrix_components(vqte_list, exact_list):
    """
    Extracts the four components (a, b, c, d) of 2x2 density matrices
    for every time step.

    Parameters:
    results_list (list of np.ndarray): List of 2x2 density matrices (rho).

    Returns:
    dict: Dictionary containing lists for 'a', 'b', 'c', 'd' components.
    """
    va_list, vb_list, vc_list, vd_list = [], [], [], []
    ea_list, eb_list, ec_list, ed_list = [], [], [], []

    for rho in vqte_list:
        va_list.append(rho[0, 0])
        vb_list.append(rho[0, 1]) # Note: b and c should generally be complex conjugates
        vc_list.append(rho[1, 0])
        vd_list.append(rho[1, 1])
        
    

    for rho in exact_list:
        ea_list.append(rho[0, 0])
        eb_list.append(rho[0, 1]) # Note: b and c should generally be complex conjugates
        ec_list.append(rho[1, 0])
        ed_list.append(rho[1, 1])

    print('VQTE A List:', va_list)
    print('Exact A List:', ea_list)
    return {'va': va_list, 'vb': vb_list, 'vc': vc_list, 'vd': vd_list, 'ea': ea_list, 'eb': eb_list, 'ec': ec_list, 'ed': ed_list}



def plot_matrix_components(all_matrix_components_by_layer, time, nt, layers):

    plt.style.use('seaborn-v0_8-talk')
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    time_axis = np.linspace(0, time, nt + 1)
    

    titles = [
        r'Diagonal Element 00 (Occupation of A)',
        r'Diagonal Element 01 (Occupation of B)',
        r'Diagonal Element 10 (Occupation of C)', 
        r'Diagonal Element 11 (Occupation of D)',
    ]
    # Keys for extraction: VQTE (v) and Exact (e)
    component_v_keys = ['va', 'vb', 'vc', 'vd']
    component_e_keys = ['ea', 'eb', 'ec', 'ed'] 


    num_layers = len(all_matrix_components_by_layer)
    vqte_colors = plt.cm.plasma(np.linspace(0.2, 0.9, num_layers))
    
    # Plotting loop
    for i, data in enumerate(all_matrix_components_by_layer):
  

        
        # Plot VQTE lines (solid, colored based on layer)
        for j in range(4):
            ax = axes[j]
            v_key = component_v_keys[j]
            # --- USE .get() with a fallback to avoid KeyError ---
            component_data = data.get(v_key, [])
            

            num_points = len(component_data)
            ax.plot(time_axis[:num_points], 
                    component_data, 
                    label=f'VQTE, {i+1} Layers', 
                    color=vqte_colors[i],
                    linestyle='-',
                    linewidth=2,
                    alpha=0.7)

        # Plot Exact lines (dashed, black/dark gray, only once per subplot)
        if i == 0:
            for j in range(4):
                ax = axes[j]
                e_key = component_e_keys[j]
                component_data = data.get(e_key, []) 
                
                if not component_data:
                     print(f"Warning: No data found for Exact key '{e_key}'.")
                     continue

                num_points = len(component_data)
                ax.plot(time_axis[:num_points], 
                        component_data, 
                        label=f'Exact', 
                        color='black',
                        linestyle='--',
                        linewidth=3,
                        alpha=0.9)


    # Finalize plots
    for ax, title in zip(axes, titles):
        ax.set_title(title, fontsize=18, pad=10)
        ax.set_xlabel('Time (t)', fontsize=16)
        ax.set_ylabel('Real Value', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(title='Method/Layers', fontsize=12, loc='best')

    plt.tight_layout()
    plt.show()
    
    
    
def plot_multiple_fidelity_vs_layers(fidelity_results, time, nt):
    plt.style.use('seaborn-v0_8-talk')
    fig, ax = plt.subplots(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt + 1)
    
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(fidelity_results)))
    
    for layer_idx, fidelities in enumerate(fidelity_results):
        num_points = len(fidelities)
        ax.plot(time_axis[:num_points], 
                fidelities, 
                label=f'{layer_idx+1} Layers', 
                color=colors[layer_idx],
                linewidth=2)
    
    ax.set_title("Fidelity vs Time for Different Ansatz Layers", fontsize=20, pad=20)
    ax.set_xlabel("Time", fontsize=16)
    ax.set_ylabel("Fidelity", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, alpha=0.3)
    # ax.set_ylim(0.0, 0.05)
    ax.legend(title='Number of Layers', fontsize=12)
    
    plt.tight_layout()
    plt.show()