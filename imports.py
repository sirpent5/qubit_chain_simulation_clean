from exact_methods import build_exact_diag_hamiltonian, perform_exact_diag
from global_methods import  build_initial_states, output_results
from vqte_methods import hamiltonian_generation, perform_vqte

from qiskit.circuit.library import EfficientSU2
from qiskit.quantum_info import Statevector
from qiskit_algorithms.time_evolvers.variational import RealMcLachlanPrinciple, ImaginaryMcLachlanPrinciple
from qiskit_algorithms import TimeEvolutionProblem, VarQRTE, VarQITE
from qiskit_algorithms.gradients import ReverseEstimatorGradient, ReverseQGT, DerivativeType
from qiskit.quantum_info import SparsePauliOp

import numpy as np
import matplotlib.pyplot as plt
import scipy  

import numpy as np
Sigma_x = np.matrix([[0, 1], [1, 0]])
Sigma_y = np.matrix([[0, -1j], [1j, 0]])
Sigma_z = np.matrix([[1, 0], [0, -1]])
Sigma_plus = (Sigma_x +1j*Sigma_y)/2
Sigma_minus = Sigma_plus.getH()
numberop = Sigma_minus@Sigma_plus