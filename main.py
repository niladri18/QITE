from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import Aer, execute
from qiskit import IBMQ, transpile, assemble
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.aer import AerSimulator
from qiskit.providers.aer.noise import NoiseModel
from qiskit.opflow import X, Y, Z, I, CX, T, H, S, PrimitiveOp
import numpy as np
from qiskit.opflow import I, X, Y, Z, H, CX, Zero, ListOp, PauliExpectation, PauliTrotterEvolution, CircuitSampler, MatrixEvolution, Suzuki
from qiskit.circuit import Parameter
from qiskit.quantum_info import Pauli
from qiskit import Aer
from qiskit.primitives import Estimator
from qiskit.opflow import I, X, Z, PauliSumOp
from qiskit.opflow import Trotter, PauliTrotterEvolution, PauliSumOp
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit import Parameter
from qiskit.opflow import StateFn, AerPauliExpectation, CircuitSampler
from qiskit.utils import QuantumInstance
from qiskit_aer import AerSimulator
from qiskit import Aer, execute
import pdb
import re
import sys
import multiprocessing
from joblib import Parallel, delayed
from joblib import parallel_backend

from qiskit_ibm_runtime import QiskitRuntimeService, Session,  Options
from qiskit.providers.fake_provider import FakeKolkata, FakeWashington
import json


'''
For DD
'''
from qiskit.circuit.library import XGate
from qiskit.transpiler import PassManager, InstructionDurations
from qiskit.transpiler.passes import ALAPScheduleAnalysis, PadDynamicalDecoupling

backend = FakeKolkata()


'''
mthree
'''
import mthree
mit = mthree.M3Mitigation(backend)

#backend = provider.get_backend('ibmq_kolkata')
#noise_model = NoiseModel.from_backend(backend)
#device = AerSimulator.from_backend(device_backend)


num_cores = multiprocessing.cpu_count()
options = Options()
options.resilience_level = 1
options.execution.shots = 10000

IBMQ.load_account()
qasm_sim = Aer.get_backend('qasm_simulator')
#provider = IBMQ.get_provider(hub='ibm-q-ornl',group='lbnl', project='chm170')
service = QiskitRuntimeService(channel="ibm_quantum")
#backend = service.backend("ibmq_qasm_simulator")
#backend = FakeKolkata()
#backend = FakeWashington()
#backend = provider.backend("ibmq_qasm_simulator")
#backend = provider.backends.ibmq_qasm_simulator


def process_all_result(result, Nq):

    states = result.keys()
    prob = [0,0]
    for state in states:
        qstate = state.split()[0]
        if qstate[-1] == "0":
            prob[0] += result[state]
        else:
            prob[1] += result[state]

    return prob

def  find_prob(prob_wf):
    n = len(prob_wf)
    p0 = 0
    for i in range(n):
        if i%2 == 0:
            p0 += prob_wf[i]
    return p0


def build_init_state(init_circuit, state):
    for i,s in enumerate(list(state)):
        if s=="1":
            init_circuit.x(i)
        if s=="+":
            init_circuit.h(i)

    return init_circuit

def makeTrotter(hamiltonian, t, reps=1):

    # reps: no of Trotter steps

    evo_op = hamiltonian.exp_i()
    #print(evo_op)
    #time = Parameter('t')
    time = t


    evolution = PauliTrotterEvolution(trotter_mode=Trotter(), reps=reps)
    evol_result = evolution.convert(time*hamiltonian.exp_i())
    evolved_state = evol_result.to_circuit()
    return evolved_state

def buildHamiltonian(fname, nq):

    #nq: number of qubits

    f0 = open(fname,"r")
    hterm = []
    ham_list = []
    while True:
        line = f0.readline().strip()
        if len(line) < 1:
            break
        pstrng = ["I"]*nq
        line = line.split("*")
        coeff = eval(line[0])
        #if len(line) < 2:
        #    pstrng = "".join(pstrng)
        #    hterm.append( ( pstrng, coeff) )
        #    continue
        pa_string = line[1]
        palist = re.findall('[A-Z]',pa_string)
        qalist = [int(nu) for nu in re.findall(r'\d+',pa_string)]
        for (p,q) in zip(palist, qalist):
            pstrng[nq-q-1] = p
        pstrng = "".join(pstrng)
        ham_list.append( [( palist, qalist), coeff ])
        hterm.append( ( pstrng, coeff))

        

    return hterm, ham_list



def buildOperator(fname, nq):

    #nq: number of qubits

    f0 = open(fname,"r")
    oterm = []
    op_list = []
    while True:
        line = f0.readline().strip()
        if len(line) < 1:
            break
        pstrng = ["I"]*nq
        line = line.split("*")
        coeff = eval(line[0])
        #if len(line) < 2:
        #    pstrng = "".join(pstrng)
        #    oterm.append( ( pstrng, coeff) )
        #    continue
        pa_string = line[1]
        palist = re.findall('[A-Z]',pa_string)
        qalist = [int(nu) for nu in re.findall(r'\d+',pa_string)]
        for (p,q) in zip(palist, qalist):
            pstrng[nq-q-1] = p
        pstrng = "".join(pstrng)
        op_list.append( ( palist, qalist) )
        #op_list.append( [( palist, qalist), coeff ])
        oterm.append( ( pstrng, coeff))

        
    f0.close()

    return oterm, op_list


def getOperatorList(fname,nq):

    # Get the operators to be measured
    f0 = open(fname,"r")
    operator_list = []
    while True:
        line = f0.readline().strip()
        print(line)

        if len(line) < 1:
            break
        op = buildOperator(line,nq)
        operator_list.append(op)

    return operator_list


#def add_controlled_P_circuit(circ,operator):
def add_controlled_P_circuit(circ,op,qb):



     if (op == 'X'):
         circ.cx(0,qb+1)
     if (op == 'Y'):
         circ.cy(0,qb+1)
     if (op == 'Z'):
         circ.cz(0,qb+1)

def measure_expectation_hadamard_qasm(new_state, operator, im):
    '''
    Uses Hadamard test to compute expectation values
    '''

    shots = 10000

    # check if identity
    if all(element == "I" for element in list(operator[0])):
        return operator[1], {}


    Nq = len(operator[0]) 
    p = QuantumCircuit(Nq+1) # one ancilla qubit
     

    ''' 
    #Nq = 1
    q = QuantumCircuit(Nq)
    p.h(0)
    q.x(0)
    p = p.compose(q, qubits = [i for i in range(1,Nq+1)])
    p.cx(0,1)
    p.h(0)
    pdb.set_trace()
    '''


    p.h(0)
    if im:
        p.sdg(0)

    p = p.compose(new_state, qubits = [i for i in range(1,Nq+1)])

    add_controlled_P_circuit(p,operator[0])

    p.h(0)

    qasm = True
    if qasm:

        p.measure_all()

        p_trans = transpile(p, backend = qasm_sim, optimization_level=3)

        result = execute(p_trans, backend=qasm_sim, shots=shots,\
                     optimization_level=3).result().get_counts()


        prob = process_all_result(result, Nq+1)
        p0 = prob[0]/sum(prob)
    else:
        simulator = Aer.get_backend('statevector_simulator')
        result = execute(p, simulator).result()
        wf = result.get_statevector()
        prob = (abs(np.array(wf)))**2
        p0 = find_prob(prob)

    return (operator[1]*(2*p0-1)), result



def measure_expectation_hadamard(new_state, operator, Nq, im):
    '''
    Uses Hadamard test to compute expectation values
    '''

    shots = 10000



    p = QuantumCircuit(Nq+1) # one ancilla qubit
     



    p.h(0)
    if im:
        p.sdg(0)

    p = p.compose(new_state, qubits = [i for i in range(1,Nq+1)])

    for pauli, qb in zip(operator[0],operator[1]):
        add_controlled_P_circuit(p,pauli,qb)

    p.h(0)
    #pdb.set_trace()

    qasm = False
    if qasm:

        p.measure_all()

        # Transpile the circuit 'qc' 10 times
        trans_circ_list = transpile([p]*10, backend = backend, optimization_level=3)

        # get the number of cnot gates
        cx_counts = np.asarray([circ.count_ops().get('cx', 0) for circ in trans_circ_list])

        # select the circuit with the lowest number
        best_idx = np.argmin(cx_counts)
        best_trans_qc = trans_circ_list[best_idx]



        #trans_circ = transpile(p, backend = qasm_sim, optimization_level=3)
        #cx_counts = trans_circ.count_ops().get('cx', 0)
        print("CNOTS:",cx_counts)


        durations = InstructionDurations.from_backend(backend)
        constraints = backend.configuration().timing_constraints
        dd_sequence = [XGate(), XGate()]
        pm = PassManager([ALAPScheduleAnalysis(durations),
                  PadDynamicalDecoupling(durations, dd_sequence,
                  pulse_alignment=constraints['pulse_alignment'])])
        opt_qc = pm.run(best_trans_qc)




        raw_counts = execute(opt_qc, backend=backend, shots=shots,\
                     optimization_level=3).result().get_counts()

        measured_qubits = mthree.utils.final_measurement_mapping(opt_qc)

        mit.cals_from_system(measured_qubits)
        noisy_result = mit.apply_correction(raw_counts, measured_qubits)



        prob = process_all_result(noisy_result, Nq+1)
        p0 = prob[0]/sum(prob)
    else:
        simulator = Aer.get_backend('statevector_simulator')
        result = execute(p, simulator).result()
        wf = result.get_statevector()
        prob = (abs(np.array(wf)))**2
        p0 = find_prob(prob)
        noisy_result = result


    return (2*p0-1)


def compute_expectation_operator(operator, new_state):

    op_sum = 0
    result_noisy_op = [] 
    for op_j in operator:
        pdb.set_trace()
        #expectation_value, result = measure_expectation_hadamard_qasm(new_state,op_j)
        expectation_value, result = measure_expectation_hadamard(new_state,op_j)
        #expectation_value = measure_expectation_hadamard_qasm(new_state,op_j)
        
        op_sum += expectation_value
        result_noisy_op.append(result)

    return op_sum, result_noisy_op

def compute_expectation_series(new_state, operator_list): 

    expectation_list = []
    result_noisy = []
    for operator in operator_list:
        op_sum, result_op = compute_expectation_operator(operator, new_state)
        result_noisy = result_noisy + result_op
        expectation_list.append(op_sum)
        print(op_sum)

    return expectation_list, result_noisy


def compute_smat(new_state, operator_list, Nq): 

    expectation_list = []
    result_noisy = []
    n = len(operator_list)
    smat = np.eye(n)
    for i,op1 in enumerate(operator_list):
        for j,op2 in enumerate(operator_list):
            if i >= j:
                continue
            operator = (op1[0]+op2[0],op1[1]+op2[1] )

            #op_sum, result_op = compute_expectation_operator(operator, new_state)
            expectation_value = measure_expectation_hadamard(new_state,operator, Nq, False)
            smat[i,j] = expectation_value
            smat[j,i] = smat[i,j] 

        #result_noisy = result_noisy + result_op
        #expectation_list.append(op_sum)
        #print(op_sum)

    return smat


def compute_bvec(new_state, operator_list, ham_term, Nq): 

    expectation_list = []
    result_noisy = []
    n = len(operator_list)
    bvec = np.zeros(n, dtype = complex)
    for i,op1 in enumerate(operator_list):
        term = 0
        for j,op2 in enumerate(ham_term):
            operator = (op1[0]+op2[0][0],op1[1]+op2[0][1] )
            #print(op1,op2)
            coeff = op2[1]

            expectation_value = measure_expectation_hadamard(new_state,operator, Nq, True)
            term -= 2*coeff*expectation_value

        #bvec[i] = 1j*term - 1j*np.conj(term)
        bvec[i] = term

        #result_noisy = result_noisy + result_op
        #expectation_list.append(op_sum)
        #print(op_sum)

    return bvec

def compute_expectation_parallel(new_state, operator_list): 


    expectation_list = []
    n_op = len(operator_list)

    print(num_cores)

    ans = Parallel(n_jobs=num_cores , prefer="threads")(delayed(compute_expectation_operator)(op, new_state) for op in operator_list)
    
    expectation_list = []
    result_noisy = []
    for item in ans:
        expectation_list.append(item[0])
        result_noisy = result_noisy + item[1]
        print(item[0])

    return expectation_list, result_noisy

def generate_operator_estimator(operator):
    '''
    Generate operator in appropriate form to 
    evaluate expectation using estimator
    '''
    pdb.set_trace()
    op_list = []
    for op_j in operator:
        op = SparsePauliOp.from_list([op_j])
        op_list.append(op)

    pdb.set_trace()
    return op_list

            
def compute_amat_(new_state, operator_list): 


    expectation_list = []
    n = len(operator_list)
    Amat = np.eye(n)
    for i in range(n):
        #op1 = generate_operator_estimator(operator_list[i])
        op1 = operator_list[i]
        for j in range(i+1,n):
            #op2 = generate_operator_estimator(operator_list[j])
            op2 = operator_list[j]
            operator = (op1[0]+op2[0],op1[1]*op2[1] )
            coeff = op1[1]*op2[1]
            op_list = []
            qb_list = []
            op = ""
            for i, (sigma1,sigma2) in enumerate(zip(op1[0],op2[0])):
                if sigma1 == "I" and sigma2 == "I":
                    continue
                if sigma1 == "I" and sigma2 != "I":
                    op += sigma2
                    qb_list.append(i) 
                if sigma1 != "I" and sigma2 == "I":
                    op += sigma1
                    qb_list.append(i) 
                if sigma1 != "I" and sigma2 != "I":
                    op += sigma1+sigma2
                    qb_list.append(i) 
                    qb_list.append(i) 

            pdb.set_trace()
            #qiskit_op = SparsePauliOp.from_list([operator])
            qiskit_op = SparsePauliOp.from_sparse_list([(op, qb_list, op1[1]*op2[1])], num_qubits = 3)
            estimator = Estimator()
            pdb.set_trace()
            tmp = estimator.run([new_state], qiskit_op).result().values.real
            pdb.set_trace()
            Amat[i,j] = tmp
            Amat[j,i] = tmp

    return Amat

def get_minv(a, delta=1.e-6):
    ap = a + delta*np.eye(a.shape[0])
    ainv = np.linalg.pinv(ap)
    return ainv

def compute_expectation_fast(new_state, operator_list): 


    expectation_value = 0

    with Session(service=service, backend="ibmq_qasm_simulator") as session:
        for operator in operator_list:

            estimator = Estimator()
            #job = estimator.run(new_state, op)
            #result = job.result()
            num_op = len(operator)
            expectation_value += estimator.run([new_state]*num_op, operator).result().values.real
    session.close()

    return expectation_value


if __name__=="__main__":

    nq = 2

    # load correlation operators
    op_fname = "n2_op.dat"
    oterm, operator_list = buildOperator(op_fname,nq)

    #load Hamiltonian
    fname = "N2_U1.0_J0.0_h0.0.dat"
    t = 0.1

    hterm, hamiltonian_list = buildHamiltonian(fname, nq)
    nterm  = len(hterm)
    print(nterm)
    hamiltonian = PauliSumOp.from_list(hterm )

    #hamiltonian = PauliSumOp.from_list([("IZI", 0.39), ("IXX", 0.5)])



    #load initial state
    state = "++" # declare the initial state
    init_circuit = QuantumCircuit(nq)
    init_circuit = build_init_state(init_circuit, state)
    #print(init_circuit)



    # time evolve state
    steps = 1
    T = 2.0
    dt = 0.1

    evolved_state = makeTrotter(hamiltonian, dt, steps) 

    ofname = f"run2_noisy_T{T}_circ_dt{dt}"+fname[:-4]+".dat"
    ftmp = open(ofname,"w")

    data_fname = f"run2_noisy_T{T}_circ_dt{dt}"+fname[:-4]+".rom"

    new_state = init_circuit
    all_data = []
    all_noisy_result = []
    for t in list(np.arange(0, T+dt,dt)):
        # add the initial state and the evolved state
        #if t == 0:
        #else:
        #    new_state = init_circuit.compose(evolved_state)
        #trans_state = transpile(new_state, backend, optimization_level=3, approximation_degree=0)



        for hterm in hamiltonian_list:

            expectation_value = compute_expectation_fast(new_state, hamiltonian) #energy
            norm = np.sqrt(1 - 2* dt*expectation_value)

            smat = compute_smat(new_state, operator_list, nq)
            Smat_ = smat + np.transpose(smat)
            bvec = compute_bvec(new_state, operator_list, [hterm], nq) 
            bvec /= norm

            amat = get_minv(Smat_).dot(bvec)

            pauli_terms = []
            for t1, t2 in zip(amat, oterm):
                pauli_terms.append((t2[0],t1)) 
            pauli_terms = PauliSumOp.from_list( pauli_terms )
            evolved_state = makeTrotter(pauli_terms, dt, steps) 
            new_state = init_circuit.compose(evolved_state)
            init_circuit = new_state

        data = [t, expectation_value[0]] 
        print(*data, sep = "\t")
        print(*data, sep = "\t", file = ftmp)

        ftmp.flush()
        #all_data.append(data)

    all_data = np.array(all_data)
    #np.savetxt(f"T{T}_circ_dt{dt}"+fname[:-4]+".dat",all_data)
    sys.exit()

