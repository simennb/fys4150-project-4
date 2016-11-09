import matplotlib.pyplot as plt
from numpy import array, var
import sys


def plotter(filename, direction, task):
    f = open(filename)

    #initialize the lists we need
    MCcycles  = []
    T         = []
    E         = []
    EE        = []
    M         = []
    MM        = []
    r_counter = []
    
    variables_name = ["T", "<E>", "sigma_E", "<|M|>", "sigma_M", "r_counter"]
    variables = [T, E, EE, M, MM, r_counter]
    for line in f:
        if str(line.split(' ')[0]) == "MCcycles":
            MCcycles.append(int(line.split(' ')[2].split('\n')[0]))
            writer = True
        for variable_name in variables_name:
            if writer == True and str(line.split(' ')[0]) == variable_name:
                variables[variables_name.index(variable_name)].append(float(line.split(' ')[2].split('\n')[0]))
        if variable_name == variables_name[-1]:
            wrtier = False
    for i in variables:
        plt.plot(array(MCcycles), array(i), label = variables_name[variables.index(i)])
        plt.title(variables_name[variables.index(i)])
        plt.xlabel('Monte Carlo Cycles')
        plt.ylabel('Eigenvalue')
        plt.legend(loc = 'best') 
        plt.savefig("../figures/task_%s/%s_T=%1.1f_dir%s.png" %(task, variables_name[variables.index(i)], T[-1], direction))
        plt.grid(True)
        plt.show()

def d_reader(filename, T, direction):
    f = open(filename)
    E_values = []
    E_counter = []

    counter = False
    values = False
    
    for line in f:
        if line.strip(' ').split(' ')[0] == "E_counter":
            counter = True
            values = False
            continue
        elif line.strip(' ').split(' ')[0] == "E_values":
            counter = False
            values = True
            continue
        
        if counter == True:
            E_counter.append(float(line.strip(' ').split(' ')[0].split('\n')[0]))
        if values == True:
            E_values.append(float(line.strip(' ').split(' ')[0].split('\n')[0]))

    E_counter_list = []
    E_values_list = []
    
    for i in range(0, len(E_values)):
        if E_counter[i] != 0.0:
            E_counter_list.append(E_counter[i])
            E_values_list.append(E_values[i])

    AE_counter = array(E_counter)#_list)
    AE_value = array(E_values)#_list)

    plt.figure(1)
    plt.bar(AE_value, AE_counter/float(sum(AE_counter))) #, align = 'center')
    plt.grid(True)
    plt.savefig("../figures/task_d/EvaluesBAR_T=%1.1f_dir%s.png" %(T, direction))
    plt.show()

    print var(AE_counter)
    print var(AE_counter/float(sum(AE_counter)))

if sys.argv[1] == "c" and sys.argv[2] == "up":
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirup_T1.0.xyz", "up", "c")
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirup_T2.4.xyz", "up", "c")
elif sys.argv[1] == "c" and sys.argv[2] == "rand":
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirrand_T1.0.xyz", "rand", "c")
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirrand_T2.4.xyz", "rand", "c")
elif sys.argv[1] == "c" and sys.argv[2] == "down":
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirdown_T1.0.xyz", "down", "c")
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirdown_T2.4.xyz", "down", "c")

if sys.argv[1] == "d" and sys.argv[2] == "up" and sys.argv[3] == 2:
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirup_T1.0.xyz", "up", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirup_T2.4.xyz", "up", "d")
elif sys.argv[1] == "d" and sys.argv[2] == "rand" and sys.argv[3] == 2:
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirrand_T1.0.xyz", "rand", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirrand_T2.4.xyz", "rand", "d")
elif sys.argv[1] == "d" and sys.argv[2] == "down" and sys.argv[3] == 2:
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirdown_T1.0.xyz", "down", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirdown_T2.4.xyz", "down", "d")

elif sys.argv[1] == "d" and sys.argv[2] == "rand" and sys.argv[3] == 20:
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim20_dirrand_T1.0.xyz", "rand", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim20_dirrand_T2.4.xyz", "rand", "d")

elif sys.argv[1] == "d" and sys.argv[2] == "rand" and sys.argv[3] == "Evalues":
    d_reader("../benchmarks/task_d/EnergyValues_T2.4_L20.0_dir_rand.txt", 2.4, "rand")
