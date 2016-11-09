import matplotlib.pyplot as plt
from numpy import array
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
    E_out = []
    
    variables_name = ["T", "<E>", "sigma_E", "<|M|>", "sigma_M", "r_counter", "E_out"]
    variables = [T, E, EE, M, MM, r_counter, E_out]
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
        plt.show()

if sys.argv[1] == "c" and sys.argv[2] == "up":
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirup_T1.0.xyz", "up", "c")
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirup_T2.4.xyz", "up", "c")
elif sys.argv[1] == "c" and sys.argv[2] == "rand":
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirrand_T1.0.xyz", "rand", "c")
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirrand_T2.4.xyz", "rand", "c")
elif sys.argv[1] == "c" and sys.argv[2] == "down":
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirdown_T1.0.xyz", "down", "c")
    plotter("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim2_dirdown_T2.4.xyz", "down", "c")

if sys.argv[1] == "d" and sys.argv[2] == "up":
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirup_T1.0.xyz", "up", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirup_T2.4.xyz", "up", "d")
elif sys.argv[1] == "d" and sys.argv[2] == "rand":
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirrand_T1.0.xyz", "rand", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirrand_T2.4.xyz", "rand", "d")
elif sys.argv[1] == "d" and sys.argv[2] == "down":
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirdown_T1.0.xyz", "down", "d")
    plotter("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim2_dirdown_T2.4.xyz", "down", "d")
