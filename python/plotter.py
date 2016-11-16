import numpy as np
import matplotlib.pyplot as plt
import sys

def the_plotter(temp, L, T):
    data = np.loadtxt("../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T%1.1f.xyz" %float(T), skiprows=1)
    L = float(L)
    MCC = np.zeros(len(data))
    T = np.zeros(len(data))
    Eexp = np.zeros(len(data))
    Mexp = np.zeros(len(data))
    sigma_E = np.zeros(len(data))
    sigma_M = np.zeros(len(data))
    accepted = np.zeros(len(data))
    
    for i in range(0, len(data)):
        MCC[i] = data[i][0]
        T[i] = data[i][1]
        Eexp[i] = data[i][2]/L/L
        Mexp[i] = data[i][3]/L/L
        sigma_E[i] = data[i][4]
        sigma_M[i] = data[i][5]
        accepted[i] = data[i][6]
   
    x = []
    for i in range(0, 7):
        x.append(10**i)
    print "Analytical %1.5f %1.5f %1.5f %1.5f" %(Energy_eig(1, 1)/L/L, Mag_eig(1, 1)/L/L, heatcap(1,1)/L/L, suscept(1,1)/L/L)
    for i in range(0, len(MCC)):
        for j in x:
            if MCC[i] == j:
                print "%5.5f %1.5f %1.5f %1.5f %1.5f" %(MCC[i], Eexp[i], Mexp[i], sigma_E[i], sigma_M[i])
                

    if str(temp) == "E" or str(temp) == "all":
        plt.figure()
        plt.plot(MCC, Eexp, label = "Numerical <E>")
        #plt.plot(MCC, Energy_eig(1,1)*MCC/MCC/L/L, label = 'Analytical')
        plt.title('Energy eigenvalues at T = %d' %T[1])
        plt.legend(loc = 'best')
        plt.ylim(-1.99, -2.002)
        plt.xlim(-5000, 500000)
        plt.savefig("../figures/task_b/energyeig_T%s.png" %(T[1]))
        plt.show()

    if str(temp) == "M" or str(temp) == "all":
        plt.figure()
        plt.plot(MCC, Mexp, label = "<M>")
        #plt.plot(MCC, Mag)
        plt.title('Magnetization eigenvalues at T = %d' %T[1])
        plt.legend(loc = 'best')
        #plt.ylim(-1.99, -2.002)
        plt.xlim(-5000, 500000)
        plt.savefig("../figures/task_b/Mageig_T%s.png" %(T[1]))
        plt.show()

    if str(temp) == "sE" or str(temp) == "all":
        plt.figure()
        plt.plot(MCC, sigma_E, label = r"$\sigma_E^2$")
        plt.title(r'$\sigma_E^2$ at T = %d' %T[1])
        plt.legend(loc = 'best')
        #plt.ylim(-1.99, -2.002)
        plt.xlim(-5000, 500000)
        plt.savefig("../figures/task_b/sigE_T%s.png" %(T[1]))
        plt.show()

    if str(temp) == "sM" or str(temp) == "all":
        plt.figure()
        plt.plot(MCC, sigma_M, label = r"$\sigma_M^2$")
        plt.title(r'$\sigma_M^2$ at T = %d' %T[1])
        plt.legend(loc = 'best')
        #plt.ylim(-1.99, -2.002)
        plt.xlim(-5000, 500000)
        plt.savefig("../figures/task_b/sigM_T%s.png" %(T[1]))
        plt.show()

    if str(temp) == "acc" or str(temp) == "all":
        plt.figure()
        plt.plot(MCC, accepted, label = "Accepted")
        plt.title('Accepted states at T = %d' %T[1])
        plt.legend(loc = 'best')
        #plt.ylim(-1.99, -2.002)
        plt.xlim(-5000, 500000)
        plt.savefig("../figures/task_b/accepted_T%s.png" %(T[1]))
        plt.show()

def Z(beta, J):
    return 2*np.exp(beta*8*J) + 2*np.exp(-beta*8*J) + 12

def Energy_eig(beta, J):
    return -(32.0*J)/Z(beta, J) * np.sinh(beta*8*J)

def Mag_eig(beta, J):
    return 1.0/Z(beta, J) * (8*np.exp(beta*8*J) + 16)

def E2(beta, J):
    return (256.0 * J**2)/Z(beta, J) * np.cosh(beta*8*J)

def M2(beta, J):
    return 32.0/Z(beta, J) * (np.exp(beta*8*J) + 1)

def heatcap(beta, J):
    return (E2(beta, J) - Energy_eig(beta, J)**2)/(1*1**2)

def suscept(beta, J):
    return (M2(beta, J) - Mag_eig(beta, J)**2)/(1*1**2)

def reader(filename, L):
    L = float(L)
    data = np.loadtxt(filename, skiprows = 1)
    MCC = np.zeros(len(data))
    T = np.zeros(len(data))
    Eexp = np.zeros(len(data))
    Mexp = np.zeros(len(data))
    sigma_E = np.zeros(len(data))
    sigma_M = np.zeros(len(data))
    accepted = np.zeros(len(data))
    
    for i in range(0, len(data)):
        MCC[i] = data[i][0]
        T[i] = data[i][1]
        Eexp[i] = data[i][2]/L/L
        Mexp[i] = data[i][3]/L/L
        sigma_E[i] = data[i][4]
        sigma_M[i] = data[i][5]
        accepted[i] = data[i][6]

    return MCC, T, Eexp, Mexp, sigma_E, sigma_M, accepted

def the_plotterc(temp):
    T1rand = reader("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.0.xyz", 20)
    T1up = reader("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.0.xyz", 20)
    T24up = reader("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.4.xyz", 20)
    T24rand = reader("../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.4.xyz", 20)

    values1 = [T1rand, T1up]
    values2 = [T24rand, T24up]

    #################################
    ####           E             ####
    #################################
    if str(temp) == "E" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.subplot(2,1,1)

        counter = 0
        for i in values1:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[2], label = 'Random')
            else:
                plt.plot(i[0], i[2], label = 'Up')
        
        plt.title('T = 1')
        plt.legend(loc = 'best')
        plt.ylabel('<E>', fontsize = 16)
        plt.ylim(-1.999, -1.98)
        plt.xlim(-5000, 100000)

        plt.subplot(2, 1, 2)

        counter = 0
        for i in values2:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[2], label = 'Random')
            else:
                plt.plot(i[0], i[2], label = 'Up')
        
        plt.title('T = 2.4')
        plt.legend(loc = 'best')
        #plt.ylim(-2.01, -1.96)
        plt.xlim(-5000, 100000)
        plt.ylabel('<E>', fontsize = 16)
        plt.xlabel('Monte Carlo cycles', fontsize = 16)
        plt.savefig("../figures/task_c/energyeig.png")
        plt.show()

    #################################
    ####           M             ####
    #################################
    if str(temp) == "M" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.subplot(2,1,1)

        counter = 0
        for i in values1:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[3], label = 'Random')
            else:
                plt.plot(i[0], i[3], label = 'Up')
        
        plt.title('Magnetic eigenvalues at T = 1')
        plt.legend(loc = 'best')
        plt.ylabel('<M>', fontsize = 16)
        plt.ylim(0.95, 1.05)
        plt.xlim(-5000, 100000)

        plt.subplot(2, 1, 2)

        counter = 0
        for i in values2:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[3], label = 'Random')
            else:
                plt.plot(i[0], i[3], label = 'Up')
        
        plt.title('T = 2.4')
        plt.legend(loc = 'best')
        #plt.ylim(-2.01, -1.96)
        plt.xlim(-5000, 100000)
        plt.ylabel('<M>', fontsize = 16)
        plt.xlabel('Monte Carlo cycles')
        plt.savefig("../figures/task_c/Mageig.png")
        plt.show()


    #################################
    ####           sE            ####
    #################################
    if str(temp) == "sE" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.subplot(2,1,1)

        counter = 0
        for i in values1:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[4], label = 'Random')
            else:
                plt.plot(i[0], i[4], label = 'Up')
        
        plt.title('T = 1')
        plt.legend()
        plt.ylabel(r'$\sigma_E^2$', fontsize = 20)
        plt.ylim(-1, 3)
        plt.xlim(-5000, 100000)

        plt.subplot(2, 1, 2)

        counter = 0
        for i in values2:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[4], label = 'Random')
            else:
                plt.plot(i[0], i[4], label = 'Up')
        
        plt.title('T = 2.4')
        plt.legend(loc = 'best')
        #plt.ylim(-2.01, -1.96)
        plt.xlim(-5000, 100000)
        plt.ylabel(r'$\sigma_E^2$', fontsize = 20)
        plt.xlabel('Monte Carlo cycles', fontsize = 16)
        plt.savefig("../figures/task_c/sigmaE.png")
        plt.show()

    
    #################################
    ####           sM            ####
    #################################
    if str(temp) == "sM" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.subplot(2,1,1)

        counter = 0
        for i in values1:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[5], label = 'Random')
            else:
                plt.plot(i[0], i[5], label = 'Up')

        plt.title('T = 1')
        plt.legend(loc = 'best')
        plt.ylabel(r'$\sigma_M^2$', fontsize = 20)
        plt.ylim(-1, 4)
        plt.xlim(-5000, 100000)

        plt.subplot(2, 1, 2)

        counter = 0
        for i in values2:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[5], label = 'Random')
            else:
                plt.plot(i[0], i[5], label = 'Up')
        
        plt.title('T = 2.4')
        plt.legend(loc = 'best')
        #plt.ylim(-2.01, -1.96)
        plt.xlim(-5000, 100000)
        plt.ylabel(r'$\sigma_M^2$', fontsize = 20)
        plt.xlabel('Monte Carlo cycles', fontsize = 16)
        plt.savefig("../figures/task_c/sigmaM.png")
        plt.show()

    #################################
    ####        Accepted         ####
    #################################
    if str(temp) == "acc" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.subplot(2,1,1)

        counter = 0
        for i in values1:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[6], label = 'Random')
            else:
                plt.plot(i[0], i[6], label = 'Up')
        
        plt.title('T = 1')
        plt.legend(loc = 'best')
        plt.ylabel('Accepted', fontsize = 16)
        #plt.ylim(-1.999, -1.98)
        #plt.xlim(-5000, 100000)

        plt.subplot(2, 1, 2)

        counter = 0
        for i in values2:
            counter += 1
            if counter == 1:
                plt.plot(i[0], i[6], label = 'Random')
            else:
                plt.plot(i[0], i[6], label = 'Up')
        
        plt.title('T = 2.4')
        plt.legend(loc = 'best')
        #plt.ylim(-2.01, -1.96)
        #plt.xlim(-5000, 100000)
        plt.ylabel('Accepted', fontsize = 16)
        plt.xlabel('Monte Carlo cycles', fontsize = 16)
        plt.savefig("../figures/task_c/accepted.png")
        plt.show()

def kindle(filename):
    "Because it's the energy reader, so this is a E-reader"
    counter = []
    values = []

    E = np.loadtxt(filename, skiprows = 1)

    for i in E:
        counter.append(i[0])
        values.append(i[1])
        
    return counter, values

def Energy_arr(fixerup):
    E = np.array([])
    for i in range(0, len(fixerup[0])):
        E = np.append(E, np.repeat(fixerup[1][i], fixerup[0][i]))
    return E

def Energy_arr2(d):
    E = []
    for i in range(0, len(d[0])):
        for j in range(0, d[0][i]):
            E.append(d[1][i])
    return E

def gaussian(dX, sig2):
    print dX**2
    print 2*sig2
    print dX**2/(2*sig2)
    return np.exp(- dX**2/(2*sig2))
        
def d_solver():
    T1rand = reader("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim20_dirrand_T1.0.xyz", 20)
    T1up = reader("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim20_dirup_T1.0.xyz", 20)
    T24up = reader("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim20_dirup_T2.4.xyz", 20)
    T24rand = reader("../benchmarks/task_d/eigenvalues_MC1.0e+06_dim20_dirrand_T2.4.xyz", 20)

    E1rand = kindle("../benchmarks/task_d/EnergyValues_T1.0_L20.0_dir_rand.txt")
    E1up = kindle("../benchmarks/task_d/EnergyValues_T1.0_L20.0_dir_up.txt")
    E24up = kindle("../benchmarks/task_d/EnergyValues_T2.4_L20.0_dir_up.txt")
    E24rand = kindle("../benchmarks/task_d/EnergyValues_T2.4_L20.0_dir_rand.txt")

    E1randE = Energy_arr(E1rand)
    E1upE = Energy_arr(E1up)
    E24upE = Energy_arr(E24up)
    E24randE = Energy_arr(E24rand)

    E1upE2 = Energy_arr(E1up)
    print len(E1upE2)

    print np.var(E24randE)/20**2

    print len(E1upE)
    
    n, bins, patch = plt.hist(E24upE, bins = 20, normed = 1)
    amp = np.max(n)

    x = np.linspace(bins[0], bins[-1], 1001)
    print T24up[2][-1], T24up[4][-1]
    analytical = amp*gaussian(x - T24up[2][-1], T24up[4][-1])

    print analytical
    plt.plot(x, analytical, 'g--', label = 'analytical')
    plt.legend(loc = 'best')

    plt.show()

    """
    plt.figure(figsize = (10, 8))
    plt.subplot(211)
    plt.title('T = 1.0')
    plt.bar(E1rand[1], E1rand[0]/np.sum(E1rand[0]), label = 'Rand')
    plt.ylabel('$P(E)$', fontsize = 16)
    plt.legend(loc = 'best')

    plt.subplot(212)
    plt.bar(E1up[1], E1up[0]/np.sum(E1up[0]), label = 'Up')
    plt.ylabel('$P(E)$', fontsize = 16)
    plt.xlabel('Energies', fontsize = 16)
    plt.legend(loc = 'best')
    plt.savefig("../figures/task_d/hist1.png")

    barwidth = 4
    print len(E24rand[1])
    E24rand = binsizer(E24rand, binsize = 10)
    print len(E24rand[1])

    plt.figure(figsize = (10, 8))
    plt.subplot(211)
    plt.title('T = 2.4')
    plt.bar(E24rand[1], E24rand[0]/np.sum(E24rand[0]), barwidth, label = 'Rand')
    plt.ylabel('$P(E)$', fontsize = 16)
    plt.legend(loc = 'best')

    plt.subplot(212)
    plt.bar(E24up[1], E24up[0]/np.sum(E24up[0]), barwidth, label = 'up')
    #plt.plot(E24up[1], gaussian(E24up[1], E24up[0], T24up[4][-1], np.mean(E24up[1])), '*r' ,label = 'Gauss')
    plt.ylabel('$P(E)$', fontsize = 16)
    plt.xlabel('Energies', fontsize = 16)
    plt.legend(loc = 'best')
    plt.savefig("../figures/task_d/hist24.png")

    print "Var(E) from program for T = 1.0 is ", T1rand[4][-1]
    print "Var(E) from program for T = 2.4 is ", T24rand[4][-1]

    print "Var(E) numerically for T = 1.0 is  ", np.var(E1rand[1])/20**2
    print "Var(E) numerically for T = 2.4 is  ", np.var(E24rand[1])/20**2

    """
    #print gaussian(E24up[1], E24up[0], T24up[4][-1], np.mean(E24up[1]))
    
    #plt.show()
    
    
def binsizer(fixerupper, binsize = 10.0):
    counter = 1
    count_size = 0
    energies_size = 0
    counter1 = []
    values1 = []
    
    for i in range(0 ,len(fixerupper[0])):
        #print i, counter
        if counter < binsize:
            energies_size += fixerupper[1][i]
            count_size += fixerupper[0][i]
            counter += 1
        elif counter == binsize:
            values1.append(energies_size/binsize)
            counter1.append(count_size/binsize)
            counter = 0
            count_size = 0
            energies_size = 0

    return counter1, values1


def the_plottere(temp):
    L20  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim20_dirrand_dt0.03333.xyz', 20)
    L40  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim40_dirrand_dt0.01034.xyz', 40)
    L60  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim60_dirrand_dt0.01034.xyz', 60)
    L80  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim80_dirrand_dt0.01034.xyz', 80)
    L100 = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim100_dirrand_dt0.01034.xyz', 100)
    L140 = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim140_dirrand_dt0.00612.xyz', 140)

    values = [L20, L40, L60, L80, L100, L140]
    
    #################################
    ####           E             ####
    #################################
    if str(temp) == "E" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Unknown')

        counter = 20
        for i in values:
            plt.plot(i[1], i[2], label = 'L%d' %counter)

            if counter == 100:
                counter += 40
            else:
                counter += 20
        
        plt.legend(loc = 'best')
        plt.ylabel('<E>', fontsize = 16)
        plt.xlabel('Temperature', fontsize = 16)
        plt.legend(loc='best' )
        plt.savefig("../figures/task_e/energyeig.png")
        plt.show()

    #################################
    ####           M             ####
    #################################
    if str(temp) == "M" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Unknown')

        counter = 20
        for i in values:
            plt.plot(i[1], i[3], label = 'L%d' %counter)

            if counter == 100:
                counter += 40
            else:
                counter += 20
        
        plt.legend(loc = 'best')
        plt.ylabel('<M>', fontsize = 16)
        plt.xlabel('Temperature', fontsize = 16)
        plt.legend(loc='best' )
        plt.savefig("../figures/task_e/Mageig.png")
        plt.show()


    #################################
    ####           sE            ####
    #################################
    if str(temp) == "sE" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Unknown')

        counter = 20
        for i in values:
            plt.plot(i[1], i[4], label = 'L%d' %counter)

            if counter == 100:
                counter += 40
            else:
                counter += 20
        
        plt.legend(loc = 'best')
        plt.ylabel(r'$\sigma_E^2$', fontsize = 20)
        plt.xlabel('Temperature', fontsize = 16)
        plt.legend(loc='best')
        plt.savefig("../figures/task_e/sigmaE.png")
        plt.show()

    
    #################################
    ####           sM            ####
    #################################
    if str(temp) == "sM" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Unknown')

        counter = 20
        for i in values:
            plt.plot(i[1], i[6], label = 'L%d' %counter)

            if counter == 100:
                counter += 40
            else:
                counter += 20
        
        plt.legend(loc = 'best')
        plt.ylabel('Accepted', fontsize = 16)
        plt.xlabel('Temperature', fontsize = 16)
        plt.legend(loc='best' )
        plt.savefig("../figures/task_e/accepted.png")
        plt.show()



if __name__ == '__main__':
    if sys.argv[1] == "b":
        the_plotter(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "c":
        the_plotterc(sys.argv[2])
    elif sys.argv[1] == "d":
        d_solver()
    elif sys.argv[1] == "e":
        the_plottere(sys.argv[2])
