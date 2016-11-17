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
        T10r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.0.xyz', 2)
        T11r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.1.xyz', 2)
        T12r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.2.xyz', 2)
        T13r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.3.xyz', 2)
        T14r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.4.xyz', 2)
        T15r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.5.xyz', 2)
        T16r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.6.xyz', 2)
        T17r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.7.xyz', 2)
        T18r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.8.xyz', 2)
        T19r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T1.9.xyz', 2)
        T20r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.0.xyz', 2)
        T21r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.1.xyz', 2)
        T22r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.2.xyz', 2)
        T23r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.3.xyz', 2)
        T24r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.4.xyz', 2)
        T25r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.5.xyz', 2)
        T26r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.6.xyz', 2)
        T27r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.7.xyz', 2)
        T28r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.8.xyz', 2)
        T29r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T2.9.xyz', 2)
        T30r = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirrand_T3.0.xyz', 2)

        T10u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.0.xyz', 2)
        T11u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.1.xyz', 2)
        T12u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.2.xyz', 2)
        T13u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.3.xyz', 2)
        T14u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.4.xyz', 2)
        T15u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.5.xyz', 2)
        T16u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.6.xyz', 2)
        T17u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.7.xyz', 2)
        T18u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.8.xyz', 2)
        T19u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T1.9.xyz', 2)
        T20u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.0.xyz', 2)
        T21u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.1.xyz', 2)
        T22u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.2.xyz', 2)
        T23u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.3.xyz', 2)
        T24u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.4.xyz', 2)
        T25u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.5.xyz', 2)
        T26u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.6.xyz', 2)
        T27u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.7.xyz', 2)
        T28u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.8.xyz', 2)
        T29u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T2.9.xyz', 2)
        T30u = reader('../benchmarks/task_b/eigenvalues_MC1.0e+06_dim2_dirup_T3.0.xyz', 2)

        valuesr = [T10r, T11r, T12r, T13r, T14r, T15r, T16r, T17r, T18r, T19r, T20r, T21r, T22r, T23r, T24r, T25r, T26r, T27r, T28r, T29r, T30r]

        valuesu = [T10u, T11u, T12u, T13u, T14u, T15u, T16u, T17u, T18u, T19u, T20u, T21u, T22u, T23u, T24u, T25u, T26u, T27u, T28u, T29u, T30u]

        xlist = []
        ylist = []

        xlistu = []
        ylistu = []
        
        plt.figure(figsize = (10,8))
        for i in valuesr:
            xlist.append(i[1][-1])
            ylist.append(i[-1][-1])
        for i in valuesu:
            xlistu.append(i[1][-1])
            ylistu.append(i[-1][-1])
            
        plt.plot(xlist, ylist, '-g', label = "Random", lw = 2)
        plt.plot(xlistu, ylistu, '--b', label = "Up", lw = 2)
        plt.title('Total Accepted States for 1000000 Monte Carlo Cycles')
        plt.legend(loc = 'best')
        plt.xlabel("Temperatures", fontsize = 16)
        plt.ylabel("Numbers of total Accepted", fontsize = 16)
        plt.savefig("../figures/task_c/accepted.png")
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
        T10r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.0.xyz', 20)
        T11r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.1.xyz', 20)
        T12r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.2.xyz', 20)
        T13r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.3.xyz', 20)
        T14r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.4.xyz', 20)
        T15r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.5.xyz', 20)
        T16r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.6.xyz', 20)
        T17r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.7.xyz', 20)
        T18r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.8.xyz', 20)
        T19r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T1.9.xyz', 20)
        T20r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.0.xyz', 20)
        T21r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.1.xyz', 20)
        T22r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.2.xyz', 20)
        T23r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.3.xyz', 20)
        T24r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.4.xyz', 20)
        T25r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.5.xyz', 20)
        T26r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.6.xyz', 20)
        #T27r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.7.xyz', 20)
        T28r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.8.xyz', 20)
        T29r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T2.9.xyz', 20)
        T30r = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirrand_T3.0.xyz', 20)

        #T10u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.0.xyz', 20)
        #T11u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.1.xyz', 20)
        #T12u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.2.xyz', 20)
        #T13u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.3.xyz', 20)
        #T14u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.4.xyz', 20)
        #T15u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.5.xyz', 20)
        #T16u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.6.xyz', 20)
        #T17u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.7.xyz', 20)
        #T18u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.8.xyz', 20)
        #T19u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T1.9.xyz', 20)
        #T20u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.0.xyz', 20)
        #T21u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.1.xyz', 20)
        #T22u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.2.xyz', 20)
        #T23u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.3.xyz', 20)
        #T24u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.4.xyz', 20)
        #T25u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.5.xyz', 20)
        ##T26u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.6.xyz', 20)
        #T27u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.7.xyz', 20)
        #T28u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.8.xyz', 20)
        #T29u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T2.9.xyz', 20)
        #T30u = reader('../benchmarks/task_c/eigenvalues_MC1.0e+06_dim20_dirup_T3.0.xyz', 20)

        valuesr = [T10r, T11r, T12r, T13r, T14r, T15r, T16r, T17r, T18r, T19r, T20r, T21r, T22r, T23r, T24r, T25r, T26r, T28r, T29r, T30r]

        #valuesu = [T10u, T11u, T12u, T13u, T14u, T15u, T16u, T17u, T18u, T19u, T20u, T21u, T22u, T23u, T24u, T25u, T26u, T27u, T28u, T29u, T30u]

        xlist = []
        ylist = []

        xlistu = []
        ylistu = []
        
        plt.figure(figsize = (10,8))
        for i in valuesr:
            xlist.append(i[1][-1])
            ylist.append(i[-1][-1])
        #for i in valuesu:
        #    xlistu.append(i[1][-1])
        #    ylistu.append(i[-1][-1])
            
        plt.plot(xlist, ylist, '-g', label = "Random", lw = 2)
        #plt.plot(xlistu, ylistu, '--b', label = "Up", lw = 2)
        plt.title('Total Accepted States for 1000000 Monte Carlo Cycles')
        plt.legend(loc = 'best')
        plt.xlabel("Temperatures", fontsize = 16)
        plt.ylabel("Numbers of total Accepted", fontsize = 16)
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


    values = [[E1randE, E24randE], [E1upE, E24upE]]
    jcounter = 0
    icounter = 0
    for i in values:
        plt.figure(figsize= (10,8))
        for j in i:
            if icounter == 0:
                if jcounter == 0:
                    plt.subplot(211)
                    plt.title('T = 1.0')
                    n, bins, patch = plt.hist(j, bins = 20, normed = 1, label = 'Random')
                    plt.ylabel('Probability', fontsize = 16)
                    print np.var(j)/20**2
                    print T1rand[4][-1]

                else:
                    plt.subplot(212)
                    plt.title('T = 2.4')
                    n, bins, patch = plt.hist(j, bins = 40, normed = 1, label = 'Random')
                    plt.ylabel('Probability', fontsize = 16)
                    print np.var(j)/20**2
                    print T24rand[4][-1]
            """
            else:
                if jcounter == 0:
                    plt.subplot(211)
                    n, bins, patch = plt.hist(j, bins = 20, normed = 1, label = 'Up')
                    plt.ylabel('Probability', fontsize = 16)
                else:
                    plt.subplot(212)
                    n, bins, patch = plt.hist(j, bins = 40, normed = 1, label = 'Up')
                    plt.ylabel('Probability', fontsize = 16)
           """
            jcounter += 1
        icounter += 1
            
        plt.xlabel('Energies', fontsize = 16)
        #plt.legend(loc = 'best')
        plt.savefig('../figures/task_d/hist%d.png' %icounter)
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
    #L20  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim20_dirrand_dt0.03333.xyz', 20)
    L40  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim40_dirrand_dt0.01034.xyz', 40)
    L60  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim60_dirrand_dt0.01034.xyz', 60)
    L80  = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim80_dirrand_dt0.01034.xyz', 80)
    L100 = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim100_dirrand_dt0.01034.xyz', 100)
    L140 = reader('../benchmarks/task_e/eigenvalues_MC1.0e+06_dim140_dirrand_dt0.00612.xyz', 140)

    values = [L40, L60, L80, L100, L140]

    f1 = 20
    f2 = 22
    
    #################################
    ####           E             ####
    #################################
    if str(temp) == "E" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Energy eigenvalue', fontsize = f1)

        counter = 40
        for i in values:
            plt.plot(i[1], i[2], label = 'L%d' %counter)

            if counter == 100:
                counter += 40
            else:
                counter += 20
        
        plt.legend(loc = 'best')
        plt.ylabel('<E>', fontsize = f1)
        plt.xlabel('Temperature', fontsize = f1)
        plt.legend(loc='best' )
        plt.savefig("../figures/task_e/energyeig.png")
        plt.show()

    #################################
    ####           M             ####
    #################################
    if str(temp) == "M" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Magnetic eigenvalue', fontsize = f1)

        counter = 40
        for i in values:
            plt.plot(i[1], i[3], label = 'L%d' %counter)

            if counter == 100:
                counter += 40
            else:
                counter += 20
        
        plt.legend(loc = 'best')
        plt.ylabel('<M>', fontsize = f1)
        plt.xlabel('Temperature', fontsize = f1)
        plt.legend(loc='best' )
        plt.savefig("../figures/task_e/Mageig.png")
        plt.show()


    #################################
    ####           sE            ####
    #################################
    if str(temp) == "sE" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Heat Capacity', fontsize = f1)

        counter = 40
        for i in values:
            plt.plot(i[1], i[4]/i[1]**2, label = 'L%d' %counter)

            print counter, i[1][np.argmax(i[4]/i[1])]

            if counter == 100:
                counter += 40
            else:
                counter += 20

        plt.legend(loc = 'best')
        plt.ylabel(r'$C_v$', fontsize = f2)
        plt.xlabel('Temperature', fontsize = f1)
        plt.legend(loc='best')
        plt.savefig("../figures/task_e/sigmaE.png")
        plt.show()

    
    #################################
    ####           sM            ####
    #################################
    if str(temp) == "sM" or str(temp) == "all":
        plt.figure(figsize = (10,8))
        plt.title('Susceptibility', fontsize = f1)

        counter = 40
        for i in values:
            plt.plot(i[1], i[5]/i[1], label = 'L%d' %counter)

            print counter, i[1][np.argmax(i[5]/i[1])]

            if counter == 100:
                counter += 40
            else:
                counter += 20
            print i[1]

                
        plt.legend(loc = 'best')
        plt.ylabel(r'$\chi$', fontsize = f2)
        plt.xlabel('Temperature', fontsize = f1)
        plt.legend(loc='best' )
        plt.savefig("../figures/task_e/suscept.png")
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
