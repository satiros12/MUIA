import numpy as np
import matplotlib.pyplot as plt
import time 

def gradualTimmeUpdate(T_initial, T_final, K_steps, T_current):
    Betta = (T_initial - T_final*1.0) / (K_steps*1.0* T_initial * T_final)
    return T_current*1.0/(1.0 + Betta*T_current)

#Simulated Annealing Continuous
def SAC(MAX_ITERATIONS, InitialPoint, EnvirnomentRange, InitialTemperature, objectiveFunction, timeUpdateFunction, restrictions=[None,None]
        ,probabilityFunction=None, randomFunction=np.random.rand, debug=False, verbose=False, plotose=False ):
    if probabilityFunction == None : probabilityFunction=lambda a,b,T : np.exp(-(objectiveFunction(a) - objectiveFunction(b))/T)
    CurrentTempreature = InitialTemperature
    Current_solution = [InitialPoint , objectiveFunction(InitialPoint)]
    Best_Solution = [InitialPoint , objectiveFunction(InitialPoint)]
    if plotose :
        Acumulative_values = [Best_Solution[-1]]
    InitialTime = time.time()
    for i in xrange(MAX_ITERATIONS):
        StepTime = time.time()
        while True:
            RandomPoint = (randomFunction(len(Current_solution[0])) - 0.5) * 2
            New_solution = Current_solution[0] + RandomPoint*EnvirnomentRange 
            #print New_solution
            if restrictions[0] == None or restrictions[0] <= New_solution:
                if restrictions[1] == None or restrictions[1] >= New_solution:
                    break
        New_value = objectiveFunction(New_solution)
        if debug: print "Sol ",i,":",New_solution,New_value, " <-> ",Current_solution[0],Current_solution[1]
        if New_value < Current_solution[1] : 
            Current_solution[0] = New_solution
            Current_solution[1] = New_value
        else:
            ExpectedProbility = min([1,probabilityFunction(Current_solution[1],New_value,CurrentTempreature)])
            ObtainedProbability = np.random.rand()
            if debug : print "Exp : ",ExpectedProbility,ObtainedProbability
            if ObtainedProbability <= ExpectedProbility:
                if Best_Solution[1] > Current_solution[1]:
                    Best_Solution[0] = Current_solution[0]
                    Best_Solution[1] = Current_solution[1]
                Current_solution[0] = New_solution
                Current_solution[1] = New_value
        
        #print "Current Temp",CurrentTempreature
        CurrentTempreature = timeUpdateFunction(i,CurrentTempreature)
        if plotose :
            Acumulative_values.append(Best_Solution[-1])
        if verbose == 1:  
            print "Step : ",i," Duration:",(time.time() - StepTime) , "Temperature:",CurrentTempreature," Solution: ",Current_solution[1]," Best solution: ",Best_Solution[1]
    if plotose :
            plt.plot(Acumulative_values)
            plt.title("Best" + str(Best_Solution[1]))
            plt.show()
    if verbose > 0:
        TotalTime = (time.time()-InitialTime)
        print "Total time:" ,TotalTime, " Time per step", TotalTime*1.0/MAX_ITERATIONS, "Best solution : ",Best_Solution[1]
    return Best_Solution

def SAC_GTU(MAX_ITERATIONS, InitialPoint, EnvirnomentRange,  InitialTemperature, FinalTemperature,objectiveFunction, restrictions=[None,None], debug=False, verbose=False, plotose=False):
    return SAC(MAX_ITERATIONS, InitialPoint, EnvirnomentRange, InitialTemperature, objectiveFunction=objectiveFunction,restrictions=restrictions, timeUpdateFunction=lambda i,T : gradualTimmeUpdate(InitialTemperature,FinalTemperature,MAX_ITERATIONS,T), debug=debug, verbose=verbose, plotose=plotose)


def test1():
    def fo(x): return np.cos(x)/x
    np.random.seed(1)
    R = SAC_GTU(1000,np.asanyarray([10.0]),5.0,100.0,1.0,fo, restrictions=[0,30], plotose=False,verbose=False)
    print R
    print "Test result : ",(int(R[1][0]*10000) == -3364) and (int(R[0][0]*10000) == 28166)