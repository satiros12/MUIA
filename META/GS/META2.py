import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import SA 
import time 
import matplotlib.pylab as pylab
pylab.rcParams['figure.figsize'] = (16.0, 10.0)
import itertools

#Default maximization
class SimulatedAnnealing:
    def __init__(self):
        self.T_initial, self.T_final, self.K_steps, self.T_current, self.step, self.alpha, slef.F_current, self.F_new = [None]*8
        self.Number_steps, self.error_percentage, self.best_F, self.Last_k_Fs, self.accepted_percentage= [None]*5
        self.k, self.acepted_in_last_k_steps, self.cutoff, self.number_cutoff_acepted_solutions = [None]*4
        self.TU = None
        self.P = 0
        self.END = []
    
def linearTimmeUpdate(T_initial, T_final, K_steps, T_current):
    Betta = (T_initial - T_final*1.0) / (K_steps*1.0* T_initial * T_final)
    return T_current*1.0/(1.0 + Betta*T_current)

def geometriclTimmeUpdate(T_initial, step, alpha=0.95):
    #Betta = #(T_initial - T_final*1.0) / (K_steps*1.0* T_initial * T_final)
    return (alpha**step) * T_initial  # T_current*1.0/(1.0 + Betta*T_current)

def probabilityOfSelection(T_current, F_current, F_new):
    A= min([1.0, np.exp(- (F_new*1.0-F_current)*1.0/T_current)]) > np.random.rand() if F_new > F_current else True
    #print A,F_new > F_current,np.exp(- (F_new*1.0-F_current)*1.0/T_current)
    return A

def endCondiction_error_modification(Number_steps, error_percentage, best_F, Last_k_Fs):
    #print (min(Last_k_Fs) - best_F),error_percentage
    return np.abs(min(Last_k_Fs) - best_F)/np.abs(best_F) > error_percentage   #error_percentage*1.0  * Number_steps     

def endCondiction_accepted_percentage(Number_steps, accpeted_percentage, acepted_in_last_k_steps):
    return accpeted_percentage*Number_steps < acepted_in_last_k_steps 

def endCondiction_cutoffs(Number_steps, cutoff, number_aceoted_solutions):
     return cutoff*Number_steps < number_aceoted_solutions 
    
def endCondiction_MaxSteps(Number_steps, step):
         return Number_steps < step 
    
def fo(x): return np.cos(x*1.0)*1.0/x


def SelectTemperaturesRange(Initial_percentage_aceptance, Final_percentage_aceptance, Test_times, Objective_function, Restrictions, Envirnoment, verbose=False):
    Initial_points = np.random.rand(Test_times) * Restrictions[1] + Restrictions[0]
    Perturbations = np.random.rand(Test_times) * Envirnoment*2 - Envirnoment
    New_points = Initial_points + Perturbations
    Correct_perturbed_1 = np.where(New_points >= Restrictions[0])[0]
    Correct_perturbed_2 = np.where(New_points <= Restrictions[1])[0]
    Correct_perturbed = np.intersect1d(Correct_perturbed_1,Correct_perturbed_2)
    #filter(lambda x : x >= Restrictions[0] and x <= Restrictions[1] , 
    #print New_points[Correct_perturbed], len(Correct_perturbed)
    Non_improvements = np.where(Objective_function(Initial_points[Correct_perturbed]) < Objective_function(New_points[Correct_perturbed]))[0] 
    #print Objective_function(Initial_points[Non_improvements]) < Objective_function(New_points[Non_improvements])
    #print len(Non_improvements)
    MEAN_DIF = np.mean( np.abs(Objective_function(Initial_points[Correct_perturbed][Non_improvements]) - Objective_function(New_points[Correct_perturbed][Non_improvements])) )
    Initial_temperature = -MEAN_DIF/ np.log(Initial_percentage_aceptance)
    Final_temperature = -MEAN_DIF/ np.log(Final_percentage_aceptance)
    if verbose:
        RES2 = SelectTemperaturesRange(Initial_percentage_aceptance, Final_percentage_aceptance, Test_times, Objective_function, Restrictions, Envirnoment, verbose=False)
        print "In ",Test_times," tests, with envirnoment : ",Envirnoment
        print "Mean difference of non improvement perturbations : ", MEAN_DIF
        print "Intial temperature for a ", Initial_percentage_aceptance, "% of aceptance : ",Initial_temperature
        print "Intial percentage optained for other Sample :", np.exp(-RES2[2] / Initial_temperature)
        print "Final temperature for a ", Final_percentage_aceptance, "% of aceptance : ",Final_temperature
        print "Final percentage optained for other Sample:", np.exp(-RES2[2] / Final_temperature)
    return Initial_temperature,Final_temperature, MEAN_DIF


#Simulated Annealing Continuous
def SIM_AN(
MAX_ITERATIONS=10000,
InitialPoint = np.random.rand()*30.0,
EnvirnomentRange = 1,
objectiveFunction = fo,
restrictions=[10**-100,30],
Intial_percentage =0.999,
Final_percentage =0.3,
InitialTemperature=100,
FinalTemperature = 1000,
USE_TEMPERATURE=False,
alpha = 0.99,
TEMPERATURE_EVOLUTION=0,
debug=False,
verbose=2,
plotose=False,
#if probabilityFunction == None : probabilityFunction=lambda a,b,T : probabilityOfSelection(T,b,a)
accpeted_percentage=1.0,
error_percentage=10000.0,
k1= 500,
k2=500,
cutoff =1.0,

UP = 0.1,
LL = 0.01,

REC = 0.01,
DES = 0.001,

Ec = 0.01,


Metodo_de_aceptacion = 0
#0 : Probabilidad
#1 : Umbral
#2 : Gran Diluvio
#3 : Recuerdo del recuerdo del viaje
#4 : Microcacanonic annealing method
):
    if not USE_TEMPERATURE : 
        InitialTemperature,FinalTemperature, Mean_RE = SelectTemperaturesRange(Intial_percentage,Final_percentage,100,fo,[10**-4,30],1, verbose=verbose)
    CurrentTempreature = InitialTemperature
    Current_solution = [InitialPoint , objectiveFunction(InitialPoint)]
    Best_Solution = [InitialPoint , objectiveFunction(InitialPoint)]
    if debug:
        OutOfRangeRandoms = 0
    BEST_Acumulative_values = [Best_Solution[-1]]
    Acumulative_values_prefered = [Current_solution[-1]]
    Acumulative_values_tested = [Current_solution[-1]]
    BEST_Acumulative_sol = [Best_Solution[0]]
    Acumulative_sol_prefered = [Current_solution[0]]
    Acumulative_sol_tested = [Current_solution[0]]
    InitialTime = time.time()
    Step = 0
    k_buffer_acceptance=[]
    k_buffer_F=[]
    Total_aceptance_count = 0
    Ecs,LLs,RECs,Ts = [],[],[],[]
    while True:
        Step += 1
        StepTime = time.time()
        Ecs += [Ec]
        LLs += [LL]
        RECs += [REC]
        Ts += [CurrentTempreature]
        while True: #Solo cogemos datos restringidos
            RandomPoint = np.random.rand()*2.0 - 1.0
            New_solution = Current_solution[0] + RandomPoint*EnvirnomentRange 
            #print New_solution
            if restrictions[0] == None or restrictions[0] <= New_solution:
                if restrictions[1] == None or restrictions[1] >= New_solution:
                    break
            if debug : OutOfRangeRandoms+=1
        New_value = objectiveFunction(New_solution)
        if debug: print "Sol ",Step,":",New_solution,New_value, " <-> ",Current_solution[0],Current_solution[1]

        #print New_value, New_solution, Current_solution
        Update_current_solution = False
        RES_TP = [(Metodo_de_aceptacion == 0 and probabilityOfSelection(CurrentTempreature, Current_solution[1], New_value))
                                    ,(Metodo_de_aceptacion == 1 and np.abs(New_value - Current_solution[1]) > CurrentTempreature)
                                    ,(Metodo_de_aceptacion == 2 and (New_value) > LL)
                                    ,(Metodo_de_aceptacion == 3 and (New_value) > (REC-DES))
                                    ,(Metodo_de_aceptacion == 4 and np.abs(New_value - Current_solution[1]) < Ec)]
        #print RES_TP
        Update_current_solution = sum(RES_TP)
        #if Metodo_de_aceptacion == 0 : Update_current_solution = probabilityOfSelection(CurrentTempreature, Current_solution[1], New_value)
        #if Metodo_de_aceptacion == 1 : Update_current_solution = np.abs(New_value - Current_solution[1]) > CurrentTempreature
        #if Metodo_de_aceptacion == 2 : Update_current_solution = (New_value) > LL
        #if Metodo_de_aceptacion == 3 : Update_current_solution = (New_value) > (REC-DES)
        #if Metodo_de_aceptacion == 4 : Update_current_solution = np.abs(New_value - Current_solution[1]) < Ec
        #if Update_current_solution == 0 : print  " +++ SI +++"
        #print Metodo_de_aceptacion,Update_current_solution
        #print probabilityOfSelection(CurrentTempreature, Current_solution[1], New_value),np.abs(New_value - Current_solution[1]) > CurrentTempreature,(New_value) > LL,(New_value) > (REC-DES),np.abs(New_value - Current_solution[1]) < Ec
        if New_value > REC:
            REC = New_value

        if Update_current_solution:
            Ec = Ec - np.abs(New_value - Current_solution[1])
            Current_solution[0] = New_solution
            Current_solution[1] = New_value
            Total_aceptance_count += 1
            LL += UP


        if Best_Solution[1] > Current_solution[1]:
            Best_Solution[0] = Current_solution[0]
            Best_Solution[1] = Current_solution[1]


        #termination conditions
        if Step > k1 : 
                k_buffer_acceptance[:-1] = k_buffer_acceptance[1:]
                k_buffer_acceptance[-1] = 1 if New_value == Current_solution[1] else 0
        else:
            k_buffer_acceptance.append(1 if New_value == Current_solution[1] else 0)


        if Step > k2 :
                k_buffer_F[:-1] = k_buffer_F[1:]
                k_buffer_F[-1] = New_value   
        else: 
            k_buffer_F.append(New_value)

        BEST_Acumulative_values += [Best_Solution[1]]
        Acumulative_values_prefered += [Current_solution[1]]
        Acumulative_values_tested += [New_value]
        BEST_Acumulative_sol += [Best_Solution[0]]
        Acumulative_sol_prefered += [Current_solution[0]]
        Acumulative_sol_tested += [New_solution]

        #print "Current Temp",CurrentTempreature
        if TEMPERATURE_EVOLUTION == 0:
            CurrentTempreature = linearTimmeUpdate(InitialTemperature, FinalTemperature, MAX_ITERATIONS, CurrentTempreature)
        if TEMPERATURE_EVOLUTION == 1:
            CurrentTempreature = geometriclTimmeUpdate(InitialTemperature, Step,alpha)
        
        if plotose :
            Acumulative_values.append(Best_Solution[-1])
        if verbose == 1:  
            print "Step : ",Step," Duration:",(time.time() - StepTime) , "Temperature:",CurrentTempreature," Solution F: ",Current_solution[1]," Best solution F: ",Best_Solution[1]
        if endCondiction_MaxSteps(MAX_ITERATIONS, Step):
            if verbose : print "End by Max Steps"
            break
        if Step > k1 :
            if endCondiction_accepted_percentage(MAX_ITERATIONS, accpeted_percentage, sum(k_buffer_acceptance)):
                if verbose : print "End by Accepted Percentage",MAX_ITERATIONS
                break
        if Step > k2 :
            if endCondiction_error_modification(MAX_ITERATIONS, error_percentage, Best_Solution[1], k_buffer_F):
                if verbose : print "End by Error Modification"
                break
        if endCondiction_cutoffs(MAX_ITERATIONS, cutoff, Total_aceptance_count):
            if verbose : print "End by CUTOFFS"
            break

    if plotose :
            plt.plot(Acumulative_values)
            plt.title("Best" + str(Best_Solution[1]))
            plt.show()
    if verbose > 0:
        TotalTime = (time.time()-InitialTime)
        print "Total time:" ,TotalTime, " Time per step", TotalTime*1.0/MAX_ITERATIONS,"Steps:",Step, "Best solution : ",Best_Solution[1]
    OTHER_OUTPUT= []
    OTHER_OUTPUT += [[BEST_Acumulative_sol,BEST_Acumulative_values]]
    OTHER_OUTPUT += [[Acumulative_sol_prefered, Acumulative_values_prefered]]
    OTHER_OUTPUT += [[Acumulative_sol_tested, Acumulative_values_tested]]
    OTHER_OUTPUT += [Ecs,LLs,RECs,Ts]
    return Best_Solution,Step, OTHER_OUTPUT


def plot_results(ALGO_RESULT,SHOW_IMG=True, SAVE_FILE=False, File_name="GS"):
    plt.plot(RES[2][0][1],"r")
    plt.plot(RES[2][1][1],"b")
    plt.plot(RES[2][2][1],"g")
    plt.xlim(0,len(RES[2][0][1]))
    plt.title("Evolution of : F(x) = cos(x)/x")
    plt.xlabel("Number of steps")
    plt.ylabel("f")
    if SAVE_FILE :plt.savefig(File_name+'_EVF.png')
    if SHOW_IMG: plt.show()
    plt.close()

    plt.plot(RES[2][0][0],"r")
    plt.plot(RES[2][1][0],"b")
    plt.plot(RES[2][2][0],"g")
    plt.xlim(0,len(RES[2][0][1]))
    plt.title("Evolution of : x")
    plt.xlabel("Number of steps")
    plt.ylabel("x")
    if SAVE_FILE :plt.savefig(File_name+'_EVX.png')
    if SHOW_IMG: plt.show()
    plt.close()

    plt.plot(RES[2][-1],"r")
    plt.xlim(0,len(RES[2][0][1]))
    #plt.hlines(2.89255277615,0,Step)
    plt.title("Evolution of : Temperature")
    plt.xlabel("Number of steps")
    plt.ylabel("Temperature")
    if SAVE_FILE :plt.savefig(File_name+'_EVT.png')
    if SHOW_IMG: plt.show()
    plt.close()

    plt.plot(RES[2][-2],"r")
    plt.xlim(0,len(RES[2][0][1]))
    plt.title("Evolution of : Records")
    plt.xlabel("Number of steps")
    plt.ylabel("Records")
    if SAVE_FILE :plt.savefig(File_name+'_EVREC.png')
    if SHOW_IMG: plt.show()
    plt.close()

    plt.plot(RES[2][-3],"r")
    plt.xlim(0,len(RES[2][0][1]))
    plt.title("Evolution of : Rain")
    plt.xlabel("Number of steps")
    plt.ylabel("Rain")
    if SAVE_FILE :plt.savefig(File_name+'_EVLL.png')
    if SHOW_IMG: plt.show()
    plt.close()

    plt.plot(RES[2][-4],"r")
    plt.xlim(0,len(RES[2][0][1]))
    plt.title("Evolution of : Ec")
    plt.xlabel("Number of steps")
    plt.ylabel("Ec")
    if SAVE_FILE :plt.savefig(File_name+'_EVEc.png')
    if SHOW_IMG: plt.show()
    plt.close()