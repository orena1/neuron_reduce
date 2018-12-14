# Test script for subtree reductor:
#  May be used to test the NeuroReduce reduction as performed by subtree_reductor() in subtree_reductor_func.py.
#  Receives as arguments a Neurolucida full morphology file, a model template file, the total approximate number of segments in the reduced model (-1 for default calculation of a segment for every 0.1 lambda), 
#  and a file with synapse locations (see format in Readme) - creates the full morphology instance, reduces it with subtree_reductor(),
#  and compares the voltage traces from the simulation of the original full model and the reduced model using the same random synaptic activity in both simulations, which is determined by parameters in this file.
#
# How to use:
#  Usage from shell: python test_script.py original_morphology_file model_template_file reduction_frequency manual_total_nsegs synapse_locations_file
#  for example: python test_script.py 2013_03_06_cell08_876_H41_05_Cell2.ASC model.hoc 38 -1 origRandomSynapses-10000
#
# - This tester supports only Neurolucida morphology files. To test NeuroReduce with other kinds of morphology files, CODE CHANGES are required in createModel() - lines 79-110.
# - The model template file given to this tester must contain a complete_full_model_creation() function (see Readme and example in the attached model.hoc file) to ensure the correct and complete creation of the full model.
# - This tester supports excitatory and inhibitory synapses - randomly splits the locations in the single given synapse locations file into excitatory and inhibitory synapses according to the distribution 
#    determined by the user in line 51 (PERCENTAGE_OF_EXCITATORY_SYNAPSES); gives each synapse's NetStim object a seed that corresponds to the synapse's index as excitatory or inhibitory.
#    To support more types of synapses, CODE CHANGES are required (lines 157-181).
# - Global simulation and stimuli parameters (line 31-61) can be changed according to the user's needs.



from neuron import h, gui

import matplotlib.pyplot as plt
import math
import sys
import random
import numpy as np
import pdb
import os
import subprocess

import neuron_reduce
subtree_reductor =  neuron_reduce.subtree_reductor

APICAL = "apical"
SOMATIC = "somatic"

#############################
# simulation parameters - to be adjusted by the user:
TSTOP = 1500  # ms
STEPS_PER_MS = 10
DT = 0.1

SPIKE_BEGIN_TIME = 10  # ms
WEIGHT = 0.0008  # gbar = 0.8 nS

ICLAMP_DURATION = 1500 # ms
ICLAMP_DELAY = 10
ICLAMP_AMPLITUDE = 0 #nA

EXCITATORY_E_SYN = 0 #mV
EXCITATORY_TAU_1 = 0.3 #ms  (like AMPA)
EXCITATORY_TAU_2 = 1.8 #ms  (like AMPA)

INHIBITORY_E_SYN = -86
INHIBITORY_TAU_1 = 1 #ms
INHIBITORY_TAU_2 = 8 #ms

PERCENTAGE_OF_EXCITATORY_SYNAPSES = 85               
FREQUENCY_OF_INHIBITORY_SYNAPSES_IS_HIGHER_BY = 10

EXCITATORY_STIMULI_INTERVAL = 1000  #ms   -  100 ms interval is 10 hz
EXCITATORY_STIMULI_NUMBER = 15           

INHIBITORY_STIMULI_INTERVAL = float(EXCITATORY_STIMULI_INTERVAL) / FREQUENCY_OF_INHIBITORY_SYNAPSES_IS_HIGHER_BY  # ms  
INHIBITORY_STIMULI_NUMBER = int(float(ICLAMP_DURATION) / INHIBITORY_STIMULI_INTERVAL)   

STIMULI_NOISE = 1
#############################

REDUCTION_FREQUENCY = None
MODEL_FILE = None
SYNAPSE_FILE = None
MORPHOLOGY_FILE = None

SEED_FOR_RANDOM_SYNAPSE_CLASSIFICATION = 1

MANUAL_TOTAL_NSEGS = -1

def createModel(morphology_file, model_obj_name, purpose, create_type):
    '''
    Creates an full morphology cell instance from the given Neurolucida morphology file and given model template filename. 
    This function is called to create the original instance to be reduced and to create the control cell (the given argument purpose equals "reduction" or "control" accordingly).
    To support non-Neurolucida morphology files, this function should be CHANGED.
    '''
    if purpose == "reduction":
        instance_name = "original_cell"
    else:
        instance_name = "control_cell"
    
    if create_type in ['basic' , 'human']:
        # creates instance
        h(instance_name + " = new " + model_obj_name + "()")
        
        # instantiates according to morphology using import3d
        nl = h.Import3d_Neurolucida3()
        nl.quiet = 1
        nl.input(morphology_file)
        imprt = h.Import3d_GUI(nl, 0)

        
        instance_created = h.original_cell if purpose == "reduction" else h.control_cell
        # associates cell instance with the morphology file
        imprt.instantiate(instance_created)    
        instance_created.complete_full_model_creation() # this function should be included in the model.hoc file given as a parameter
        
    elif create_type in ['bbp', 'bbpactive']:
        h(instance_name + " = new " + model_obj_name + '(1,"' + morphology_file +'")')
        h(instance_name + '=' + instance_name + '.CellRef')
        instance_created = h.original_cell if purpose == "reduction" else h.control_cell
    elif create_type == 'bbpnew':
        cwd = os.getcwd()
        os.chdir(morphology_file[:morphology_file.rindex('/')])
        h(instance_name + " = new " + model_obj_name + '(0)')
        os.chdir(cwd)
        instance_created = h.original_cell if purpose == "reduction" else h.control_cell

    elif create_type in ['hay', 'almog','allen']:
        h(instance_name + " = new " + model_obj_name + '("' + morphology_file +'")')
        instance_created = h.original_cell if purpose == "reduction" else h.control_cell


    return instance_created 




def create_synapses(cell_instance):
    ''' 
    Creates the synapses according to the synapse locations file given to the program and attaces them to the given cell instance. 
    Returns the new synapses_list, netstims_list, netcons_list, and randoms_list (the activity generators for the NetStims).
    Supports two types of synapses, to support more types CODE CHANGES are required (lines 156-180).
    '''
    f = open(SYNAPSE_FILE)
    synapses_list = []
    netstims_list = []
    netcons_list = []
    randoms_list = []
    
    random_classifier = random.Random(SEED_FOR_RANDOM_SYNAPSE_CLASSIFICATION)  # the seed used to randomly classify synapses from the file into different types
    
    num_of_excitatory_syns = 0
    num_of_inhibitory_syns = 0
    
    # iterates over synapse locations file and creates synapses, NetStims and NetCons
    orig_synapse_index = 0
    for line in iter(f):
        line = line.strip()
        if line:  # for non-empty lines
            # extracting the synapse location
            line = line.strip("[,]")
            line_content = line.split(",")
            section_num = int(line_content[1].strip())
            x =  float(line_content[2].strip())  # synapse's relative location in the section

            # adding the synapse to its section in the model according to the extracted location
            if line_content[0].strip() == APICAL:
                cell_instance.apic[section_num].push()
            elif line_content[0].strip() == SOMATIC:
                cell_instance.soma[0].push()
            else:  # for synapses on basal subtrees
                cell_instance.dend[section_num].push()

            synapses_list.append(h.Exp2Syn(x))

            h.pop_section()
            
            ##########################################################################################################################
            # This code should be changed in order to enable more synaptic types:
            
            # randomly classifies the synapse as excitatory or inhibitory (according to the Random object instantiated above)
            excitatory_syn = False
            if random_classifier.random() <= float(PERCENTAGE_OF_EXCITATORY_SYNAPSES) / 100:  # classifies according to the user-decided distribution of synapses
                excitatory_syn = True
            
            # decides on the synapse's properties according to its type
            if excitatory_syn:
                e_syn = EXCITATORY_E_SYN
                tau1 = EXCITATORY_TAU_1
                tau2 = EXCITATORY_TAU_2
                stimuli_interval = EXCITATORY_STIMULI_INTERVAL
                stimuli_number = EXCITATORY_STIMULI_NUMBER
                num_of_excitatory_syns = num_of_excitatory_syns + 1
                randoms_list.append(h.Random(num_of_excitatory_syns-1))  # to appoint its NetStim's random-activity generating Random object a seed according to its index as an excitatory synapse

            else: # inhibitory synapse
                e_syn = INHIBITORY_E_SYN
                tau1 = INHIBITORY_TAU_1
                tau2 = INHIBITORY_TAU_2
                stimuli_interval = INHIBITORY_STIMULI_INTERVAL
                stimuli_number = INHIBITORY_STIMULI_NUMBER
                num_of_inhibitory_syns = num_of_inhibitory_syns + 1
                randoms_list.append(h.Random(num_of_inhibitory_syns-1))  # to appoint its NetStim's random-activity generating Random object a seed according to its index as an inhibitory synapse
            ##########################################################################################################################
            
            # assigns the properties to the synapse
            synapses_list[orig_synapse_index].e = e_syn
            synapses_list[orig_synapse_index].tau1 = tau1
            synapses_list[orig_synapse_index].tau2 = tau2
            
            # creates a NetStim for the synapse and assigns it properties
            netstims_list.append(h.NetStim())
            netstims_list[orig_synapse_index].interval = stimuli_interval #ms
            netstims_list[orig_synapse_index].number = stimuli_number
            netstims_list[orig_synapse_index].start = SPIKE_BEGIN_TIME
            netstims_list[orig_synapse_index].noise = STIMULI_NOISE
            randoms_list[orig_synapse_index].negexp(1)  # must specify negexp distribution with mean = 1
            netstims_list[orig_synapse_index].noiseFromRandom(randoms_list[orig_synapse_index]) # assigns the NetStim the Random object created above - to determine its activity           
            
            # creates a NetCon and connects it to the synapse and NetStim
            netcons_list.append(h.NetCon(netstims_list[orig_synapse_index], synapses_list[orig_synapse_index]))
            netcons_list[orig_synapse_index].weight[0] = WEIGHT  # the synaptic weight
            netcons_list[orig_synapse_index].delay = 0

            orig_synapse_index = orig_synapse_index + 1
    
    f.close()
    return synapses_list, netstims_list, netcons_list, randoms_list
  
def generate_random_synapses_locations(cell,syn_filename,num_syns=10000,seed=2):
    basal_L = 0
    tot_L = 0
    syn_L = []
    random_obj = random.Random(seed)

    for sec in cell.basal:
        basal_L += sec.L
        tot_L += sec.L
    for sec in cell.apical:
        tot_L += sec.L
    print "basal_L,tot_L",basal_L,tot_L
    for i in range(num_syns):
        L_counter = 0
        x_L = random_obj.uniform(0,tot_L)
        if x_L<basal_L:
            for ix,sec in enumerate(list(cell.basal)):
                if L_counter+sec.L>x_L:
                    x = (x_L - L_counter)/sec.L
                    syn_L.append(["basal",ix,x])
                    break
                else:
                    L_counter+=sec.L
        else:
            x_L -=basal_L
            for ix,sec in enumerate(list(cell.apical)):
                if L_counter+sec.L>x_L:
                    x = (x_L - L_counter)/sec.L
                    syn_L.append(["apical",ix,x])
                    break
                else:
                    L_counter+=sec.L

    print syn_filename
    with open(syn_filename,'w') as f:
        for syn in syn_L:
            f.write("[%s,%s,%s]"%(syn[0],syn[1],syn[2])+"\n")

def set_up_parameters(orig_morphology_file, model_file, frequency, manual_total_nsegs, synapse_file):
    ''' Sets the program's arguments as global parameters for the entire testing process '''



def loadtemplate(template_file):
    ''' Check if the template is already defined, if it does not defined it load it, if it is defined it does not reload the template'''
    
    model_obj_name = template_file.split(".")[0].split('/')[-1]
    if h.name_declared(model_obj_name)==0:
        print("loading template '" + model_obj_name +"'")
        h.load_file(template_file)
    else:
        print("WARNING The template '" +  model_obj_name + "' is already defined... not loading.")


def change_cell_location(control_cell, soma_list):
    if soma_list:
        control_cell.soma[0].push() 
    else:
        control_cell.soma.push()
    
    for i in range(0,int(h.n3d())):
        h.pt3dchange(i, 150+h.x3d(i), 150+h.y3d(i), 150+h.z3d(i), h.diam3d(i))


def RunTest(orig_morphology_file, orig_model_file, reduced_model_file, frequency, manual_total_nsegs, synapse_file, voltage_file, 
            write_unit_test_vectors=False, plot_voltages=False, create_type='basic',celsius=34):
    '''
    Run unit test
    variables : = ?
    '''
    global REDUCTION_FREQUENCY, MODEL_FILE, SYNAPSE_FILE, MORPHOLOGY_FILE, MANUAL_TOTAL_NSEGS
    
    #First Delete all sections
    h("forall delete_section()")
    h.celsius = celsius
    
    #Set variables
    REDUCTION_FREQUENCY = int(frequency)
    MODEL_FILE = orig_model_file
    MODEL_FILE_REDUCED = reduced_model_file
    SYNAPSE_FILE = synapse_file
    MORPHOLOGY_FILE = orig_morphology_file
    MANUAL_TOTAL_NSEGS = manual_total_nsegs
    VOLTAGEFILE = voltage_file
    print('---------------sim data --------------------------')
    print("REDUCTION_FREQUENCY " + str(REDUCTION_FREQUENCY))
    print("MODEL_FILE          " + str(MODEL_FILE))
    print("MODEL_FILE_REDUCED  " + str(MODEL_FILE_REDUCED))
    print("SYNAPSE_FILE        " + str(SYNAPSE_FILE))
    print("MORPHOLOGY_FILE     "  + str(MORPHOLOGY_FILE))
    print("MANUAL_TOTAL_NSEGS  " + str(MANUAL_TOTAL_NSEGS))
    print("VOLTAGEFILE         "  + str(VOLTAGEFILE))
    print('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n\n')

    if create_type in ['almog','hay','bbpnew', 'bbpactive','allen','human']:
        cwd = os.getcwd()
        os.chdir(MODEL_FILE[:MODEL_FILE.rindex('/')])
        process = subprocess.Popen(['nrnivmodl','mod'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('compiling mod files')
        stdout, stderr = process.communicate()
        #print(stdout)
        #print(stderr)
        h.nrn_load_dll("x86_64/.libs/libnrnmech.so.0")
        h.nrn_load_dll("x86_64/.libs/libnrnmech.0.so")
        os.chdir(cwd)

    if create_type == 'hay':
        loadtemplate(MODEL_FILE[:MODEL_FILE.rindex('/')]+'/L5PCbiophys3.hoc')
    if create_type =='allen':
        loadtemplate(MODEL_FILE[:MODEL_FILE.rindex('/')]+'/AllenBiophys.hoc')
    h.load_file("import3d.hoc")
    loadtemplate(MODEL_FILE)
    

    # creates original cell with synapses, and reduces it with NeuroReduce
    h("{objref original_cell}")
    
    model_obj_name = MODEL_FILE.split(".")[0].split('/')[-1]
    original_cell = createModel(MORPHOLOGY_FILE, model_obj_name, "reduction", create_type)
    
    # generate_random_synapses_locations(original_cell,SYNAPSE_FILE)


   

    
    # creates control cell with synapses
    h("{objref control_cell}")
    control_cell = createModel(MORPHOLOGY_FILE, model_obj_name, "control", create_type)
    synapses_list_control, netstims_list_control, netcons_list_control, random_list_control = create_synapses(control_cell)

    # simulates control and reduced cell instances together
    synapses_list, netstims_list, netcons_list, random_list = create_synapses(original_cell)
    reduced_cell, synapses_list, netcons_list = subtree_reductor(original_cell, synapses_list, netcons_list, REDUCTION_FREQUENCY, MODEL_FILE_REDUCED, MANUAL_TOTAL_NSEGS)
    
    h.steps_per_ms = STEPS_PER_MS
    h.dt = DT
    
    # check soma type
    soma_list = True if control_cell.soma.hname()[-1]==']' else False
    
   
    # sets iclamps
    iclamp_reduced_cell = h.IClamp(0.5, sec = reduced_cell.soma[0] if soma_list else  reduced_cell.soma)
    iclamp_reduced_cell.delay = ICLAMP_DELAY
    iclamp_reduced_cell.dur = ICLAMP_DURATION
    iclamp_reduced_cell.amp = ICLAMP_AMPLITUDE
    
    iclamp_control_cell = h.IClamp(0.5, sec = control_cell.soma[0] if soma_list else control_cell.soma)
    iclamp_control_cell.delay = ICLAMP_DELAY
    iclamp_control_cell.dur = ICLAMP_DURATION
    iclamp_control_cell.amp = ICLAMP_AMPLITUDE
    
    # sets recording vectors
    recording_vec_reduced = h.Vector()  # used to record the voltage throughout the simulation
    recording_vec_control = h.Vector()  # used to record the voltage throughout the simulation
    if soma_list:
        recording_vec_reduced.record(reduced_cell.soma[0](.5)._ref_v)
        recording_vec_control.record(control_cell.soma[0](.5)._ref_v)
    else:
        recording_vec_reduced.record(reduced_cell.soma(.5)._ref_v)
        recording_vec_control.record(control_cell.soma(.5)._ref_v)
    
    
    h.tstop = TSTOP 
    h.v_init = reduced_cell.soma[0].e_pas if soma_list else  reduced_cell.soma.e_pas
    h.celsius = celsius
    print('Running simulations, temperature is ' + str(celsius) + 'c')
    # running the simulation
    h.run()
    
    change_cell_location(control_cell, soma_list)

    np_recording_vec_control = np.array(recording_vec_control)
    np_recording_vec_reduced = np.array(recording_vec_reduced)
  
    if write_unit_test_vectors:
        voltage_vector_file = open(VOLTAGEFILE, 'w')
        voltage_vector_file.write("This file represents the control and reduced voltage vectors for the following parameters " +
                                  ": frequency = " + str(REDUCTION_FREQUENCY) + " Hz, synapse file = " + SYNAPSE_FILE  + 
                                   ", EXCITATORY_E_SYN = " + str(EXCITATORY_E_SYN) + " mV, EXCITATORY_TAU_1 = " +str(EXCITATORY_TAU_1) + " ms, " +
                                   " EXCITATORY_TAU_2 = " + str(EXCITATORY_TAU_2)  + " ms, EXCITATORY_STIMULI_INTERVAL = " + str(EXCITATORY_STIMULI_INTERVAL) + " ms (" + str(EXCITATORY_STIMULI_INTERVAL/1000.0) + "hz)," +
                                   " EXCITATORY_STIMULI_NUMBER = " + str(EXCITATORY_STIMULI_NUMBER) + " , STIMULI_NOISE = " + str(STIMULI_NOISE) + ", weight =" +  str(WEIGHT) +" microSiemens," +
                                   " ICLAMP_AMPLITUDE = " + str(ICLAMP_AMPLITUDE)  + "nA, ICLAMP_DELAY = " + str(ICLAMP_DELAY) + " ms, ICLAMP_DURATION = " +str(ICLAMP_DURATION) +" ms," +
                                   " tstop = " + str(TSTOP) + " ms, dt = " + str(DT)+ ", SPIKE_BEGIN_TIME = " + str(SPIKE_BEGIN_TIME)+ "ms\n")
        for i in xrange(int(recording_vec_control.size())):
            voltage_vector_file.write(repr(recording_vec_control.x[i])+ "\n")
        for i in xrange(int(recording_vec_reduced.size())):
            voltage_vector_file.write(repr(recording_vec_reduced.x[i])+ "\n")
        voltage_vector_file.close()
    # The two vectors are printed consecutively without extra spaces, first the control vector and after it the reduced vector. Each voltage value of the vector is in a new line.

    
    # unit test - compares voltage vectors of this run with the correct voltage vectors for a certain simulation/synapses configuration (detailed at the top of voltage_vectors_for_unit_test.txt):
    orig_recording_vec_control = []  # the original control voltage vector that is printed in voltage_vectors_for_unit_test
    orig_recording_vec_reduced = []  # the original reduced voltage vector that is printed in voltage_vectors_for_unit_test
    voltage_vector_file = open(VOLTAGEFILE, 'r')
    voltage_vector_file.readline()   # skips first line with the simulation configuration details
    for i in range(int(recording_vec_control.size())):
        orig_recording_vec_control.append(float(voltage_vector_file.readline())) 
    for i in range(int(recording_vec_reduced.size())):
        orig_recording_vec_reduced.append(float(voltage_vector_file.readline())) 
    voltage_vector_file.close()

    np_orig_recording_vec_reduced = np.array(orig_recording_vec_reduced)
    np_orig_recording_vec_control = np.array(orig_recording_vec_control)
    
    control_vecs_equal = True
    for i in range(int(recording_vec_control.size())):
        if abs(recording_vec_control.x[i] - orig_recording_vec_control[i])>1e-8:
            control_vecs_equal = False
            #print(i, recording_vec_control.x[i], orig_recording_vec_control[i])
        else:
            continue
    print('\n\n--------------- unit test results --------------')
    if control_vecs_equal == False:
        print("UNIT TEST FAILED: control voltage vector has been changed")
        print('Complex model: unit test V - current sim V :' + str(sum(abs(np.array(recording_vec_control)-np.array(orig_recording_vec_control)))))
        
    reduced_vecs_equal = True
    for i in range(int(recording_vec_reduced.size())):
        if abs(recording_vec_reduced.x[i] - orig_recording_vec_reduced[i])>1e-8:
            reduced_vecs_equal = False
            break
        else:
            continue

    if reduced_vecs_equal == False:
        print("UNIT TEST FAILED: reduced voltage vector has been changed")
        print('Reduced model: unit test V - current sim V :' + str(sum(abs(np.array(recording_vec_reduced)-np.array(orig_recording_vec_reduced)))))
        
    if control_vecs_equal and reduced_vecs_equal:
        print("UNIT TEST PASSED: no significant changes to voltage vectors\n")
        print('Complex model: unit test V - current sim V :' + str(sum(abs(np.array(recording_vec_control)-np.array(orig_recording_vec_control)))))
        print('Reduced model: unit test V - current sim V :' + str(sum(abs(np.array(recording_vec_reduced)-np.array(orig_recording_vec_reduced)))))
    #############################
     
     
    # calculates rmsd (root mean square deviation) between the voltage trace of the full model and that of the reduced model ( rmsd = sqrt(E[(x2 - x1)^2])  )
    number_of_recordings = int(recording_vec_reduced.size())
    sqrd_err_tot = 0
    for j in range(number_of_recordings):
        curr_sqrd_err = (recording_vec_reduced.x[j]- recording_vec_control.x[j])**(2)
        sqrd_err_tot = sqrd_err_tot + curr_sqrd_err
    rmsd = math.sqrt(float(sqrd_err_tot)/ number_of_recordings)
     
     
    # plotting graphs of both voltage traces
    time = []
    original_control_voltage = []
    reduced_cell_voltage = []
    time = np.arange(DT,DT*number_of_recordings+DT,DT)
    original_control_voltage = np.array(recording_vec_control)
    reduced_cell_voltage = np.array(recording_vec_reduced)
    rmsd_new = np.sqrt(np.sum((np_recording_vec_reduced-np_recording_vec_control)**2)/np_recording_vec_reduced.size)
    rmsd_old = np.sqrt(np.sum((np_orig_recording_vec_reduced-np_orig_recording_vec_control)**2)/np_orig_recording_vec_reduced.size)
    rmsd_two_reduced_vecs = np.sqrt(np.sum((np_recording_vec_reduced-np_orig_recording_vec_reduced)**2)/np_orig_recording_vec_reduced.size)
    print "current sim: rmsd between reduced vs complex model:", rmsd_new
    print "unit test: rmsd between reduced vs complex model:", rmsd_old
    print "rmsd between current sim reduced and unit test reduced :", rmsd_two_reduced_vecs
    if plot_voltages:
        plt.figure(1)
        plt.plot(time, original_control_voltage, label = 'original model ',c='b')
        plt.plot(time, reduced_cell_voltage, label = 'subtree reduction ',c='g')
        plt.plot(time,np_orig_recording_vec_reduced,label = 'orgiinal subtree reduction ',c='r')

        plt.xlabel('time (ms)')
        plt.ylabel('voltage at soma (mV)')

        # plt.title("Subtree Reduction - " + str(REDUCTION_FREQUENCY) + " Hz (rmsd = " + str(rmsd) + ")")
        

        plt.title("new Reduction RMSD "+str(rmsd_new) + " old rmsd " + str(rmsd_old) +")")
        plt.legend(loc = 'lower right')
        plt.figure(2)
        plt.plot(time, np_recording_vec_control, label = 'new control - voltage at soma (mV)',c='b')
        plt.plot(time, np_orig_recording_vec_control   , label = 'original control - voltage at soma (mV)',c='r')
        plt.title("differences in control -  Hz (sum = " + str(np.sum(np.abs(np_recording_vec_control-np_orig_recording_vec_control))) + ")")
        plt.legend(loc = 'lower right')
        plt.figure(3)
        plt.plot(time, np_orig_recording_vec_reduced   , label = 'original reduction - voltage at soma (mV)',c='b')
        plt.plot(time, np_recording_vec_reduced, label = 'new reduction - voltage at soma (mV)',c='g')
        plt.title("differences in Reduction - " + str(REDUCTION_FREQUENCY) + " Hz (sum = " + str(rmsd_two_reduced_vecs)+")")
        plt.legend(loc = 'lower right')
        plt.show()
    h("forall delete_section()")
    print('--------------------------------------- END -----------------------------------')
    if control_vecs_equal and reduced_vecs_equal:
        return(True)
    else:
        return(False)


def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError # evil ValueError that doesn't tell you what the wrong value was
     

if __name__ == "__main__":
    '''
    Run unit tests
    '''
    
    
    print sys.argv
    
    Path = 'TestsFiles/Test_1/'
    orig_morphology_file = sys.argv[1]
    orig_model_file = sys.argv[2]
    reduced_model_file = sys.argv[3]
    frequency = float(sys.argv[4])
    manual_total_nsegs = int(sys.argv[5])
    synapse_file = sys.argv[6]
    voltage_file = sys.argv[7]
    write_unit_test_vectors = str_to_bool(sys.argv[8])
    plot_voltages = str_to_bool(sys.argv[9])
    create_type = sys.argv[10]
    celsius = float(sys.argv[11])
    RunTest(orig_morphology_file, orig_model_file, reduced_model_file, frequency, manual_total_nsegs, synapse_file, voltage_file, 
            write_unit_test_vectors=write_unit_test_vectors, plot_voltages=plot_voltages, create_type=create_type, celsius=celsius)
