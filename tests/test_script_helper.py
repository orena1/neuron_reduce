'''Test script for subtree reductor:
 May be used to test the NeuroReduce reduction as performed by
 subtree_reductor() in subtree_reductor_func.py.  Receives as arguments a
 Neurolucida full morphology file, a model template file, the total approximate
 number of segments in the reduced model (-1 for default calculation of a
 segment for every 0.1 lambda), and a file with synapse locations (see format
 in Readme) - creates the full morphology instance, reduces it with
 subtree_reductor(), and compares the voltage traces from the simulation of the
 original full model and the reduced model using the same random synaptic
 activity in both simulations, which is determined by parameters in this file.

 How to use: Usage from shell: python test_script.py original_morphology_file
 model_template_file reduction_frequency manual_total_nsegs
 synapse_locations_file for example: python test_script.py
 2013_03_06_cell08_876_H41_05_Cell2.ASC model.hoc 38 -1
 origRandomSynapses-10000

 - This tester supports only Neurolucida morphology files. To test NeuroReduce
   with other kinds of morphology files, CODE CHANGES are required in
   create_model()
 - The model template file given to this tester must contain a
   complete_full_model_creation() function (see Readme and example in the
                                            attached model.hoc file) to ensure
   the correct and complete creation of the full model.
 - This tester supports excitatory and inhibitory synapses - randomly splits
   the locations in the single given synapse locations file into excitatory and
   inhibitory synapses according to the distribution determined by the user in
   PERCENTAGE_OF_EXCITATORY_SYNAPSES; gives each synapse's NetStim
   object a seed that corresponds to the synapse's index as excitatory or
   inhibitory.  To support more types of synapses, CODE CHANGES are required
   in create_synapses()

 - Global simulation and stimuli parameters can be changed
   according to the user's needs.
'''
from __future__ import print_function

from collections import namedtuple
from contextlib import contextmanager
import sys
import random
import os
import logging
import subprocess
from neuron import h
import numpy as np
import matplotlib.pyplot as plt

from neuron_reduce import subtree_reductor

logging.basicConfig(level=os.environ.get("LOGLEVEL", "DEBUG"))

#############################
# simulation parameters - to be adjusted by the user
SIM_PARAMS = {'tstop': 1500,  # ms
              'steps_per_ms': 10,
              'dt': 0.1,
              }

SPIKE_BEGIN_TIME = 10  # ms
WEIGHT = 0.0008  # gbar = 0.8 nS

ICLAMP_PARAMS = {'duration': 1500,  # ms
                 'delay': 10,
                 'amplitude': 0,  # nA
                 }

EXCITATORY_STIMULI_INTERVAL = 1000  # ms   -  100 ms interval is 10 hz
EXCITATORY_STIMULI_NUMBER = 15
FREQUENCY_OF_INHIBITORY_SYNAPSES_IS_HIGHER_BY = 10

INHIBITORY_STIMULI_INTERVAL = (float(EXCITATORY_STIMULI_INTERVAL) /
                               FREQUENCY_OF_INHIBITORY_SYNAPSES_IS_HIGHER_BY)  # ms
INHIBITORY_STIMULI_NUMBER = int(float(ICLAMP_PARAMS['duration']) / INHIBITORY_STIMULI_INTERVAL)

SYNAPSE_PARAMS = {'excitatory': {'e_syn': 0,  # mV
                                 'tau_1': 0.3,  # ms (like AMPA)
                                 'tau_2': 1.8,  # ms (like AMPA)
                                 'stimuli_interval': EXCITATORY_STIMULI_INTERVAL,
                                 'stimuli_number': EXCITATORY_STIMULI_NUMBER,
                                 },
                  'inhibitory': {'e_syn': -86,
                                 'tau_1': 1.,  # ms
                                 'tau_2': 8.,  # ms
                                 'stimuli_interval': INHIBITORY_STIMULI_INTERVAL,
                                 'stimuli_number': INHIBITORY_STIMULI_NUMBER,
                                 },
                  }

PERCENTAGE_OF_EXCITATORY_SYNAPSES = 85.
STIMULI_NOISE = 1
SEED_FOR_RANDOM_SYNAPSE_CLASSIFICATION = 1


def abs_diff(a, b):
    '''return summed absolute difference'''
    return np.sum(np.abs(a - b))


def rmsd(reduced, control):
    '''find root mean squared deviation'''
    diff = reduced - control
    return np.sqrt(np.dot(diff, diff) / reduced.size)


@contextmanager
def chdir(path):
    '''change to path and change back'''
    cwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(cwd)


def plot_traces(time, reduction_frequency,
                original_control_voltage, reduced_cell_voltage,
                np_orig_recording_vec_reduced, np_recording_vec_reduced,
                rmsd_new, rmsd_old, rmsd_two_reduced_vecs,
                np_recording_vec_control, np_orig_recording_vec_control):
    '''plot traces'''
    plt.figure(1)
    plt.plot(time, original_control_voltage, label='original model ', c='b')
    plt.plot(time, reduced_cell_voltage, label='subtree reduction ', c='g')
    plt.plot(time, np_orig_recording_vec_reduced, label='orginal subtree reduction ', c='r')

    plt.xlabel('time (ms)')
    plt.ylabel('voltage at soma (mV)')

    plt.title("new Reduction RMSD %s old %s )" % (rmsd_new, rmsd_old))
    plt.legend(loc='lower right')

    plt.figure(2)
    plt.plot(time, np_recording_vec_control, label='new control - voltage at soma (mV)', c='b')
    plt.plot(time, np_orig_recording_vec_control,
             label='original control - voltage at soma (mV)', c='r')
    plt.title("differences in control -  Hz (sum = %s)" %
              abs_diff(np_recording_vec_control, np_orig_recording_vec_control))
    plt.legend(loc='lower right')

    plt.figure(3)
    plt.plot(time, np_orig_recording_vec_reduced,
             label='original reduction - voltage at soma (mV)', c='b')
    plt.plot(time, np_recording_vec_reduced,
             label='new reduction - voltage at soma (mV)', c='g')
    plt.title("differences in Reduction - %s Hz (sum = %s)" %
              (reduction_frequency, rmsd_two_reduced_vecs))
    plt.legend(loc='lower right')
    plt.show()


def create_model(morphology_file, model_obj_name, instance_name, create_type):
    '''Creates a full morphology cell instance from Neurolucida morphology model template filenames

    This function is called to create the original instance to be reduced and
    to create the control cell

    TODO: To support non-Neurolucida morphology files, this function should be CHANGED.
    '''
    assert instance_name in ('original_cell', 'control_cell'), \
        "name must be one of ('original_cell', 'control_cell')"

    h("{objref %s}" % instance_name)
    model = dict(instance_name=instance_name, model_obj_name=model_obj_name)
    if create_type in ('basic', 'human'):
        h("{instance_name} = new {model_obj_name}()".format(**model))

        # instantiates according to morphology using import3d
        nl = h.Import3d_Neurolucida3()
        nl.quiet = 1
        nl.input(morphology_file)
        import_3d = h.Import3d_GUI(nl, 0)

        instance_created = getattr(h, instance_name)
        # associates cell instance with the morphology file
        import_3d.instantiate(instance_created)
        # this function should be included in the model.hoc file given as a parameter
        instance_created.complete_full_model_creation()
    elif create_type in ('bbp', 'bbpactive'):
        h('{instance_name} = new {model_obj_name}(1, "{morphology_file}")'.format(
            morphology_file=morphology_file, **model))
        h('{instance_name} = {instance_name}.CellRef'.format(**model))

    elif create_type in ('bbpnew', ):
        with chdir(morphology_file[:morphology_file.rindex('/')]):
            h('{instance_name} = new {model_obj_name}(0)'.format(**model))
    elif create_type in ('hay', 'almog', 'allen'):
        h('{instance_name} = new {model_obj_name}("{morphology_file}")'.format(
            morphology_file=morphology_file, **model))

    return getattr(h, instance_name)


SynapseLocation = namedtuple('SynapseLocation', 'type, section_num, x')


def create_synapses(cell_instance,
                    synapse_file,
                    seed_for_random_synapse_classification=SEED_FOR_RANDOM_SYNAPSE_CLASSIFICATION,
                    percentage_of_excitatory_synapses=PERCENTAGE_OF_EXCITATORY_SYNAPSES):
    '''Creates the synapses according to the synapse locations file

    attaches them to the given cell instance.
    Returns the new synapses_list, netstims_list, netcons_list, and
    randoms_list (the activity generators for the NetStims).

    TODO: Supports two types of synapses, to support more types CODE CHANGES are required
    '''
    synapses = []
    with open(synapse_file) as fd:
        for line in fd:
            line = line.strip()
            if not line:
                continue

            # extracting the synapse location
            line_content = line.strip("[,]").split(",")
            type_ = line_content[0].strip()
            section_num = int(line_content[1].strip())
            x = float(line_content[2].strip())  # synapse's relative location in the section
            synapses.append(SynapseLocation(type_, section_num, x))

    synapses_list, netstims_list, netcons_list, randoms_list = [], [], [], []
    num_of_inhibitory_syns = num_of_excitatory_syns = 0

    # the seed used to randomly classify synapses from the file into different types
    rand = random.Random(seed_for_random_synapse_classification).random

    # iterates over synapse locations file and creates synapses, NetStims and NetCons
    for i, synapse in enumerate(synapses):
        if synapse.type == 'apical':
            cell_instance.apic[synapse.section_num].push()
        elif synapse.type == 'somatic':
            cell_instance.soma[0].push()
        else:  # for synapses on basal subtrees
            cell_instance.dend[synapse.section_num].push()

        synapse = h.Exp2Syn(synapse.x)

        h.pop_section()

        # This code should be changed in order to enable more synaptic types:
        # randomly classifies the synapse as excitatory or inhibitory
        # (according to the Random object instantiated above)
        # classifies according to the user-decided distribution of synapses
        if rand() <= percentage_of_excitatory_synapses / 100.:  # excitatory
            params = SYNAPSE_PARAMS['excitatory']
            # to appoint its NetStim's random-activity generating Random object
            # a seed according to its index as an excitatory synapse
            randoms_list.append(h.Random(num_of_excitatory_syns))
            num_of_excitatory_syns += 1
        else:  # inhibitory synapse
            params = SYNAPSE_PARAMS['inhibitory']
            # to appoint its NetStim's random-activity generating Random object
            # a seed according to its index as an inhibitory synapse
            randoms_list.append(h.Random(num_of_inhibitory_syns))
            num_of_inhibitory_syns += 1

        synapse.e = params['e_syn']
        synapse.tau1 = params['tau_1']
        synapse.tau2 = params['tau_2']
        stimuli_interval = params['stimuli_interval']
        stimuli_number = params['stimuli_number']
        synapses_list.append(synapse)

        # creates a NetStim for the synapse and assigns it properties
        netstim = h.NetStim()
        netstim.interval = stimuli_interval  # ms
        netstim.number = stimuli_number
        netstim.start = SPIKE_BEGIN_TIME
        netstim.noise = STIMULI_NOISE

        randoms_list[i].negexp(1)  # must specify negexp distribution with mean = 1
        # assigns the NetStim the Random object created above - to determine its activity
        netstim.noiseFromRandom(randoms_list[i])
        netstims_list.append(netstim)

        # creates a NetCon and connects it to the synapse and NetStim
        netcons_list.append(h.NetCon(netstims_list[i], synapses_list[i]))
        netcons_list[i].weight[0] = WEIGHT  # the synaptic weight
        netcons_list[i].delay = 0

    return synapses_list, netstims_list, netcons_list, randoms_list


def generate_random_synapses_locations(cell, syn_filename, num_syns=10000, seed=2):
    basal_L = 0
    tot_L = 0
    syn_L = []
    random_obj = random.Random(seed)

    for sec in cell.basal:
        basal_L += sec.L
        tot_L += sec.L

    for sec in cell.apical:
        tot_L += sec.L

    for _ in range(num_syns):
        L_counter = 0
        x_L = random_obj.uniform(0, tot_L)
        if x_L < basal_L:
            for ix, sec in enumerate(list(cell.basal)):
                if L_counter + sec.L > x_L:
                    x = (x_L - L_counter) / sec.L
                    syn_L.append(["basal", ix, x])
                    break
                else:
                    L_counter += sec.L
        else:
            x_L -= basal_L
            for ix, sec in enumerate(list(cell.apical)):
                if L_counter + sec.L > x_L:
                    x = (x_L - L_counter) / sec.L
                    syn_L.append(["apical", ix, x])
                    break
                else:
                    L_counter += sec.L

    print(syn_filename)
    with open(syn_filename, 'w') as f:
        for syn in syn_L:
            f.write("[%s,%s,%s]\n" % (syn[0], syn[1], syn[2]))


def loadtemplate(template_file):
    '''load template if not already defined'''
    model_obj_name = template_file.split(".hoc")[0].split('/')[-1]
    if h.name_declared(model_obj_name) == 0:
        print("loading template '%s'" % model_obj_name)
        h.load_file(template_file)
    else:
        print("WARNING The template '%s' is already defined... not loading." % model_obj_name)


def change_cell_location(soma, offset=150):
    soma.push()

    for i in range(int(h.n3d())):
        h.pt3dchange(i,
                     offset + h.x3d(i),
                     offset + h.y3d(i),
                     offset + h.z3d(i),
                     h.diam3d(i))


def write_test_vectors(voltage_file,
                       reduction_frequency, synapse_file,
                       recording_vec_control, recording_vec_reduced):
    '''
    The two vectors are printed consecutively without extra spaces, first the
    control vector and after it the reduced vector. Each voltage value of the
    vector is in a new line.
    '''

    with open(voltage_file, 'w') as fd:
        fd.write("This file represents the control and reduced voltage vectors "
                 "for the following parameters: "
                 "frequency = {REDUCTION_FREQUENCY} Hz, "
                 "synapse file = {SYNAPSE_FILE} , "
                 "EXCITATORY_E_SYN = {EXCITATORY_E_SYN} mV, "
                 "EXCITATORY_TAU_1 = {EXCITATORY_TAU_1} ms, "
                 "EXCITATORY_TAU_2 = {EXCITATORY_TAU_2} ms, "
                 "EXCITATORY_STIMULI_INTERVAL = {EXCITATORY_STIMULI_INTERVAL} ms "
                 "EXCITATORY_STIMULI_NUMBER = {EXCITATORY_STIMULI_NUMBER}, "
                 "STIMULI_NOISE = {STIMULI_NOISE}, "
                 "weight = {WEIGHT} microSiemens, "
                 "ICLAMP_AMPLITUDE = {ICLAMP_AMPLITUDE} nA, "
                 "ICLAMP_DELAY = {ICLAMP_DELAY} ms, "
                 "ICLAMP_DURATION = {ICLAMP_DURATION} ms, "
                 "tstop = {TSTOP} ms, "
                 "dt = {DT}, "
                 "SPIKE_BEGIN_TIME = {SPIKE_BEGIN_TIME}ms\n".format(
                     REDUCTION_FREQUENCY=reduction_frequency,
                     SYNAPSE_FILE=synapse_file,
                     EXCITATORY_E_SYN=SYNAPSE_PARAMS['excitatory']['e_syn'],
                     EXCITATORY_TAU_1=SYNAPSE_PARAMS['excitatory']['tau_1'],
                     EXCITATORY_TAU_2=SYNAPSE_PARAMS['excitatory']['tau_2'],
                     EXCITATORY_STIMULI_INTERVAL=EXCITATORY_STIMULI_INTERVAL,
                     EXCITATORY_STIMULI_NUMBER=EXCITATORY_STIMULI_NUMBER,
                     STIMULI_NOISE=STIMULI_NOISE,
                     WEIGHT=WEIGHT,
                     ICLAMP_AMPLITUDE=ICLAMP_PARAMS['amplitude'],
                     ICLAMP_DELAY=ICLAMP_PARAMS['delay'],
                     ICLAMP_DURATION=ICLAMP_PARAMS['duration'],
                     TSTOP=SIM_PARAMS['tstop'],
                     DT=SIM_PARAMS['dt'],
                     SPIKE_BEGIN_TIME=SPIKE_BEGIN_TIME,
                 ))
        for i in range(int(recording_vec_control.size())):
            fd.write(repr(recording_vec_control.x[i]) + "\n")
        for i in range(int(recording_vec_reduced.size())):
            fd.write(repr(recording_vec_reduced.x[i]) + "\n")


def load_voltages(voltage_file, number_of_recordings):
    '''load_voltages'''
    with open(voltage_file) as fd:
        fd.readline()   # skips first line with the simulation configuration details
        voltages = [float(v) for v in fd]

    assert len(voltages) == number_of_recordings * 2
    np_orig_recording_vec_control = np.array(voltages[:number_of_recordings])
    np_orig_recording_vec_reduced = np.array(voltages[number_of_recordings:])
    return np_orig_recording_vec_control, np_orig_recording_vec_reduced


def setup_iclamps(sec, params):
    '''sets iclamp for `sec` with `params`'''
    iclamp_reduced_cell = h.IClamp(0.5, sec=sec)
    iclamp_reduced_cell.delay = params['delay']
    iclamp_reduced_cell.dur = params['duration']
    iclamp_reduced_cell.amp = params['amplitude']
    return iclamp_reduced_cell


def run_test(morphology_file,
             model_file,
             model_file_reduced,
             reduction_frequency,
             manual_total_nsegs,
             synapse_file,
             voltage_file,
             write_unit_test_vectors,
             plot_voltages,
             create_type,
             celsius):
    '''Run a test reduction

    Note: This must be run from a new python instance so that neuron is reset
    '''
    h.celsius = celsius  # needs to be set for reduction to work properly

    reduction_frequency = int(reduction_frequency)
    print('---------------sim data --------------------------')
    print("REDUCTION_FREQUENCY %s" % reduction_frequency)
    print("MODEL_FILE          %s" % model_file)
    print("MODEL_FILE_REDUCED  %s" % model_file_reduced)
    print("SYNAPSE_FILE        %s" % synapse_file)
    print("MORPHOLOGY_FILE     %s" % morphology_file)
    print("MANUAL_TOTAL_NSEGS  %s" % manual_total_nsegs)
    print("VOLTAGEFILE         %s" % voltage_file)
    print('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n\n')

    if create_type in ('almog', 'hay', 'bbpnew', 'bbpactive', 'allen', 'human'):
        with chdir(model_file[:model_file.rindex('/')]):
            process = subprocess.Popen(['nrnivmodl', 'mod'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            print('compiling mod files')
            stdout, stderr = process.communicate()
            #print(stdout)
            #print(stderr)
            h.nrn_load_dll("x86_64/.libs/libnrnmech.so.0")
            h.nrn_load_dll("x86_64/.libs/libnrnmech.0.so")
            h.nrn_load_dll("x86_64/.libs/libnrnmech.so")
            

    if create_type == 'hay':
        loadtemplate(model_file[:model_file.rindex('/')] + '/L5PCbiophys3.hoc')
    elif create_type == 'allen':
        loadtemplate(model_file[:model_file.rindex('/')] + '/AllenBiophys.hoc')

    h.load_file("import3d.hoc")
    loadtemplate(model_file)

    # creates original cell with synapses, and reduces it with neuron_reduce
    model_obj_name = os.path.basename(model_file.split(".hoc")[0])
    original_cell = create_model(morphology_file, model_obj_name, "original_cell", create_type)

    # creates control cell with synapses
    control_cell = create_model(morphology_file, model_obj_name, "control_cell", create_type)
    synapses_list_control, netstims_list_control, netcons_list_control, random_list_control = \
        create_synapses(control_cell, synapse_file)

    # simulates control and reduced cell instances together
    synapses_list, netstims_list, netcons_list, random_list = \
        create_synapses(original_cell, synapse_file)
    reduced_cell, synapses_list, netcons_list = subtree_reductor(original_cell,
                                                                 synapses_list,
                                                                 netcons_list,
                                                                 reduction_frequency,
                                                                 model_file_reduced,
                                                                 manual_total_nsegs)

    if control_cell.soma.hname()[-1] == ']':
        control_soma = control_cell.soma[0]
        reduced_soma = reduced_cell.soma[0]
    else:
        control_soma = control_cell.soma
        reduced_soma = reduced_cell.soma

    h.steps_per_ms = SIM_PARAMS['steps_per_ms']
    h.dt = SIM_PARAMS['dt']
    h.tstop = SIM_PARAMS['tstop']
    h.v_init = reduced_soma.e_pas

    # set iclamps
    setup_iclamps(reduced_soma, ICLAMP_PARAMS)
    setup_iclamps(control_soma, ICLAMP_PARAMS)

    # sets recording vectors
    recording_vec_reduced = h.Vector()
    recording_vec_control = h.Vector()

    recording_vec_reduced.record(reduced_soma(.5)._ref_v)
    recording_vec_control.record(control_soma(.5)._ref_v)

    print('Running simulations, temperature is ' + str(celsius) + 'c')
    h.run()

    # for debugging, it helps if the two cells don't overlap
    change_cell_location(control_soma)

    np_recording_vec_control = np.array(recording_vec_control)
    np_recording_vec_reduced = np.array(recording_vec_reduced)

    if write_unit_test_vectors:
        write_test_vectors(voltage_file,
                           reduction_frequency, synapse_file,
                           recording_vec_control, recording_vec_reduced)

    number_of_recordings = int(recording_vec_control.size())
    np_orig_recording_vec_control, np_orig_recording_vec_reduced = \
        load_voltages(voltage_file, number_of_recordings)

    print('\n\n--------------- unit test results --------------')
    control_vecs_equal = np.allclose(np_recording_vec_control,
                                     np_orig_recording_vec_control,
                                     atol=1e-8)
    if not control_vecs_equal:
        print("UNIT TEST FAILED: control voltage vector has been changed")
        print('Complex model: unit test V - current sim V : %s' %
              abs_diff(np_recording_vec_control, np_orig_recording_vec_control))

    reduced_vecs_equal = np.allclose(np_recording_vec_reduced,
                                     np_orig_recording_vec_reduced,
                                     atol=1e-8)
    if not reduced_vecs_equal:
        print("UNIT TEST FAILED: reduced voltage vector has been changed")
        print('Reduced model: unit test V - current sim V : %s' %
              abs_diff(np_recording_vec_reduced, np_orig_recording_vec_reduced))

    if control_vecs_equal and reduced_vecs_equal:
        print("UNIT TEST PASSED: no significant changes to voltage vectors\n")
        print('Complex model: unit test V - current sim V : %s' %
              abs_diff(np_recording_vec_control, np_orig_recording_vec_control))
        print('Reduced model: unit test V - current sim V : %s' %
              abs_diff(np_recording_vec_reduced, np_orig_recording_vec_reduced))

    # plotting graphs of both voltage traces
    dt = SIM_PARAMS['dt']
    time = np.arange(dt, dt * number_of_recordings + dt, dt)
    rmsd_new = rmsd(np_recording_vec_reduced, np_recording_vec_control)
    rmsd_old = rmsd(np_orig_recording_vec_reduced, np_orig_recording_vec_control)
    rmsd_two_reduced_vecs = rmsd(np_recording_vec_reduced, np_orig_recording_vec_reduced)
    print("current sim: rmsd between reduced vs complex model:", rmsd_new)
    print("unit test: rmsd between reduced vs complex model:", rmsd_old)
    print("rmsd between current sim reduced and unit test reduced :", rmsd_two_reduced_vecs)

    if plot_voltages:
        plot_traces(time, reduction_frequency,
                    np_recording_vec_control, np_recording_vec_reduced,
                    np_orig_recording_vec_reduced, np_recording_vec_reduced,
                    rmsd_new, rmsd_old, rmsd_two_reduced_vecs,
                    np_recording_vec_control, np_orig_recording_vec_control)
    print('--------------------------------------- END -----------------------------------')
    return control_vecs_equal and reduced_vecs_equal


def str_to_bool(s):
    s = s.lower()
    assert s in ('true', 'false'), '%s Must be a either true of false' % s
    return s == 'true'


def main(argv):
    '''main'''
    orig_morphology_file = argv[1]
    orig_model_file = argv[2]
    reduced_model_file = argv[3]
    frequency = float(argv[4])
    manual_total_nsegs = int(argv[5])
    synapse_file = argv[6]
    voltage_file = argv[7]
    write_unit_test_vectors = str_to_bool(argv[8])
    plot_voltages = str_to_bool(argv[9])
    create_type = argv[10]
    celsius = float(argv[11])
    return run_test(orig_morphology_file,
                    orig_model_file,
                    reduced_model_file,
                    frequency,
                    manual_total_nsegs,
                    synapse_file,
                    voltage_file,
                    write_unit_test_vectors,
                    plot_voltages,
                    create_type,
                    celsius)


if __name__ == "__main__":
    '''Run unit tests'''
    print(sys.argv)
    sys.exit(not main(sys.argv))
