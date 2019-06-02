# The main file:
# Contains function subtree_reductor() - which reduces a morphologically detailed cell instance into a morphologically simplified cell instance, according to NeuroReduce and merges synapses of the same type (same reverse potential, tau1, and tau2) that are mapped to the same segment. (see more in Readme on tool and usage)
#
# Function usage: subtree_reductor(original_cell instance, synapses_list, netcons_list, model_template_filename, reduction_frequency, total_segments_manual = -1)  - For details, see comments in function
# Function outputs a reduced cell instance, a new synapses_list, and the netcons_list, which now corresponds to the new synapses.
#
# - The model template file must have an init() function (see example in the attached model.hoc file) and the following public definitions specifying sections and section lists accordingly: 
#    public soma, dend, apic ; public all, somatic, apical, basal
# - Supports numerous types of synapses (two synapses are considered to be of different types if they are different from each other in at least one of the following values: reverse potential, tau1, tau2)


from neuron import h, gui
from .reducing_methods import reduce_subtree, set_up_for_reduction, reduce_synapse, measure_input_impedance_of_subtree
import math
import numpy as np
import cmath
import time
import re

SOMA_LABEL = "soma"

    
def create_sections_in_hoc(type_of_section, num, instance_as_str):
    ''' creates sections in the hoc world according to the given section type and number of sections - in the instance whose name is given as a string '''
    h.string = type_of_section
    h('{sprint(string, "create %s[%d]", string, ' + str(num) +') }')
    h("{execute(string, " + instance_as_str+ ")}" )    

    
def connect_to(child_section, parent_section, loc_on_parent_to_connect, loc_on_child_to_connect):
    ''' connects the given child section in the given loc_on_child_to_connect to the given parent section in the given loc_on_parent_to_connect '''
    child_section.connect(parent_section, loc_on_parent_to_connect, loc_on_child_to_connect)


def append_to_section_lists(section, type_of_sectionlist, instance_as_str):
    ''' appends given section to the sectionlist of the given type and to the "all" sectionlist - in the hoc world in the instance whose name is given as a string '''
    h.string = section + " " + type_of_sectionlist + ".append()"
    h("{execute(string, " + instance_as_str+ ")}" )    
    h.string = section + " all.append()"
    h("{execute(string, " + instance_as_str+ ")}" )    


def find_section_number(section):
    ''' extracts and returns the section number from the given section object '''
    sec_name = h.secname(sec = section)
    ints_in_name = re.findall('\d+', sec_name)
    sec_num = ints_in_name[len(ints_in_name) - 1]  # extracts section number
    return sec_num
    

def calculate_nsegs_from_manual_arg(new_dends_electrotonic_length,dends_nsegs, total_segments_wanted):
    '''
    Calculates the number of segments for each section in the reduced model, according to the given total_segments_wanted and the given new_dends_electrotonic_length (the electrotonic lengths of all the new sections).
    Stores calculations in the given list dends_nsegs.
    Called when the user chooses to give to the program the approximate total number of segments that the reduced model should have (non-default calculation).
    '''
    sum_of_lengths = 0   # total electrotonic length of reduced dendritic cables
    total_segments_in_dendrites = total_segments_wanted - 1 # minus one for the one segment of the soma
    for i in range(len(new_dends_electrotonic_length)):
        sum_of_lengths = sum_of_lengths + new_dends_electrotonic_length[i]
    # the num of segments assigned to each section is in proportion to the section's relative contribution to the total electrotonic length in the model
    for i in range(len(dends_nsegs)):
        new_nseg_to_put = int( (float(new_dends_electrotonic_length[i]) / sum_of_lengths) * total_segments_in_dendrites )
        if new_nseg_to_put < 1:
            new_nseg_to_put = 1
        dends_nsegs[i] = new_nseg_to_put
    
    
def calculate_nsegs_from_lambda(new_cables_dimensions, dends_nsegs):
    '''
    The default calculation for the number of segments for each section in the reduced model, according to the length (in microns) and space constant (= lambda - in microns) that were previously calculated for each section and are given in subtree_dimensions.
    According to this calculation, a segment is formed for every 0.1 * lambda in a section. (lambda = space constant = electrotonic length unit).
    '''  
    for i in range(len(dends_nsegs)):
        length = new_cables_dimensions[i][0]
        space_const_in_microns = new_cables_dimensions[i][2]
        dends_nsegs[i] = int((float(length)/space_const_in_microns)*10/2)*2+1  # for every unit of electrotonic length (length/space_constant such units), ~10 segments are formed

SectionsToDelete = [] # A list of all section that will be deleted in the end, it is easier with a global.
section_per_subtree_index = {}
def mark_subtree_sections_with_subtree_index(root_sec_of_subtree, mapping_sections_to_subtree_index, section_type, subtree_index):
    '''
    Recursively marks all sections in the subtree of the given root section as belonging to the given subtree_index using the given dict mapping_sections_to_subtree_index, as follows:
    mapping_sections_to_subtree_index[(<section_type>, <section_number>)] = given subtree_index
    '''
    SectionsToDelete.append(root_sec_of_subtree)
    section_per_subtree_index.setdefault(subtree_index,[])
    section_per_subtree_index[subtree_index].append(root_sec_of_subtree)
    
    section_number = find_section_number(root_sec_of_subtree)
    if len(root_sec_of_subtree.children()) == 0:
        mapping_sections_to_subtree_index[(section_type, section_number)] = subtree_index
    else:
        for i in range(len(root_sec_of_subtree.children())):
            child = root_sec_of_subtree.children()[i]
            mark_subtree_sections_with_subtree_index(child, mapping_sections_to_subtree_index, section_type, subtree_index)
        mapping_sections_to_subtree_index[(section_type, section_number)] = subtree_index
  

def find_synapse_loc(synapse_or_segment, mapping_sections_to_subtree_index, segment=False):
    ''' Returns the location (<dend_index> , <section index>, <relative location in section: x>) of the given synapse object
        if segment == True, it finds a segment location and not a synapse location
    '''
    if segment==False:
        x = synapse_or_segment.get_loc()  # extracts x: the relative location of the synapse in the section (0<=x<=1)
    elif segment==True:
        synapse_or_segment.sec.push()
        x = synapse_or_segment.x

    # extracts the section type ("soma", "apic", "dend") and the section number out of the section name
    full_sec_name = h.secname()
    sec_name_as_list = full_sec_name.split(".")
    short_sec_name = sec_name_as_list[len(sec_name_as_list)-1]
    section_type = short_sec_name.split("[")[0]
    section_num = re.findall('\d+', short_sec_name)[0]
  
    # finds the index of the subtree that this synapse belongs to using the given mapping_sections_to_subtree_index which maps sections to the subtree indexes that they belong to
    if section_type == "apic":
        subtree_index = mapping_sections_to_subtree_index[("apic", section_num)]
    elif section_type == "dend":
        subtree_index = mapping_sections_to_subtree_index[("basal", section_num)]
    else: # somatic synapse
        h.pop_section()
        return SOMA_LABEL, 0, 0

    h.pop_section()
    
    return subtree_index, int(section_num), x


def find_and_diconnect_Axon(soma_ref):
    '''Searching for an axon, it can be a child of the soma or a parent of the soma.'''
    AxonSection = []
    AxonParent = False
    soma_childrens = [soma_ref.child[i] for i in range(int(soma_ref.nchild()))]
    for sec in soma_childrens:
        if 'axon' in sec.hname() or 'Axon' in sec.hname() or 'hill' in sec.hname() :
            AxonSection.append(sec)
            #disconnect axon
            sec.push()
            h.disconnect()
            h.define_shape()

    if soma_ref.has_parent():
        if 'axon' in soma_ref.parent().hname() or 'Axon' in soma_ref.parent().hname() or 'hill' in soma_ref.parent().hname():
            AxonSection.append(soma_ref.parent())
            AxonParent = True
            soma_ref.push()
            h.disconnect()
        else:
            raise Exception('Soma has a parent which is not an axon')
    if len(AxonSection)>1:
        raise Exception('Soma has a two axons')
    return(AxonSection, AxonParent)




def create_segments_to_mech_dic(remove_mechs=True, exclude=['pas', 'na_ion', 'k_ion', 'ca_ion','h_ion']):
    '''This function copy the create a mapping between a dictionary and the mechanisms that it have
       plus the values of those mechanisms. It also remove the mechanisms from the model in order to
       create a passive model
       
       Arguments:
                remove_mechs - False|True
                               if True remove the mechs after creating the mapping, False - keep the mechs
                               
                exclude - List of all the mechs name that should not be removed
       '''

    segment_to_mech_vals ={}
    mech_names = []
    for sec in SectionsToDelete:
        for seg in sec:
            segment_to_mech_vals[seg] = {}
            for mech in seg:
                segment_to_mech_vals[seg][mech.name()] = {}
                for n in dir(mech):
                    if '__' not in n and n not in ['next','name', 'is_ion']:
                        if not n.endswith('_' + mech.name()) and not mech.name().endswith('_ion'):
                            n += '_' + mech.name()
                        tm = getattr(seg,n)
                        #exec('tm = seg.' + n)
                        segment_to_mech_vals[seg][mech.name()][n] = tm
                        mech_names.append(mech.name())
    
    if remove_mechs==True: # Remove all the mechs from the sections
        for sec in SectionsToDelete:
            sec.push()
            for mech in set(mech_names):
                if mech not in exclude:
                    h("uninsert " + mech)
            h.pop_section()
    return(segment_to_mech_vals)


def create_seg_to_seg(original_cell, roots_of_subtrees, mapping_sections_to_subtree_index, new_cables_nsegs, reduced_cables_electrotonic_length,  
                      reduced_cables_electrotonic_length_as_complex_nums, has_apical, apic, basals, subtree_ind_to_q, mapping_type):
    '''create mapping between segments in the original model to segments in the reduced model
       the only significant argumnet is mapping_type
       if mapping_type == impedance the mapping will be a response to the transfer impedance of each segment to the soma (like the synapses)
       if mapping_type == distance  the mapping will be a response to the distance of each segment to the soma (like the synapses) NOT IMPLEMENTED YET
       
       '''
    
    if mapping_type == 'impedance':
        original_seg_to_reduced_seg = {} # the keys are the segments of the original model, the values are the segments of the reduced model.
        reduced_seg_to_original_seg = {}
        for sub_tree_ind in section_per_subtree_index:
            for sec in section_per_subtree_index[sub_tree_ind]:
                for seg in sec:
                    curr_subtree_index, curr_section_num, x = find_synapse_loc(seg, mapping_sections_to_subtree_index, segment=True)

                    imp_obj, subtree_input_impedance = measure_input_impedance_of_subtree(roots_of_subtrees[curr_subtree_index])
                    if has_apical and curr_subtree_index == 0: # if synapse is on the apical subtree
                        on_basal_subtree = False
                    else:
                        on_basal_subtree = True
                    
                    mid_of_segment_loc = reduce_synapse(original_cell, curr_section_num, x, on_basal_subtree, new_cables_nsegs[curr_subtree_index], imp_obj, subtree_input_impedance,
                                                        reduced_cables_electrotonic_length[curr_subtree_index], reduced_cables_electrotonic_length_as_complex_nums[curr_subtree_index],
                                                        subtree_ind_to_q[sub_tree_ind])
                    reduced_cable_index = curr_subtree_index
                    if on_basal_subtree:
                        if has_apical:
                            new_section_for_synapse = basals[reduced_cable_index-1]
                        else:
                            new_section_for_synapse = basals[reduced_cable_index]
                    else:
                        new_section_for_synapse = apic

                    reduced_seg = new_section_for_synapse(mid_of_segment_loc)
                    original_seg_to_reduced_seg[seg] = reduced_seg
                    if reduced_seg not in reduced_seg_to_original_seg: reduced_seg_to_original_seg[reduced_seg] =[]
                    reduced_seg_to_original_seg[reduced_seg].append(seg)
    elif mapping_type == 'distance':
        raise Exception('Not impleMENTED yet!')
    return(original_seg_to_reduced_seg, reduced_seg_to_original_seg)


def copy_dendritic_mech(original_cell, roots_of_subtrees, mapping_sections_to_subtree_index, new_cables_nsegs, reduced_cables_electrotonic_length,
                      reduced_cables_electrotonic_length_as_complex_nums, has_apical, apic, basals, segment_to_mech_vals, subtree_ind_to_q, mapping_type='impedance'):
    ''' copies the mechanisms from the original model to the reduced model'''
    
    #create segment to segment mapping
    original_seg_to_reduced_seg, reduced_seg_to_original_seg = create_seg_to_seg(original_cell, roots_of_subtrees, mapping_sections_to_subtree_index, new_cables_nsegs, reduced_cables_electrotonic_length,  
                                                                                 reduced_cables_electrotonic_length_as_complex_nums, has_apical, apic, basals, subtree_ind_to_q, mapping_type)
    
    #copy mechanisms
    mech_names_per_segment = {} # this is needed for the case where some segements were not been mapped
    vals_per_mech_per_segment = {}
    for reduced_seg in reduced_seg_to_original_seg:
        vals_per_mech_per_segment[reduced_seg] = {}
        mech_names_per_segment[reduced_seg] = []

        vals_per_mech = {}
        mech_names = []
        for original_seg in reduced_seg_to_original_seg[reduced_seg]:
            for mech in segment_to_mech_vals[original_seg]:
                for n in segment_to_mech_vals[original_seg][mech]:
                    if n not in vals_per_mech_per_segment[reduced_seg]: vals_per_mech_per_segment[reduced_seg][n] = []
                    vals_per_mech_per_segment[reduced_seg][n].append(segment_to_mech_vals[original_seg][mech][n]) 

                mech_names_per_segment[reduced_seg].append(mech)
        for mech in mech_names_per_segment[reduced_seg]:
            reduced_seg.sec.insert(mech)
        for n in vals_per_mech_per_segment[reduced_seg]:
            exec('reduced_seg.' + n +' = np.mean(vals_per_mech_per_segment[reduced_seg][n])')
    
    
    
    #count number of segments for debugging
    all_segments = []
    if apic!=[]: 
        all_segments.extend(list(apic))
    for bas in basals:
        all_segments.extend(list(bas))

    if len(all_segments)!=len(reduced_seg_to_original_seg.keys()):
        print('There is no segment to segment copy, it means that some segments in the reduced model did not receive channels from the original cell')
        #print('there are {:} orphan segments out of {:}'.format(len(all_segments)-len(reduced_seg_to_original_seg.keys()),len(all_segments)))
        print('trying to compensate by copying channels from neighboring segments')
        handle_orphan_segments(original_seg_to_reduced_seg, all_segments, vals_per_mech_per_segment, mech_names_per_segment)
        #print('finish orphan segments fix')

def handle_orphan_segments(original_seg_to_reduced_seg, all_segments, vals_per_mech_per_segment, mech_names_per_segment):
    ''' This function handle reduced segments that did not had original segments mapped to them'''
    all_mapped_control_segments = original_seg_to_reduced_seg.values() # Get all reduced segments that have been mapped by a original model segment
    non_mapped_segments = set(all_segments) - set(all_mapped_control_segments)

    for reduced_seg in non_mapped_segments:
        seg_secs = list(reduced_seg.sec)
        #find valid parent
        parent_seg_index = seg_secs.index(reduced_seg)-1
        parent_seg = None
        while parent_seg_index>-1:
            if seg_secs[parent_seg_index] in all_mapped_control_segments:
                parent_seg = seg_secs[parent_seg_index]
                break
            else:
                parent_seg_index-=1
        
        #find valid child
        child_seg_index = seg_secs.index(reduced_seg)+1
        child_seg = None
        while child_seg_index<len(seg_secs):
            if seg_secs[child_seg_index] in all_mapped_control_segments:
                child_seg = seg_secs[child_seg_index]
                break
            else:
                child_seg_index+=1

        if not parent_seg and not child_seg:
            raise Exception("no child seg nor parent seg, with active channels, was found")
        
        if parent_seg and not child_seg:
            for mech in mech_names_per_segment[parent_seg]:
                reduced_seg.sec.insert(mech)
            for n in vals_per_mech_per_segment[parent_seg]:
                exec('reduced_seg.' + n +' = np.mean(vals_per_mech_per_segment[parent_seg][n])')
  
        if not parent_seg and child_seg:
            for mech in mech_names_per_segment[child_seg]:
                reduced_seg.sec.insert(mech)
            for n in vals_per_mech_per_segment[child_seg]:
                exec('reduced_seg.' + n +' = np.mean(vals_per_mech_per_segment[child_seg][n])')

        ## if both parent and child were found, we add to the segment all the mech in both
        ## this is just a decision
        
        if parent_seg and child_seg:
            for mech in set(mech_names_per_segment[child_seg]) & set(mech_names_per_segment[parent_seg]):
                reduced_seg.sec.insert(mech)
            
            for n in vals_per_mech_per_segment[child_seg]:
                if n in vals_per_mech_per_segment[parent_seg]:
                    exec('reduced_seg.' + n +' = (np.mean(vals_per_mech_per_segment[child_seg][n]) + ' + \
                            'np.mean(vals_per_mech_per_segment[parent_seg][n]))/2')
                else:
                    exec('reduced_seg.' + n +' = np.mean(vals_per_mech_per_segment[child_seg][n])')
            for n in vals_per_mech_per_segment[parent_seg]:
                if n in vals_per_mech_per_segment[child_seg]:
                    exec('reduced_seg.' + n +' = (np.mean(vals_per_mech_per_segment[child_seg][n]) + ' + \
                            'np.mean(vals_per_mech_per_segment[parent_seg][n]))/2')
                else:
                    exec('reduced_seg.' + n +' = np.mean(vals_per_mech_per_segment[parent_seg][n])')

    
    
    
        
def add_PP_properties_to_dict(PP,PP_params_dict):
    '''
    add the propeties of a point process to PP_params_dict. 
    The only propeties added to the dictionary are those worth comparing 
    '''  
    PP_params = []
    for param in dir(PP):
        if param[:2] =="__":
            continue
        if param in ["Section","allsec","baseattr","cas","g","get_loc","has_loc","hname",
                    'hocobjptr',"i","loc","next","ref","same","setpointer","state" , "get_segment"]:
            continue
        PP_params.append(param)
    PP_params_dict[type_of_point_process(PP)] = PP_params


def type_of_point_process(PP):
    s = PP.hname()
    ix = PP.hname().find("[")
    return s[:ix]

def subtree_reductor(original_cell, synapses_list, netcons_list, reduction_frequency, model_filename='model.hoc', total_segments_manual = -1,PP_params_dict = {}):
    '''
    Receives an instance of a cell with a loaded full morphology, a list of synapse objects, a list of NetCon objects (the i'th netcon in the list should correspond to the i'th synapse), 
    the filename (string) of the model template hoc file that the cell was instantiated from, the desired reduction frequency as a float, optional parameter for the approximate desired number of segments in the new model (if this parameter is empty, the number of segments will be such that there is a segment for every 0.1 lambda),
    and an optional param for the point process to be compared before deciding on whethet to merge a synapse or not
    and reduces the cell (using the given reduction_frequency). Creates a reduced instance using the model template in the file whose filename is given as a parameter, and merges synapses of the same type that get mapped to the same segment (same "reduced" synapse object for them all, but different NetCon objects).
    Returns the new reduced cell, a list of the new synapses, and the list of the inputted netcons which now have connections with the new synapses.
    note #0: a default template is available, one can use: model_filename=model.hoc
    note #1: The original cell instance, synapses and Netcons given as arguments are altered by the function and cannot be used outside of it in their original context.
    note #2: Synapses are determined to be of the same type and mergeable if their reverse potential, tau1 and tau2 values are identical.
    note #3: Merged synapses are assigned a single new synapse object that represents them all, but keep their original NetCon objects. Each such NetCon now connects the original synapse's NetStim with the reduced synapse.
    '''  
    h.init()
    start_time = time.time()
    global SectionsToDelete, section_per_subtree_index
    SectionsToDelete = []
    
    
    model_obj_name = model_filename.split(".")[0].split('/')[-1]
    if h.name_declared(model_obj_name)==0:
        print("loading template '" + model_obj_name +"'")
        if model_filename =='model.hoc':
            print("loading default reduced model")
            load_default_model()
        else:
            h.load_file(model_filename)
    else:
        print("WARNING The template '" +  model_obj_name + "' is already defined... not loading.")
    

    # finds soma properties
    soma = original_cell.soma[0] if original_cell.soma.hname()[-1]==']' else original_cell.soma
    soma_ref = h.SectionRef(sec=soma)
    
    SOMA_CM = soma.cm
    SOMA_RM = 1.0 / soma.g_pas
    SOMA_RA = soma.Ra
    SOMA_E_PAS= soma.e_pas
    soma_length = soma.L
    soma_diam = soma.diam
  
    has_apical = True
    if len(list(original_cell.apical)) == 0:
        has_apical = False

    AxonSection, AxonIsParent = find_and_diconnect_Axon(soma_ref) # This function find the axon, if the axon is the parent of the soma, AxonParent will be True


    roots_of_subtrees = [] # will hold the root sections of each of the soma's subtrees
    num_of_subtrees = [] # adding correctly the number of subtrees, excluding the axon
    for i in range(int(soma_ref.nchild())):
        num_of_subtrees.append(i)
        roots_of_subtrees.append(soma_ref.child[i])
    

    # assuming up to one apical tree
    ix_of_apical = None
    for i in num_of_subtrees:
        if 'apic' in roots_of_subtrees[i].hname():
            ix_of_apical = i
            
    if ix_of_apical:
        roots_of_subtrees = [roots_of_subtrees[ix_of_apical]]+roots_of_subtrees[:ix_of_apical]+roots_of_subtrees[ix_of_apical+1:]


    h("num_of_subtrees = -1")
    h.num_of_subtrees = len(num_of_subtrees)
    h("strdef string")



    mapping_sections_to_subtree_index = {}  # dict that maps section indexes to the subtree index they are in: keys are string tuples: ("apic"/"basal", orig_section_index) , values are ints: subtree_instance_index 
    for i,soma_child in enumerate(roots_of_subtrees):
        # inserts each section in this subtree into the above dict, which maps it to the subtree index
        if 'apic' in soma_child.hname():
            if i!=0:
                raise Exception('The apical is not the first child of the soma! a code refactoring is needed in order to accept it')
            mark_subtree_sections_with_subtree_index(soma_child, mapping_sections_to_subtree_index, "apic", i)
        elif 'dend' in soma_child.hname() or 'basal' in soma_child.hname():
            mark_subtree_sections_with_subtree_index(soma_child, mapping_sections_to_subtree_index, "basal", i)


    # mapping_sections_to_subtree_index = {}  # dict that maps section indexes to the subtree index they are in: keys are string tuples: ("apic"/"basal", orig_section_index) , values are ints: subtree_instance_index 
    # for i in num_of_subtrees:
    #     # inserts each section in this subtree into the above dict, which maps it to the subtree index
    #     if 'apic' in soma_ref.child[i].hname():
    #         if i!=0:
    #             raise Exception('The apical is not the first child of the soma! a code refactoring is needed in order to accept it')
    #         mark_subtree_sections_with_subtree_index(soma_ref.child[i], mapping_sections_to_subtree_index, "apic", i)
    #     elif 'dend' in soma_ref.child[i].hname() or 'basal' in soma_ref.child[i].hname():
    #         mark_subtree_sections_with_subtree_index(soma_ref.child[i], mapping_sections_to_subtree_index, "basal", i)

  
    # preparing for reduction
    
    # remove active conductances and get seg_to_mech dictionary
    segment_to_mech_vals = create_segments_to_mech_dic()
    
    # disconnects all the subtrees from the soma  
    for subtree_root in roots_of_subtrees:
        h.disconnect(sec = subtree_root)
    
    reduced_cables_dimensions = []  # the i'th element is a list that will hold the i'th new cable's length, diam, and space_const (in micron) 
    for i in num_of_subtrees:
        reduced_cables_dimensions.append([])
        
    biophysical_cable_properties = [] # the i'th element is a list that will hold the i'th subtree's trunk's cm, rm, ra, and e_pas  
    for i in num_of_subtrees:
        biophysical_cable_properties.append([])
        
    reduced_cables_electrotonic_length = [None] * len(num_of_subtrees)
    reduced_cables_electrotonic_length_as_complex_nums = [None] * len(num_of_subtrees)
  
  
    # reducing the subtrees
    set_up_for_reduction(reduction_frequency, has_apical)
    for i in num_of_subtrees:
        reduced_cables_dimensions[i] = reduce_subtree(roots_of_subtrees[i], i, reduced_cables_electrotonic_length, reduced_cables_electrotonic_length_as_complex_nums, biophysical_cable_properties[i])
    
    # calculating the new number of segments for each reduced dendritic cable
    new_cables_nsegs = [None] * len(num_of_subtrees)
    if total_segments_manual != -1:  # if a value for the total_num_of_segments was given as an argument to the program  (optional arg)
        calculate_nsegs_from_manual_arg(reduced_cables_electrotonic_length, new_cables_nsegs, total_segments_manual)
    else:
        calculate_nsegs_from_lambda(reduced_cables_dimensions, new_cables_nsegs) # the default nseg calculation - places a segment for every 0.1 lambda
    

    # creating reduced_cell
    h("objref reduced_cell")
    h("reduced_cell = new " + model_obj_name + "()")

    reduced_cell = h.reduced_cell
    create_sections_in_hoc("soma", 1, "reduced_cell")
    # create cell python template
    cell = Neuron(reduced_cell)
    cell.set_soma(original_cell.soma)   

    soma = original_cell.soma[0] if original_cell.soma.hname()[-1]==']' else original_cell.soma

    if has_apical: # creates reduced apical cable if apical subtree existed
        create_sections_in_hoc("apic", 1, "reduced_cell")
        apic = reduced_cell.apic[0]
        num_of_basal_subtrees = len(num_of_subtrees) - 1
    else:
        num_of_basal_subtrees = len(num_of_subtrees)

    # creates reduced basal cables 
    create_sections_in_hoc("dend", num_of_basal_subtrees, "reduced_cell")
    basals = []
    for i in range(num_of_basal_subtrees):
        basals.append(reduced_cell.dend[i])
    
    append_to_section_lists("soma[0]", "somatic", "reduced_cell")
    
    if has_apical:  
        apic.L = reduced_cables_dimensions[0][0]  # the apical cable is index 0 in reduced_cables_dimensions and in the below biophysical_cable_properties
        apic.diam = reduced_cables_dimensions[0][1]
        apic.nseg = new_cables_nsegs[0]
        append_to_section_lists("apic[0]", "apical", "reduced_cell")
        connect_to(apic, soma, 1, 0)
        apic.insert('pas')
        apic.cm = biophysical_cable_properties[0][0]
        apic.g_pas= 1.0 / biophysical_cable_properties[0][1] 
        apic.Ra = biophysical_cable_properties[0][2]
        apic.e_pas = biophysical_cable_properties[0][3]
 
    for i in range(num_of_basal_subtrees):
        if has_apical:
            index_in_reduced_cables_dimensions = i + 1
        else:  
            index_in_reduced_cables_dimensions = i
      
        basals[i].L = reduced_cables_dimensions[index_in_reduced_cables_dimensions][0]
        basals[i].diam = reduced_cables_dimensions[index_in_reduced_cables_dimensions][1]
        basals[i].nseg = new_cables_nsegs[index_in_reduced_cables_dimensions]
        curr_dend_as_str = "dend["+ str(i) + "]"
        append_to_section_lists(curr_dend_as_str, "basal", "reduced_cell")
        connect_to(basals[i], soma, 0, 0)
        basals[i].insert('pas')
        basals[i].cm = biophysical_cable_properties[index_in_reduced_cables_dimensions][0]
        basals[i].g_pas= 1.0 / biophysical_cable_properties[index_in_reduced_cables_dimensions][1] 
        basals[i].Ra = biophysical_cable_properties[index_in_reduced_cables_dimensions][2]
        basals[i].e_pas = biophysical_cable_properties[index_in_reduced_cables_dimensions][3]



        
    
   # dividing the original synapses into baskets, so that all synapses that are on the same subtree will be together in the same basket
    baskets = []  # a list of baskets of synapses, each basket in the list will hold the synapses of the subtree of the corresponding basket index
    somatic_synapses = [] # each somatic or axonal synapse will be stored here
    soma_synapses_syn_to_netcon = {}
    location_in_basket_to_orig_index_synapse_map = {}  # dict that maps synapses in baskets to their orig synapse index: keys are tuples - (basket_index, index_in_basket), values are the original synapse index
    for i in num_of_subtrees:
        baskets.append([])
    for syn_index, synapse in enumerate(synapses_list):
        curr_subtree_index, curr_section_num, x = find_synapse_loc(synapse, mapping_sections_to_subtree_index)
        if curr_subtree_index == SOMA_LABEL or curr_subtree_index == 'axon':  # for a somatic synapse
            somatic_synapses.append(synapse)
            soma_synapses_syn_to_netcon[synapse] = netcons_list[syn_index]
        else:
            baskets[curr_subtree_index].append(synapse)
            location_in_basket_to_orig_index_synapse_map[(curr_subtree_index,len(baskets[curr_subtree_index]) -1)] = syn_index
    
    # mapping (non-somatic) synapses to their new location on the reduced model (the new location is the exact location of the middle of the segment they were mapped to, in order to enable merging)
    new_merged_synapse_locs = {}  # holds merged loc of non-somatic synapses and their "type" signature - (reverse potential, tau1, tau2) [as a tuple holding two tuples]
    new_synapses_list = []
    new_synapse_counter = 0
  
    subtree_ind_to_q = {}
    for i in num_of_subtrees:
        imp_obj, subtree_input_impedance = measure_input_impedance_of_subtree(roots_of_subtrees[i])
        #calculating the q of this subtree:
        rm = 1.0/roots_of_subtrees[i].g_pas
        rc = rm * (float(roots_of_subtrees[i].cm) / 1000000)
        angular_freq = 2 * math.pi * reduction_frequency 
        q_imaginary = angular_freq * rc
        q_subtree = complex(1, q_imaginary)   # q=1+iwRC
        q_subtree = cmath.sqrt(q_subtree)
        subtree_ind_to_q[i] = q_subtree
        # iterates over the synapses in the curr basket
        for i_in_basket, synapse in enumerate(baskets[i]):
            curr_subtree_index, curr_section_num, x = find_synapse_loc(synapse, mapping_sections_to_subtree_index)
            if has_apical and curr_subtree_index == 0: # if synapse is on the apical subtree
                on_basal_subtree = False
            else:
                on_basal_subtree = True
        
            # "reduces" the synapse - finds this synapse's new "merged" location on its corresponding reduced cable
            x = reduce_synapse(original_cell, curr_section_num, x, on_basal_subtree, new_cables_nsegs[curr_subtree_index], imp_obj, 
                               subtree_input_impedance, reduced_cables_electrotonic_length[curr_subtree_index],
                               reduced_cables_electrotonic_length_as_complex_nums[curr_subtree_index], q_subtree)
            
            reduced_cable_index = curr_subtree_index
            # find the section of the synapse
            if on_basal_subtree:
                if has_apical:
                    section_for_synapse = basals[reduced_cable_index-1]
                else:
                    section_for_synapse = basals[reduced_cable_index]
            else:
                section_for_synapse = apic
            # find the segment of the synapse
            seg_pointer = section_for_synapse(x)
            PPs_on_seg = seg_pointer.point_processes()
            # go over all point processes in this seg,ent and see whther one of them has the same proporties of this synapse
            # If there's such a synapse link the original NetCon with this point processes
            # If not, move the synapse to this segment.
            for PP in PPs_on_seg:
                if type_of_point_process(PP) not in PP_params_dict:
                    add_PP_properties_to_dict(PP,PP_params_dict)
                    
                try:
                    for param in PP_params_dict[type_of_point_process(PP)]:
                        if param not in ['rng']: # See this ticket :https://github.com/neuronsimulator/nrn/issues/136
                            if str(type(getattr(PP,param)))!="<type 'hoc.HocObject'>": #ignore hoc objects
                                if getattr(PP,param)!=getattr(synapse,param):
                                    break
                    else: # same synapse propeties: merge them.
                        netcons_list[location_in_basket_to_orig_index_synapse_map[(curr_subtree_index, i_in_basket)]].setpost(PP)
                        break

                except:
                    continue
            else: #If for finish the loop, it means that it is the first appearance of this synapse
                synapse.loc(x, sec = section_for_synapse)
                new_synapses_list.append(synapse)
                new_synapse_counter = new_synapse_counter + 1
    
    # merging somatic and axonal synapses
    synapses_per_seg = {}
    for synapse in somatic_synapses:
        seg_pointer = synapse.get_segment()
        x = seg_pointer.x
        section_for_synapse = seg_pointer.sec
        if seg_pointer not in synapses_per_seg: 
            synapses_per_seg[seg_pointer] = []
        PPs_on_seg = synapses_per_seg[seg_pointer]
        for PP in PPs_on_seg:
            if type_of_point_process(PP) not in PP_params_dict:
                add_PP_properties_to_dict(PP,PP_params_dict)
                
            try:
                for param in PP_params_dict[type_of_point_process(PP)]:
                    if param not in ['rng']: # See this ticket :https://github.com/neuronsimulator/nrn/issues/136
                        if str(type(getattr(PP,param)))!="<type 'hoc.HocObject'>": #ignore hoc objects
                            if getattr(PP,param)!=getattr(synapse,param):
                                break
                else: # same synapse propeties: merge them.
                    soma_synapses_syn_to_netcon[synapse].setpost(PP)
                    break

            except:
                continue
        else: #If for finish the loop, it means that it is the first appearance of this synapse
            synapse.loc(x, sec = section_for_synapse)
            new_synapses_list.append(synapse)
            new_synapse_counter = new_synapse_counter + 1
            synapses_per_seg[seg_pointer].append(synapse)

    # copy active mechanisms
    if has_apical==False: apic=[]
    copy_dendritic_mech(original_cell, roots_of_subtrees, mapping_sections_to_subtree_index, new_cables_nsegs, reduced_cables_electrotonic_length,  
                      reduced_cables_electrotonic_length_as_complex_nums, has_apical, apic, basals, segment_to_mech_vals, subtree_ind_to_q, mapping_type ='impedance')
    
    # Connect axon back to the soma
    if len(AxonSection)>0:
        if AxonIsParent ==True:
            soma.connect(AxonSection[0])
        else:
            AxonSection[0].connect(soma)
    
    ## Now we delete the original model, ToDelete is a global.
    for i in range(len(SectionsToDelete)):
        SectionsToDelete[i].push()
        h.delete_section()
        h.pop_section()

    #set axon of cell to the original axon
    cell.set_axon(AxonSection)
    
    #set dend of cell to the simplified dend
    cell.set_dend(reduced_cell.dend)
    
    #set apic of cell to the simplified apic
    if apic!=[]:
        cell.set_apic(reduced_cell.apic)
    
    # delete reduced soma which was created 
    reduced_cell.soma[0].push()
    h.delete_section()
    
    # reset globals
    SectionsToDelete = [] # A list of all section that will be deleted in the end, it is easier with a global.
    section_per_subtree_index = {}
    
    secs_elapsed = time.time() - start_time
    print("reduction time in seconds: ", secs_elapsed)
    return cell, new_synapses_list, netcons_list


class Neuron:
    'Python neuron class for hoc models'
    def __init__(self, model):
        self.hoc_model = model
    
    def set_soma(self, soma):
        self.soma = soma
        
    def set_dend(self, dend):
        self.dend = dend
    
    def set_apic(self, apic):
        self.apic = apic
    
    def set_axon(self, axon):
        self.axon = axon



def load_default_model():
    h('''begintemplate model

public init, biophys, geom_nseg, delete_axon, finish_creating_model_after_loading_morphology

public soma, dend, apic, axon  // sections
public all, somatic, apical, axonal, basal // section lists
objref all, somatic, apical, axonal, basal, this

proc init() {
    all = new SectionList()
    somatic = new SectionList()
    basal = new SectionList()
    apical = new SectionList()
    axonal = new SectionList()
    
    forall delete_section()
    StepDist = 60 // human cells have no spines in their first 60 um
                                // from soma - see Benavides-Piccione 2013
    F_Spines = 1.9       //As calculated - see detailes in Eyal 2015
    
    CM =0.45	// uF/cm2
    RM = 38907
    RA = 203
    E_PAS =  -86
    
}

create soma[1], dend[1], apic[1], axon[1]   

//external lambda_f
proc geom_nseg() {
    soma distance() 
    
    forsec all {
        RA_calc = RA
        RM_calc = RM*F_Spines
        if (distance(1)>StepDist){
            RA_calc = RA
            RM_calc = RM*F_Spines
        }
        d = diam
        lambda = sqrt(RM_calc/RA_calc*d/10000/4)*10000
        nseg = int(L/lambda*10/2)*2+1		
    }
}


proc biophys() {
    forsec all {
        insert pas
        cm =CM
        g_pas=1/RM
        Ra = RA
        e_pas = E_PAS
    }
  
    soma distance()
  
    forsec basal {	
        if (distance(0.5)>StepDist) {
            L = L*F_Spines^(2/3)
            diam = diam*(F_Spines^(1/3))
        }	
    }
    forsec apical {
        if (distance(0.5)>StepDist) {
            L = L*F_Spines^(2/3)
            diam = diam*(F_Spines^(1/3))
        }	
    }
}


proc delete_axon(){
    forsec axonal{delete_section()}
}


proc complete_full_model_creation() {
    geom_nseg()      		             // calculates num of segments
    delete_axon()		                     // deletes the axon
    biophys()			             // increases cell dimensions to account for spines
}

endtemplate model''')
