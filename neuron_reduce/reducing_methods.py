# This file contains the reduction algorithm itself
# added the method by Guy to find L and X
import cmath
import math
from fractions import Fraction
from neuron import h
import os

h('''obfunc lowest_impedance_recursive() { local lowest_impedance, lowest_phase, i   localobj curr_subtree_root, sref1, lowest_imp_vec, lowest_child_subtree_impedance, imp_obj
    curr_subtree_root = $o1  // in the first call to the function, this is a root section of a dendritic trunk
    imp_obj = $o2
    curr_subtree_root.sec {
        lowest_impedance = imp_obj.transfer(1) // farthest tip of the the curr root section
        lowest_phase = imp_obj.transfer_phase(1)                                        
    }
  
    if (curr_subtree_root.nchild != 0) { // if the curr section has child sections
        for i=0, curr_subtree_root.nchild-1 curr_subtree_root.child[i] {  // for each child of the root, finds the lowest impedance within the subtree whose root is the curr child (in relation to the proximal tip in the curr root child)
            curr_subtree_root.child[i] sref1 = new SectionRef()
            lowest_child_subtree_impedance = lowest_impedance_recursive(sref1, imp_obj) // recursively returns the lowest transfer impedance and transfer phase within the curr subtree as a vector
            if (lowest_child_subtree_impedance.x[0] < lowest_impedance) {                             
                lowest_impedance = lowest_child_subtree_impedance.x[0]
                lowest_phase = lowest_child_subtree_impedance.x[1]
            }
        }
    }
    lowest_imp_vec = new Vector(2)
    lowest_imp_vec.x[0] = lowest_impedance
    lowest_imp_vec.x[1] = lowest_phase
    return lowest_imp_vec
}''')

RC = None
RM = None
RA = None
q = None

has_apical = None
frequency = None


def set_up_for_reduction(freq, has_apic_subtree):
    ''' sets global parameters to be used during the entire reduction process '''
    global frequency, has_apical
    frequency = freq  # in Hz
    has_apical = has_apic_subtree


def set_subtree_biophysical_properties(subtree_root_ref, biophysical_properties_list):
    ''' sets the biophysical cable properties (Rm, Ra, Rc) and q for the subtree to be reduced according to the properties of the root section of the subtree '''
    subtree_root_section = subtree_root_ref.sec
    
    global RC, q, RM, RA
    RM = 1.0 / subtree_root_section.g_pas  # in ohm * cm^2
    RA = subtree_root_section.Ra           # in ohm * cm
    RC = RM * (float(subtree_root_section.cm) / 1000000)   # in secs, with conversion of the capacitance from uF/cm2 to F/cm2  
    
    # defining q=sqrt(1+iwRC))
    angular_freq = 2 * math.pi * frequency   # = w                  
    q_imaginary = angular_freq * RC
    q = complex(1, q_imaginary)   # q=1+iwRC
    q = cmath.sqrt(q)		# q = sqrt(1+iwRC)
    
    # stores this subtree's biophysical properties
    biophysical_properties_list.append(subtree_root_section.cm)
    biophysical_properties_list.append(RM)
    biophysical_properties_list.append(RA)
    biophysical_properties_list.append(subtree_root_section.e_pas)


def find_lowest_subtree_impedance(subtree_root_ref, imp_obj): 
    ''' 
    finds the segment in the subtree with the lowest transfer impedance in relation to the proximal-to-soma end of the given subtree root section, using a recursive hoc function, 
    and returns the lowest impedance in Ohms
    '''
    lowest_impedance = h.lowest_impedance_recursive(subtree_root_ref, imp_obj) # returns [lowest subtree transfer impedance in Mohms, transfer phase]
    curr_lowest_subtree_imp = cmath.rect(lowest_impedance.x[0] * 1000000, lowest_impedance.x[1])  # impedance saved as a complex number after converting Mohms to ohms
    return curr_lowest_subtree_imp




def compute_zl_polar(Z0,L,q):
    '''
    given Z0 , L and q computes the polar represntation of ZL accoridng to equation 2.9 in Gals thesis
    '''
    ZL = Z0*1.0/cmath.cosh(q*L)
    ZL = cmath.polar(ZL)
    return ZL

def find_best_real_L(Z0,ZL_goal,q,max_L=10.0):
    '''
    This function finds the best real L s.t. the modulus part of the impedance of ZL in eq 2.9 will be correct
    Since the modulus is a decreasing function of L, it is easy to find it using binary search. 
    '''
    min_L = 0.0
    current_L = (min_L+max_L)/2.0
    ZL_goal_A = cmath.polar(ZL_goal)[0]
    Z_current_L_A = compute_zl_polar(Z0,current_L,q)[0]
    depth = 0 
    while abs(ZL_goal_A-Z_current_L_A)>0.001 and depth<50: # Z are in Ohms , normal values are >10^6 
        if ZL_goal_A>Z_current_L_A:
            current_L,max_L = (min_L+current_L)/2.0,current_L
        if ZL_goal_A<Z_current_L_A:
            current_L,min_L = (max_L+current_L)/2.0,current_L
        Z_current_L_A = compute_zl_polar(Z0,current_L,q)[0]
        depth += 1
    if depth == 50:
        print("The difference between L and the goal L is larger than 0.001")
    return current_L


def compute_zx_polar(Z0,L,q,x):
    '''
    given Z0 , L , q and x computes the polar represntation of Zx accoridng to equation 2.8 in Gals thesis
    '''
    try:
        ZX = Z0*cmath.cosh(q*(L-x))/cmath.cosh(q*L)
        ZX = cmath.polar(ZX)
    except:
        pdb.set_trace()
    return ZX


def find_best_real_X(Z0,ZX_goal,q,L):
    '''
    This function finds the best location of a suncapse (X) s.t. the modulus part of the impedance of ZX
    in eq 2.8 will be correct.
    Since the modulus is a decreasing function of L, it is easy to find it using binary search. 
    '''


    min_X = 0.0
    max_x = L
    current_x = (min_X+max_x)/2.0
    ZX_goal = cmath.polar(ZX_goal)[0]
    Z_current_X_A = compute_zx_polar(Z0,L,q,current_x)[0]
    depth = 0 
    while abs(ZX_goal-Z_current_X_A)>0.001 and depth<50:

        if ZX_goal>Z_current_X_A:
            current_x,max_x = (min_X+current_x)/2.0,current_x
        if ZX_goal<Z_current_X_A:
            current_x,min_X = (max_x+current_x)/2.0,current_x
        Z_current_X_A = compute_zx_polar(Z0,L,q,current_x)[0]
        depth += 1
    if depth == 50:
        print("The difference between X and the goal X is larger than 0.001")
    return current_x

def find_subtree_new_electrotonic_length(root_input_impedance, lowest_subtree_impedance):
    '''
    finds the subtree's reduced cable's electrotonic length according to the given complex impedance values, based on the following equation:
    lowest_subtree_impedance = subtree_root_input_impedance/cosh(q*L)
    '''
    
    # this equation could be solved analytically using L = 1/q * arcosh(subtree_root_input_impedance/lowest_subtree_impedance),
    # But since L in this equation is complex number and we chose to focus on finding the correct attenuation
    # we decided to search the L that will result with correct attenuation from the tip of the dendrite to the soma.
    # We chose to use only real L (without a complex part)

    L = find_best_real_L(root_input_impedance,lowest_subtree_impedance,q)
    return L, complex(L,0)


def find_subtree_new_diam_in_cm(Root_input_impedance, electrotonic_length_as_complex):
    '''
    finds the subtree's new cable's diameter (in cm) according to the given complex input impedance at the segment in the original subtree that is closest to the soma (the tip), and the given cable electrotonic length,
    with the following equation: d (in cm) = (2/PI * (sqrt(RM*RA)/(q*subtree_root_input_impedance)) * (coth(q * NewCableElectrotonicLength)) )^(2/3) , derived from Rall's cable theory for dendrites (Gal Eliraz)
    '''  
    
    diam_in_cm = (2.0/math.pi * (math.sqrt(RM*RA)/(q*Root_input_impedance)) * (1/cmath.tanh(q * electrotonic_length_as_complex)) ) ** (2.0/3)   # coth = 1/tanh
    
    '''
    # for debugging inaccuracies:
    if diam_in_cm.imag != 0:
        if abs(diam_in_cm.imag) > 0.03:
        print "PROBLEM - DIAM HAS SUBSTANTIAL IMAGINARY PART"
        print "\n"
    '''
    
    new_subtree_dend_diam_in_cm = cmath.polar(diam_in_cm)[0]  # the radius of the complex number received from the equation
    return new_subtree_dend_diam_in_cm


def find_space_const_in_cm(currNewCableDiameterInCm):
    ''' returns space constant (lambda) in cm, according to: space_const = sqrt(rm/(ri+r0)) '''
    diameter = currNewCableDiameterInCm              # receives diam in cm
    rm = float(RM) / (math.pi * diameter)            # rm = Rm/(PI * diam), diam is in cm and Rm is in ohm * cm^2
    ri = float(4 * RA) / (math.pi * (diameter**2))   # ri = 4*Ra/ (PI * diam^2), diam is in cm and Ra is in ohm * cm
    space_const = math.sqrt(rm/ri)                   # r0 is negligible
    return space_const   
    

def reduce_subtree(subtree_root, subtree_index, new_cables_electrotonic_lengths, new_cables_electrotonic_length_as_complex_nums, biophysical_properties_list):        
    '''
    Reduces the subtree with the given subtree_index from the original_cell into one single section (cable).
    The reduction is done by finding the length and diameter of the cable (a single solution) that preserves the subtree's input impedance at the somatic end, and the transfer impedance in the subtree from the distal end to the proximal somatic end (between the new cable's two tips).
      - the new electrotonic length of the cable is stored as a real number and as a complex number in the given lists new_cables_electrotonic_lengths and new_cables_electrotonic_length_as_complex_nums, according to the cable's subtree index.
    '''

    subtree_root_ref = h.SectionRef(sec = subtree_root)
    set_subtree_biophysical_properties(subtree_root_ref, biophysical_properties_list)
    
    # finds the subtree's input impedance (at the somatic-proximal end of the subtree root section) and the lowest transfer impedance in the subtree in relation to the somatic-proximal end (see more in Readme on NeuroReduce)
    imp_obj, root_input_impedance = measure_input_impedance_of_subtree(subtree_root)
    curr_lowest_subtree_imp = find_lowest_subtree_impedance(subtree_root_ref, imp_obj)   # in Ohms (a complex number)    
    ##print "input impedance and lowest transfer impedance: ", subtree_index, repr(root_input_impedance), repr(curr_lowest_subtree_imp)  
    
    # reducing the whole subtree into one section:   L = 1/q * arcosh(ZtreeIn(f)/min(ZtreeX,0(f)), d = ( (2/pi * (sqrt(Rm*Ra)/q*ZtreeIn(f)) * coth(qL) )^(2/3) - from Gal Eliraz's thesis 1999
    new_cable_electrotonic_length, cable_electrotonic_length_as_complex = find_subtree_new_electrotonic_length(root_input_impedance, curr_lowest_subtree_imp)
    new_cable_diameter_in_cm = find_subtree_new_diam_in_cm(root_input_impedance, cable_electrotonic_length_as_complex)
    new_cable_diameter = new_cable_diameter_in_cm * 10000   #in microns
     
    #calculating the space constant, in order to find the cylinder's length:  space_const = sqrt(rm/(ri+r0)) 
    curr_space_const_in_cm = find_space_const_in_cm(new_cable_diameter_in_cm) 
    curr_space_const_in_micron = 10000 * curr_space_const_in_cm
    new_cable_length = curr_space_const_in_micron * new_cable_electrotonic_length  # in microns 	

    # storing the cable's electrotonic length
    new_cables_electrotonic_lengths[subtree_index] = new_cable_electrotonic_length
    new_cables_electrotonic_length_as_complex_nums[subtree_index] = cable_electrotonic_length_as_complex    
    
    newCableDimensions = [new_cable_length, new_cable_diameter, curr_space_const_in_micron] # all in microns
    return newCableDimensions
    
    
def find_merged_loc(cable_nseg, relative_loc):
    '''
    Returns a synapse's merged relative location (x) on the cable, according to its given relative location on the cable and the given number of segments in the cable.
    The merged location is the relative location of the middle of the segment the synapse is in (or 0 or 1 if it is at one of the tips of the cable).
    '''
    
    if relative_loc == 0 or relative_loc == 1:
        return relative_loc
      
    mapped_segment_for_curr_syn = int(relative_loc*cable_nseg)+1  # finds the segment that the synapse is in, according to its relative location and the num of segments in the cable (1 through nseg)
    return ( (float(mapped_segment_for_curr_syn)/cable_nseg) + (float(mapped_segment_for_curr_syn-1)/cable_nseg) ) / 2   # location of middle of segment = the average between relative location of end of segment and relative location of beginning of segment
        

def measure_input_impedance_of_subtree(subtree_root_section):
    ''' measures the input impedance of the subtree with the given root section (at the "0" tip, the soma-proximal end), returns the Impedance hoc object and the input impedance as a complex value'''
    
    imp_obj = h.Impedance()	
    CLOSE_TO_SOMA_EDGE = 0
    imp_obj.loc(CLOSE_TO_SOMA_EDGE, sec = subtree_root_section)  # sets origin for impedance calculations (soma-proximal end of root section)
    imp_obj.compute(frequency+1/9e9, 0)                                # computes transfer impedance from every segment in the model in relation to the origin location above
    
    root_input_impedance = imp_obj.input(CLOSE_TO_SOMA_EDGE, sec = subtree_root_section) * 1000000  # in Ohms (impedance measured at soma-proximal end of root section)
    root_input_phase = imp_obj.input_phase(CLOSE_TO_SOMA_EDGE, sec = subtree_root_section)                              
    root_input_impedance = cmath.rect(root_input_impedance, root_input_phase)                       # creates a complex impedance value out of the given polar coordinates
    return imp_obj, root_input_impedance
    
def reduce_synapse(cell_instance, section_num, x, on_basal, nseg_of_this_cable, imp_obj, root_input_impedance, 
                   new_cable_electrotonic_length, new_cable_electrotonic_length_as_complex, q_subtree):
    '''
    Receives an instance of a cell, the location (section + relative location(x)) of a synapse to be reduced, a boolean on_basal that is True if the synapse is on a basal subtree, 
    the number of segments in the reduced cable that this synapse is in, an Impedance calculating Hoc object, the input impedance at the root of this subtree, and the electrotonic length of the reduced cable that represents the current subtree (as a real and as a complex number) -
    and maps the given synapse to its new location on the reduced cable according to the NeuroReduce algorithm.
    Returns the new "post-merging" relative location of the synapse on the reduced cable (x, 0<=x<=1), that represents the middle of the segment that this synapse is located at in the new reduced cable.
    '''
    
    synapse_section_num = section_num
    synapse_relative_loc_in_section = x
  
    # measures the original transfer impedance from the synapse to the somatic-proximal end in the subtree root section
    if not on_basal:  # apical subtree
        cell_instance.apic[synapse_section_num].push()
    else:             # basal subtree
        cell_instance.dend[synapse_section_num].push()	  	
    synapse_orig_transfer_imp = imp_obj.transfer(synapse_relative_loc_in_section) * 1000000  # in ohms
    synapse_orig_transfer_phase = imp_obj.transfer_phase(synapse_relative_loc_in_section)
    orig_synapse_transfer_impedance = cmath.rect(synapse_orig_transfer_imp, synapse_orig_transfer_phase) # creates a complex Impedance value with the given polar coordinates
    h.pop_section()
    
    # synapse location could be calculated using: X = L - (1/q) * arcosh( (Zx,0(f) / ZtreeIn(f)) * cosh(q*L) ), derived from Rall's cable theory for dendrites (Gal Eliraz)
    # but we chose to find the X that will give the correct modulus. See comment about L values
    
    synapse_new_electrotonic_location = find_best_real_X(root_input_impedance,orig_synapse_transfer_impedance,q_subtree,new_cable_electrotonic_length)
    new_relative_loc_in_section = float(synapse_new_electrotonic_location) / new_cable_electrotonic_length
  
    if new_relative_loc_in_section>1:  # PATCH
        new_relative_loc_in_section = 0.999999
  
    return new_relative_loc_in_section
    

   
