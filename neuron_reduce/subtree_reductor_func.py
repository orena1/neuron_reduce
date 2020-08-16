'''
function subtree_reductor():
 which reduces a morphologically detailed cell instance into a morphologically
 simplified cell instance, according to NeuroReduce and merges synapses of the
 same type (same reverse potential, tau1, and tau2) that are mapped to the same
 segment. (see more in Readme on tool and usage)

 usage: For details, see comments in function

 outputs: reduced cell instance, a new synapses_list, and the netcons_list,
          which now corresponds to the new synapses.

- The model template file must have an init() function (see example in the
  attached model.hoc file) and the following public definitions specifying
  sections and section lists accordingly:
   public soma, dend, apic ; public all, somatic, apical, basal

- Supports numerous types of synapses (two synapses are considered to be of
  different types if they are different from each other in at least one of the
  following values: reverse potential, tau1, tau2)
'''
import collections
import itertools as it
import logging
import math
import re
import cmath

import numpy as np
import neuron
from neuron import h
h.load_file("stdrun.hoc")

from .reducing_methods import (reduce_subtree,
                               reduce_synapse,
                               measure_input_impedance_of_subtree,
                               CableParams,
                               SynapseLocation,
                               push_section,
                               )

logger = logging.getLogger(__name__)
SOMA_LABEL = "soma"
EXCLUDE_MECHANISMS = ('pas', 'na_ion', 'k_ion', 'ca_ion', 'h_ion', 'ttx_ion', )


def create_sections_in_hoc(type_of_section, num, instance_as_str):
    '''creates sections in the hoc world according to the given section type and number of sections

    in the instance whose name is given as a string
    '''
    h("strdef string")
    h.string = type_of_section
    h('{sprint(string, "create %s[%d]", string, ' + str(num) + ') }')
    h("{execute(string, " + instance_as_str + ")}")


def append_to_section_lists(section, type_of_sectionlist, instance_as_str):
    ''' appends given section to the sectionlist of the given type and to the "all" sectionlist

    in the hoc world in the instance whose name is given as a string
    '''
    h("strdef string")
    h.string = section + " " + type_of_sectionlist + ".append()"
    h("{execute(string, " + instance_as_str + ")}")
    h.string = section + " all.append()"
    h("{execute(string, " + instance_as_str + ")}")


def find_section_number(section):
    ''' extracts and returns the section number from the given section object '''
    sec_name = h.secname(sec=section)
    ints_in_name = re.findall(r'\d+', sec_name)
    sec_num = ints_in_name[len(ints_in_name) - 1]  # extracts section number
    return sec_num


def calculate_nsegs_from_manual_arg(new_cable_properties, total_segments_wanted):
    '''Calculates the number of segments for each section in the reduced model

    according to the given total_segments_wanted and the given
    new_dends_electrotonic_length (the electrotonic lengths of all the new
    sections).  Called when the user chooses to give to the program the
    approximate total number of segments that the reduced model should have
    (non-default calculation).
    '''
    # minus one for the one segment of the soma:
    total_segments_in_dendrites = total_segments_wanted - 1

    # total electrotonic length of reduced dendritic cables
    sum_of_lengths = sum(prop.electrotonic_length
                         for prop in new_cable_properties)

    # the num of segments assigned to each section is in proportion to the
    # section's relative contribution to the total electrotonic length in the
    # model
    dends_nsegs = []
    for prop in new_cable_properties:
        new_length = prop.electrotonic_length
        new_nseg_to_put = int(round((float(new_length) / sum_of_lengths) *
                              total_segments_in_dendrites))
        if new_nseg_to_put < 1:
            new_nseg_to_put = 1
        dends_nsegs.append(new_nseg_to_put)
    return dends_nsegs


def calculate_nsegs_from_lambda(new_cable_properties):
    '''calculate the number of segments for each section in the reduced model

    according to the length (in microns) and space constant (= lambda - in
    microns) that were previously calculated for each section and are given in
    subtree_dimensions.  According to this calculation, a segment is formed for
    every 0.1 * lambda in a section. (lambda = space constant = electrotonic length unit).
    '''
    dends_nsegs = []
    for cable in new_cable_properties:
        # for every unit of electronic length (length/space_constant such units)
        # ~10 segments are formed
        dends_nsegs.append(int((float(cable.length) / cable.space_const) * 10 / 2) * 2 + 1)
    return dends_nsegs


def mark_subtree_sections_with_subtree_index(sections_to_delete,
                                             section_per_subtree_index,
                                             root_sec_of_subtree,
                                             mapping_sections_to_subtree_index,
                                             section_type,
                                             subtree_index):
    '''Recursively marks all sections in the subtree as belonging to the given subtree_index

    using the given dict mapping_sections_to_subtree_index, as follows:
    mapping_sections_to_subtree_index[(<section_type>, <section_number>)] = given subtree_index
    '''
    sections_to_delete.append(root_sec_of_subtree)
    section_per_subtree_index.setdefault(subtree_index, [])
    section_per_subtree_index[subtree_index].append(root_sec_of_subtree)

    section_num = find_section_number(root_sec_of_subtree)

    for child in root_sec_of_subtree.children():
        mark_subtree_sections_with_subtree_index(sections_to_delete,
                                                 section_per_subtree_index,
                                                 child,
                                                 mapping_sections_to_subtree_index,
                                                 section_type,
                                                 subtree_index)
    mapping_sections_to_subtree_index[(section_type, section_num)] = subtree_index


def find_synapse_loc(synapse_or_segment, mapping_sections_to_subtree_index):
    ''' Returns the location  of the given synapse object'''

    if not isinstance(synapse_or_segment, neuron.nrn.Segment):
        synapse_or_segment = synapse_or_segment.get_segment()

    x = synapse_or_segment.x

    with push_section(synapse_or_segment.sec):
        # extracts the section type ("soma", "apic", "dend") and the section number
        # out of the section name
        full_sec_name = h.secname()
        sec_name_as_list = full_sec_name.split(".")
        short_sec_name = sec_name_as_list[len(sec_name_as_list) - 1]
        section_type = short_sec_name.split("[")[0]
        section_num = re.findall(r'\d+', short_sec_name)[0]

    # finds the index of the subtree that this synapse belongs to using the
    # given mapping_sections_to_subtree_index which maps sections to the
    # subtree indexes that they belong to
    if section_type == "apic":
        subtree_index = mapping_sections_to_subtree_index[("apic", section_num)]
    elif section_type == "dend":
        subtree_index = mapping_sections_to_subtree_index[("basal", section_num)]
    else:  # somatic synapse
        subtree_index, section_num, x = SOMA_LABEL, 0, 0

    return SynapseLocation(subtree_index, int(section_num), x)


def find_and_disconnect_axon(soma_ref):
    '''Searching for an axon, it can be a child of the soma or a parent of the soma.'''
    axon_section, axon_parent, soma_axon_x  = [], False, None

    for sec in soma_ref.child:
        name = sec.hname().lower()
        if 'axon' in name or 'hill' in name:
            axon_section.append(sec)
            # disconnect axon
            soma_axon_x = sec.parentseg().x
            sec.push()
            h.disconnect()
            h.define_shape()

    if soma_ref.has_parent():
        name = soma_ref.parent().sec.hname().lower()
        if 'axon' in name or 'hill' in name:
            axon_section.append(soma_ref.parent())
            axon_parent = True
            soma_axon_x = None
            soma_ref.push()
            h.disconnect()
        else:
            raise Exception('Soma has a parent which is not an axon')

    if len(axon_section) > 1:
        raise Exception('Soma has a two axons')

    return axon_section, axon_parent, soma_axon_x


def create_segments_to_mech_vals(sections_to_delete,
                                 remove_mechs=True,
                                 exclude=EXCLUDE_MECHANISMS):
    '''This function copy the create a mapping between a dictionary and the mechanisms that it have
       plus the values of those mechanisms. It also remove the mechanisms from the model in order to
       create a passive model

       Arguments:
           remove_mechs - False|True
               if True remove the mechs after creating the mapping, False - keep the mechs
           exclude - List of all the mechs name that should not be removed
       '''
    exclude = set(exclude)
    segment_to_mech_vals, mech_names = {}, set()

    for seg in it.chain.from_iterable(sections_to_delete):
        segment_to_mech_vals[seg] = {}
        for mech in seg:
            mech_name = mech.name()
            segment_to_mech_vals[seg][mech_name] = {}
            for n in dir(mech):
                if n.startswith('__') or n in ('next', 'name', 'is_ion', 'segment', ):
                    continue

                if not n.endswith('_' + mech_name) and not mech_name.endswith('_ion'):
                    n += '_' + mech_name

                segment_to_mech_vals[seg][mech_name][n] = getattr(seg, n)
                mech_names.add(mech_name)

    mech_names -= exclude

    if remove_mechs:  # Remove all the mechs from the sections
        for sec in sections_to_delete:
            with push_section(sec):
                for mech in mech_names:
                    h("uninsert " + mech)

    return segment_to_mech_vals


def create_seg_to_seg(original_cell,
                      section_per_subtree_index,
                      roots_of_subtrees,
                      mapping_sections_to_subtree_index,
                      new_cable_properties,
                      has_apical,
                      apic,
                      basals,
                      subtree_ind_to_q,
                      mapping_type,
                      reduction_frequency):
    '''create mapping between segments in the original model to segments in the reduced model

       if mapping_type == impedance the mapping will be a response to the
       transfer impedance of each segment to the soma (like the synapses)

       if mapping_type == distance  the mapping will be a response to the
       distance of each segment to the soma (like the synapses) NOT IMPLEMENTED
       YET

       '''

    assert mapping_type == 'impedance', 'distance mapping not implemented yet'
    # the keys are the segments of the original model, the values are the
    # segments of the reduced model
    original_seg_to_reduced_seg = {}
    reduced_seg_to_original_seg = collections.defaultdict(list)
    for subtree_index in section_per_subtree_index:
        for sec in section_per_subtree_index[subtree_index]:
            for seg in sec:
                synapse_location = find_synapse_loc(seg, mapping_sections_to_subtree_index)
                imp_obj, subtree_input_impedance = measure_input_impedance_of_subtree(
                    roots_of_subtrees[subtree_index], reduction_frequency)

                # if synapse is on the apical subtree
                on_basal_subtree = not (has_apical and subtree_index == 0)

                mid_of_segment_loc = reduce_synapse(
                    original_cell,
                    synapse_location,
                    on_basal_subtree,
                    imp_obj,
                    subtree_input_impedance,
                    new_cable_properties[subtree_index].electrotonic_length,
                    subtree_ind_to_q[subtree_index])

                if on_basal_subtree:
                    if has_apical:
                        new_section_for_synapse = basals[subtree_index - 1]
                    else:
                        new_section_for_synapse = basals[subtree_index]
                else:
                    new_section_for_synapse = apic

                reduced_seg = new_section_for_synapse(mid_of_segment_loc)
                original_seg_to_reduced_seg[seg] = reduced_seg
                reduced_seg_to_original_seg[reduced_seg].append(seg)

    return original_seg_to_reduced_seg, dict(reduced_seg_to_original_seg)


def copy_dendritic_mech(original_seg_to_reduced_seg,
                        reduced_seg_to_original_seg,
                        apic,
                        basals,
                        segment_to_mech_vals,
                        mapping_type='impedance'):
    ''' copies the mechanisms from the original model to the reduced model'''

    # copy mechanisms
    # this is needed for the case where some segements were not been mapped
    mech_names_per_segment = collections.defaultdict(list)
    vals_per_mech_per_segment = {}
    for reduced_seg, original_segs in reduced_seg_to_original_seg.items():
        vals_per_mech_per_segment[reduced_seg] = collections.defaultdict(list)

        for original_seg in original_segs:
            for mech_name, mech_params in segment_to_mech_vals[original_seg].items():
                for param_name, param_value in mech_params.items():
                    vals_per_mech_per_segment[reduced_seg][param_name].append(param_value)

                mech_names_per_segment[reduced_seg].append(mech_name)
                reduced_seg.sec.insert(mech_name)

        for param_name, param_values in vals_per_mech_per_segment[reduced_seg].items():
            setattr(reduced_seg, param_name, np.mean(param_values))

    all_segments = []
    if apic is not None:
        all_segments.extend(list(apic))

    for bas in basals:
        all_segments.extend(list(bas))

    if len(all_segments) != len(reduced_seg_to_original_seg):
        logger.warning('There is no segment to segment copy, it means that some segments in the'
                    'reduced model did not receive channels from the original cell.'
                    'Trying to compensate by copying channels from neighboring segments')
        handle_orphan_segments(original_seg_to_reduced_seg,
                               all_segments,
                               vals_per_mech_per_segment,
                               mech_names_per_segment)


def handle_orphan_segments(original_seg_to_reduced_seg,
                           all_segments,
                           vals_per_mech_per_segment,
                           mech_names_per_segment):
    ''' This function handle reduced segments that did not had original segments mapped to them'''
    # Get all reduced segments that have been mapped by a original model segment
    all_mapped_control_segments = original_seg_to_reduced_seg.values()
    non_mapped_segments = set(all_segments) - set(all_mapped_control_segments)

    for reduced_seg in non_mapped_segments:
        seg_secs = list(reduced_seg.sec)
        # find valid parent
        parent_seg_index = seg_secs.index(reduced_seg) - 1
        parent_seg = None
        while parent_seg_index > -1:
            if seg_secs[parent_seg_index] in all_mapped_control_segments:
                parent_seg = seg_secs[parent_seg_index]
                break
            else:
                parent_seg_index -= 1

        # find valid child
        child_seg_index = seg_secs.index(reduced_seg) + 1
        child_seg = None
        while child_seg_index < len(seg_secs):
            if seg_secs[child_seg_index] in all_mapped_control_segments:
                child_seg = seg_secs[child_seg_index]
                break
            else:
                child_seg_index += 1

        if not parent_seg and not child_seg:
            raise Exception("no child seg nor parent seg, with active channels, was found")

        if parent_seg and not child_seg:
            for mech in mech_names_per_segment[parent_seg]:
                reduced_seg.sec.insert(mech)
            for n in vals_per_mech_per_segment[parent_seg]:
                setattr(reduced_seg, n, np.mean(vals_per_mech_per_segment[parent_seg][n]))

        if not parent_seg and child_seg:
            for mech in mech_names_per_segment[child_seg]:
                reduced_seg.sec.insert(mech)
            for n in vals_per_mech_per_segment[child_seg]:
                setattr(reduced_seg, n, np.mean(vals_per_mech_per_segment[child_seg][n]))

        # if both parent and child were found, we add to the segment all the mech in both
        # this is just a decision

        if parent_seg and child_seg:
            for mech in set(mech_names_per_segment[child_seg]) & set(mech_names_per_segment[parent_seg]):
                reduced_seg.sec.insert(mech)

            for n in vals_per_mech_per_segment[child_seg]:
                child_mean = np.mean(vals_per_mech_per_segment[child_seg][n])
                if n in vals_per_mech_per_segment[parent_seg]:
                    parent_mean = np.mean(vals_per_mech_per_segment[parent_seg][n])
                    setattr(reduced_seg, n, (child_mean + parent_mean) / 2)
                else:
                    setattr(reduced_seg, n, child_mean)

            for n in vals_per_mech_per_segment[parent_seg]:
                parent_mean = np.mean(vals_per_mech_per_segment[parent_seg][n])
                if n in vals_per_mech_per_segment[child_seg]:
                    child_mean = np.mean(vals_per_mech_per_segment[child_seg][n])
                    setattr(reduced_seg, n, (child_mean + parent_mean) / 2)
                else:
                    setattr(reduced_seg, n, parent_mean)


def add_PP_properties_to_dict(PP, PP_params_dict):
    '''
    add the propeties of a point process to PP_params_dict.
    The only propeties added to the dictionary are those worth comparing
    '''
    skipped_params = {"Section", "allsec", "baseattr", "cas", "g", "get_loc", "has_loc", "hname",
                      'hocobjptr', "i", "loc", "next", "ref", "same", "setpointer", "state",
                      "get_segment",
                      }
    PP_params = []
    for param in dir(PP):
        if param.startswith("__") or param in skipped_params:
            continue
        PP_params.append(param)
    PP_params_dict[type_of_point_process(PP)] = PP_params


def type_of_point_process(PP):
    s = PP.hname()
    ix = PP.hname().find("[")
    return s[:ix]


def apply_params_to_section(name, type_of_sectionlist, instance_as_str, section, cable_params, nseg):
    section.L = cable_params.length
    section.diam = cable_params.diam
    section.nseg = nseg

    append_to_section_lists(name, type_of_sectionlist, instance_as_str)

    section.insert('pas')
    section.cm = cable_params.cm
    section.g_pas = 1.0 / cable_params.rm
    section.Ra = cable_params.ra
    section.e_pas = cable_params.e_pas


def calculate_subtree_q(root, reduction_frequency):
    rm = 1.0 / root.g_pas
    rc = rm * (float(root.cm) / 1000000)
    angular_freq = 2 * math.pi * reduction_frequency
    q_imaginary = angular_freq * rc
    q_subtree = complex(1, q_imaginary)   # q=1+iwRC
    q_subtree = cmath.sqrt(q_subtree)
    return q_subtree


def synapse_properties_match(synapse, PP, PP_params_dict):
    if PP.hname()[:PP.hname().rindex('[')] != synapse.hname()[:synapse.hname().rindex('[')]:
        return False
    for param in PP_params_dict[type_of_point_process(PP)]:
        if(param not in ['rng'] and  # https://github.com/neuronsimulator/nrn/issues/136
           str(type(getattr(PP, param))) != "<type 'hoc.HocObject'>" and  # ignore hoc objects
           getattr(PP, param) != getattr(synapse, param)):
            return False
    return True


def load_model(model_filename):
    model_obj_name = model_filename.split(".")[0].split('/')[-1]
    if h.name_declared(model_obj_name) == 0:
        logger.debug("loading template '%s'" % model_obj_name)
        if model_filename == 'model.hoc':
            logger.debug("loading default reduced model")
            load_default_model()
        else:
            h.load_file(model_filename)
    else:
        logger.info("The template '%s' is already defined... not loading." % model_obj_name)
    return model_obj_name


def gather_subtrees(soma_ref):
    '''get all the subtrees of the soma

    assumes the axon is already disconnected
    return (list(roots_of_subtrees), list(num_of_subtrees))
    where:
      roots_of_subtrees holds the root sections of each of the soma's subtrees
        note: The apical, if it exists, has been moved to the front
      num_of_subtrees correctly the number of subtrees, excluding the axon
    '''

    roots_of_subtrees = []
    num_of_subtrees = []
    for i in range(int(soma_ref.nchild())):
        if 'soma' in str(soma_ref.child[i]):
            logger.warning("soma is child, ignore - not tested yet")
            continue
        num_of_subtrees.append(i)
        roots_of_subtrees.append(soma_ref.child[i])

    # assuming up to one apical tree
    ix_of_apical = None
    for i in num_of_subtrees:
        if 'apic' in roots_of_subtrees[i].hname():
            assert ix_of_apical is None, 'Multiple apical dendrites not suppored'
            ix_of_apical = i

    if ix_of_apical is not None:
        roots_of_subtrees = ([roots_of_subtrees[ix_of_apical]] +
                             roots_of_subtrees[:ix_of_apical] +
                             roots_of_subtrees[ix_of_apical + 1:])
    return roots_of_subtrees, num_of_subtrees


def gather_cell_subtrees(roots_of_subtrees):
    # dict that maps section indexes to the subtree index they are in: keys are
    # string tuples: ("apic"/"basal", orig_section_index) , values are ints:
    # subtree_instance_index
    sections_to_delete = []
    section_per_subtree_index = {}
    mapping_sections_to_subtree_index = {}
    for i, soma_child in enumerate(roots_of_subtrees):
        # inserts each section in this subtree into the above dict, which maps
        # it to the subtree index
        if 'apic' in soma_child.hname():
            assert i == 0, ('The apical is not the first child of the soma! '
                            'a code refactoring is needed in order to accept it')
            mark_subtree_sections_with_subtree_index(sections_to_delete,
                                                     section_per_subtree_index,
                                                     soma_child,
                                                     mapping_sections_to_subtree_index,
                                                     "apic",
                                                     i)
        elif 'dend' in soma_child.hname() or 'basal' in soma_child.hname():
            mark_subtree_sections_with_subtree_index(sections_to_delete,
                                                     section_per_subtree_index,
                                                     soma_child,
                                                     mapping_sections_to_subtree_index,
                                                     "basal",
                                                     i)

    return sections_to_delete, section_per_subtree_index, mapping_sections_to_subtree_index


def create_reduced_cell(soma_cable,
                        has_apical,
                        original_cell,
                        model_obj_name,
                        new_cable_properties,
                        new_cables_nsegs,
                        subtrees_xs):
    h("objref reduced_cell")
    h("reduced_cell = new " + model_obj_name + "()")

    create_sections_in_hoc("soma", 1, "reduced_cell")

    soma = original_cell.soma[0] if original_cell.soma.hname()[-1] == ']' else original_cell.soma
    append_to_section_lists("soma[0]", "somatic", "reduced_cell")

    if has_apical:  # creates reduced apical cable if apical subtree existed
        create_sections_in_hoc("apic", 1, "reduced_cell")
        apic = h.reduced_cell.apic[0]
        num_of_basal_subtrees = len(new_cable_properties) - 1

        cable_params = new_cable_properties[0]
        nseg = new_cables_nsegs[0]
        apply_params_to_section("apic[0]", "apical", "reduced_cell",
                                apic, cable_params, nseg)
        apic.connect(soma, subtrees_xs[0], 0)
    else:
        apic = None
        num_of_basal_subtrees = len(new_cable_properties)

    # creates reduced basal cables
    create_sections_in_hoc("dend", num_of_basal_subtrees, "reduced_cell")
    basals = [h.reduced_cell.dend[i] for i in range(num_of_basal_subtrees)]

    for i in range(num_of_basal_subtrees):
        if has_apical:
            index_in_reduced_cables_dimensions = i + 1
        else:
            index_in_reduced_cables_dimensions = i

        cable_params = new_cable_properties[index_in_reduced_cables_dimensions]
        nseg = new_cables_nsegs[index_in_reduced_cables_dimensions]

        apply_params_to_section("dend[" + str(i) + "]", "basal", "reduced_cell",
                                basals[i], cable_params, nseg)

        basals[i].connect(soma, subtrees_xs[index_in_reduced_cables_dimensions], 0)

    # create cell python template
    cell = Neuron(h.reduced_cell)
    cell.soma = original_cell.soma
    cell.apic = apic

    return cell, basals


def merge_and_add_synapses(num_of_subtrees,
                           new_cable_properties,
                           PP_params_dict,
                           synapses_list,
                           mapping_sections_to_subtree_index,
                           netcons_list,
                           has_apical,
                           roots_of_subtrees,
                           original_cell,
                           basals,
                           cell,
                           reduction_frequency):
    # dividing the original synapses into baskets, so that all synapses that are
    # on the same subtree will be together in the same basket

    # a list of baskets of synapses, each basket in the list will hold the
    # synapses of the subtree of the corresponding basket index
    baskets = [[] for _ in num_of_subtrees]
    soma_synapses_syn_to_netcon = {}

    for syn_index, synapse in enumerate(synapses_list):
        synapse_location = find_synapse_loc(synapse, mapping_sections_to_subtree_index)

        # for a somatic synapse
        # TODO: 'axon' is never returned by find_synapse_loc...
        if synapse_location.subtree_index in (SOMA_LABEL, 'axon'):
            soma_synapses_syn_to_netcon[synapse] = netcons_list[syn_index]
        else:
            baskets[synapse_location.subtree_index].append((synapse, synapse_location, syn_index))

    # mapping (non-somatic) synapses to their new location on the reduced model
    # (the new location is the exact location of the middle of the segment they
    # were mapped to, in order to enable merging)
    new_synapses_list, subtree_ind_to_q = [], {}
    for subtree_index in num_of_subtrees:
        imp_obj, subtree_input_impedance = measure_input_impedance_of_subtree(
            roots_of_subtrees[subtree_index], reduction_frequency)
        subtree_ind_to_q[subtree_index] = calculate_subtree_q(
            roots_of_subtrees[subtree_index], reduction_frequency)

        # iterates over the synapses in the curr basket
        for synapse, synapse_location, syn_index in baskets[subtree_index]:
            on_basal_subtree = not (has_apical and subtree_index == 0)

            # "reduces" the synapse - finds this synapse's new "merged"
            # location on its corresponding reduced cable
            x = reduce_synapse(original_cell,
                               synapse_location,
                               on_basal_subtree,
                               imp_obj,
                               subtree_input_impedance,
                               new_cable_properties[subtree_index].electrotonic_length,
                               subtree_ind_to_q[subtree_index])

            # find the section of the synapse
            if on_basal_subtree:
                if has_apical:
                    section_for_synapse = basals[subtree_index - 1]
                else:
                    section_for_synapse = basals[subtree_index]
            else:
                section_for_synapse = cell.apic

            # go over all point processes in this segment and see whether one
            # of them has the same proporties of this synapse
            # If there's such a synapse link the original NetCon with this point processes
            # If not, move the synapse to this segment.
            for PP in section_for_synapse(x).point_processes():
                if type_of_point_process(PP) not in PP_params_dict:
                    add_PP_properties_to_dict(PP, PP_params_dict)

                if synapse_properties_match(synapse, PP, PP_params_dict):
                    netcons_list[syn_index].setpost(PP)
                    break
            else:  # If for finish the loop -> first appearance of this synapse
                synapse.loc(x, sec=section_for_synapse)
                new_synapses_list.append(synapse)

    # merging somatic and axonal synapses
    synapses_per_seg = collections.defaultdict(list)
    for synapse in soma_synapses_syn_to_netcon:
        seg_pointer = synapse.get_segment()

        for PP in synapses_per_seg[seg_pointer]:
            if type_of_point_process(PP) not in PP_params_dict:
                add_PP_properties_to_dict(PP, PP_params_dict)

            if synapse_properties_match(synapse, PP, PP_params_dict):
                soma_synapses_syn_to_netcon[synapse].setpost(PP)
                break
        else:  # If for finish the loop -> first appearance of this synapse
            synapse.loc(seg_pointer.x, sec=seg_pointer.sec)
            new_synapses_list.append(synapse)
            synapses_per_seg[seg_pointer].append(synapse)

    return new_synapses_list, subtree_ind_to_q

def textify_seg_to_seg(segs):
    '''convert segment dictionary to text'''
    ret = {str(k): str(v) for k, v in segs.items()}
    return ret
   
def subtree_reductor(original_cell,
                     synapses_list,
                     netcons_list,
                     reduction_frequency,
                     model_filename='model.hoc',
                     total_segments_manual=-1,
                     PP_params_dict=None,
                     mapping_type='impedance',
                     return_seg_to_seg=False
                     ):

    '''
    Receives an instance of a cell with a loaded full morphology, a list of
    synapse objects, a list of NetCon objects (the i'th netcon in the list
    should correspond to the i'th synapse), the filename (string) of the model
    template hoc file that the cell was instantiated from, the desired
    reduction frequency as a float, optional parameter for the approximate
    desired number of segments in the new model (if this parameter is empty,
    the number of segments will be such that there is a segment for every 0.1
    lambda), and an optional param for the point process to be compared before
    deciding on whether to merge a synapse or not and reduces the cell (using
    the given reduction_frequency). Creates a reduced instance using the model
    template in the file whose filename is given as a parameter, and merges
    synapses of the same type that get mapped to the same segment
    (same "reduced" synapse object for them all, but different NetCon objects).



    model_filename : model.hoc  will use a default template
    total_segments_manual: sets the number of segments in the reduced model
                           can be either -1, a float between 0 to 1, or an int
                           if total_segments_manual = -1 will do automatic segmentation
                           if total_segments_manual>1 will set the number of segments
                           in the reduced model to total_segments_manual
                           if 0>total_segments_manual>1 will automatically segment the model
                           but if the automatic segmentation will produce a segment number that
                           is lower than original_number_of_segments*total_segments_manual it
                           will set the number of segments in the reduced model to:
                           original_number_of_segments*total_segments_manual
    return_seg_to_seg: if True the function will also return a textify version of the mapping
                       between the original segments to the reduced segments 


    Returns the new reduced cell, a list of the new synapses, and the list of
    the inputted netcons which now have connections with the new synapses.

    Notes:
    1) The original cell instance, synapses and Netcons given as arguments are altered
    by the function and cannot be used outside of it in their original context.
    2) Synapses are determined to be of the same type and mergeable if their reverse
    potential, tau1 and tau2 values are identical.
    3) Merged synapses are assigned a single new synapse object that represents them
    all, but keep their original NetCon objects. Each such NetCon now connects the
    original synapse's NetStim with
    the reduced synapse.
    '''
    if PP_params_dict is None:
        PP_params_dict = {}

    h.init()

    model_obj_name = load_model(model_filename)

    # finds soma properties
    soma = original_cell.soma[0] if original_cell.soma.hname()[-1] == ']' else original_cell.soma

    soma_cable = CableParams(length=soma.L, diam=soma.diam, space_const=None,
                             cm=soma.cm, rm=1.0 / soma.g_pas, ra=soma.Ra, e_pas=soma.e_pas,
                             electrotonic_length=None)

    has_apical = len(list(original_cell.apical)) != 0

    soma_ref = h.SectionRef(sec=soma)
    axon_section, axon_is_parent, soma_axon_x = find_and_disconnect_axon(soma_ref)
    roots_of_subtrees, num_of_subtrees = gather_subtrees(soma_ref)

    sections_to_delete, section_per_subtree_index, mapping_sections_to_subtree_index = \
        gather_cell_subtrees(roots_of_subtrees)

    # preparing for reduction

    # remove active conductances and get seg_to_mech dictionary
    segment_to_mech_vals = create_segments_to_mech_vals(sections_to_delete)

    # disconnects all the subtrees from the soma
    subtrees_xs = []
    for subtree_root in roots_of_subtrees:
        subtrees_xs.append(subtree_root.parentseg().x)
        h.disconnect(sec=subtree_root)

    # reducing the subtrees
    new_cable_properties = [reduce_subtree(roots_of_subtrees[i], reduction_frequency)
                            for i in num_of_subtrees]

    if total_segments_manual > 1:
        new_cables_nsegs = calculate_nsegs_from_manual_arg(new_cable_properties,
                                                           total_segments_manual)
    else:
        new_cables_nsegs = calculate_nsegs_from_lambda(new_cable_properties)
        if total_segments_manual > 0:
            original_cell_seg_n = (sum(i.nseg for i in list(original_cell.basal)) +
                                   sum(i.nseg for i in list(original_cell.apical))
                                   )
            min_reduced_seg_n = int(round((total_segments_manual * original_cell_seg_n)))
            if sum(new_cables_nsegs) < min_reduced_seg_n:
                logger.debug("number of segments calculated using lambda is {}, "
                      "the original cell had {} segments.  "
                      "The min reduced segments is set to {}% of reduced cell segments".format(
                          sum(new_cables_nsegs),
                          original_cell_seg_n,
                          total_segments_manual * 100))
                logger.debug("the reduced cell nseg is set to %s" % min_reduced_seg_n)
                new_cables_nsegs = calculate_nsegs_from_manual_arg(new_cable_properties,
                                                                   min_reduced_seg_n)

    cell, basals = create_reduced_cell(soma_cable,
                                       has_apical,
                                       original_cell,
                                       model_obj_name,
                                       new_cable_properties,
                                       new_cables_nsegs,
                                       subtrees_xs)

    new_synapses_list, subtree_ind_to_q = merge_and_add_synapses(
        num_of_subtrees,
        new_cable_properties,
        PP_params_dict,
        synapses_list,
        mapping_sections_to_subtree_index,
        netcons_list,
        has_apical,
        roots_of_subtrees,
        original_cell,
        basals,
        cell,
        reduction_frequency)

    # create segment to segment mapping
    original_seg_to_reduced_seg, reduced_seg_to_original_seg = create_seg_to_seg(
        original_cell,
        section_per_subtree_index,
        roots_of_subtrees,
        mapping_sections_to_subtree_index,
        new_cable_properties,
        has_apical,
        cell.apic,
        basals,
        subtree_ind_to_q,
        mapping_type,
        reduction_frequency)

    # copy active mechanisms
    copy_dendritic_mech(original_seg_to_reduced_seg,
                        reduced_seg_to_original_seg,
                        cell.apic,
                        basals,
                        segment_to_mech_vals,
                        mapping_type)
    
    if return_seg_to_seg:
        original_seg_to_reduced_seg_text = textify_seg_to_seg(original_seg_to_reduced_seg)

    # Connect axon back to the soma
    if len(axon_section) > 0:
        if axon_is_parent:
            soma.connect(axon_section[0])
        else:
            axon_section[0].connect(soma, soma_axon_x)

    # Now we delete the original model
    for section in sections_to_delete:
        with push_section(section):
            h.delete_section()

    cell.axon = axon_section
    cell.dend = cell.hoc_model.dend

    with push_section(cell.hoc_model.soma[0]):
        h.delete_section()
    if return_seg_to_seg:
        return cell, new_synapses_list, netcons_list, original_seg_to_reduced_seg_text
    else:
        return cell, new_synapses_list, netcons_list


class Neuron(object):
    'Python neuron class for hoc models'
    def __init__(self, model):
        self.hoc_model = model
        self.soma = None
        self.dend = None
        self.apic = None
        self.axon = None


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
