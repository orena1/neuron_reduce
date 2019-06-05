#!/usr/bin/env python
import os

WRITE_UNIT_TEST_VECTORS = False
PLOT_VOLTAGES = False

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
TESTDATA_PATH = os.path.join(BASE_PATH, 'TestsFiles')


def run_reduce(morphology_file,
               model_file,
               frequency,
               synapse_file,
               voltage_file,
               create_type,
               celsius,
               write_unit_test_vectors=WRITE_UNIT_TEST_VECTORS,
               plot_voltages=PLOT_VOLTAGES,
               reduced_model_file='model.hoc',
               manual_total_nsegs=-1):
    args = ("python " +
            os.path.join(BASE_PATH, "test_script_helper.py ") +
            ' '.join([str(a) for a in (morphology_file,
                                       model_file,
                                       reduced_model_file,
                                       frequency,
                                       manual_total_nsegs,
                                       synapse_file,
                                       voltage_file,
                                       write_unit_test_vectors,
                                       plot_voltages,
                                       create_type,
                                       celsius)]))
    assert os.system(args) == 0


def test1():
    '''Test 1 passive neuron'''
    path = os.path.join(TESTDATA_PATH, 'Test_1')
    kwargs = dict(morphology_file=os.path.join(path, "2013_03_06_cell08_876_H41_05_Cell2.ASC"),
                  model_file=os.path.join(path, "model.hoc"),
                  synapse_file=os.path.join(path, "origRandomSynapses-10000"),
                  create_type='basic',
                  celsius=37)

    for frequency in (0, 10, 38, 200):
        kwargs['frequency'] = frequency
        kwargs['voltage_file'] = os.path.join(path, "voltage_vectors_for_unit_test_%s.txt" % frequency)
        run_reduce(**kwargs)


def test2():
    '''Test 2 passive neuron not deleting the axon'''
    path = os.path.join(TESTDATA_PATH, 'Test_2')

    run_reduce(morphology_file=os.path.join(path, "2013_03_06_cell08_876_H41_05_Cell2.ASC"),
               model_file=os.path.join(path, "model.hoc"),
               frequency=38,
               synapse_file=os.path.join(path, "origRandomSynapses-10000"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='basic',
               celsius=37,
               )


def test3_a71075_passive():
    path = os.path.join(TESTDATA_PATH, 'Test_3')
    run_reduce(morphology_file=os.path.join(path, "dend-C050800E2_cor_axon-C120398A-P2_-_Scale_x1.000_y1.050_z1.000_-_Clone_81.asc"),
               model_file=os.path.join(path, "cADpyr230_L4_SS_4_dend_C050800E2_cor_axon_C120398A_P2___Scale_x1_000_y1_050_z1_000___Clone_81.hoc"),
               frequency=38,
               synapse_file=os.path.join(path, "synapse_fromh5a71075.txt"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='bbp',
               celsius=34,
               )


def test4_amsalem_2016():
    path = os.path.join(TESTDATA_PATH, 'Test_4_LBC_amsalem/')
    run_reduce(morphology_file=os.path.join(path, "C230300D1.asc"),
               model_file=os.path.join(path, "cNAC187_L23_LBC_3_C230300D1_new_new_fit.hoc"),
               frequency=9,
               synapse_file=os.path.join(path, "synapse_fromh5a71075.txt"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='bbpactive',
               celsius=34
               )


def test5_Hay_2011_active_dendrite():
    path = os.path.join(TESTDATA_PATH, 'Test_5_Hay_2011/')
    run_reduce(morphology_file=os.path.join(path, "cell1.asc"),
               model_file=os.path.join(path, "L5PCtemplate.hoc"),
               frequency=38,
               synapse_file=os.path.join(path, "origRandomSynapses-10000"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='hay',
               celsius=37)


def test6_L4_LBC_cNAC187_5_for_run():
    path = os.path.join(TESTDATA_PATH, 'L4_LBC_cNAC187_5_for_run/')
    run_reduce(morphology_file=os.path.join(path, "2013_03_06_cell08_876_H41_05_Cell2.ASC"),
               model_file=os.path.join(path, "cNAC187_L4_LBC_8e834c24cb.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "1487081844_732516.txt"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test_0.txt"),
               create_type='bbpnew',
               celsius=34,
               )


def test7_Almog_Korngreen_2014():
    path = os.path.join(TESTDATA_PATH, 'Test_7_Almog/')
    run_reduce(morphology_file=path,
               model_file=os.path.join(path, "A140612_1.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "origRandomSynapses-10000"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test_0.txt"),
               create_type='almog',
               celsius=34,
               )


def test8_Marasco_Limongiello_Migliore_2012():
    '''Test 8 Marasco Limongiello Migliore 2012'''
    path = os.path.join(TESTDATA_PATH, 'Test_8_C1_Marasco/')
    run_reduce(morphology_file=path,
               model_file=os.path.join(path, "geo5038801modMod.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "1487081844_732516.txt"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test_0.txt"),
               create_type='almog',
               celsius=34,
               )


def test9_model_48310820():
    '''Test 9 model 48310820 (L5PC) from the Allen celltypes data base'''
    path = os.path.join(TESTDATA_PATH, 'Test_9_Allen_483108201/')
    run_reduce(morphology_file=os.path.join(path, "reconstruction.swc"),
               model_file=os.path.join(path, "AllenTemplate.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "origRandomSynapses-10000"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='allen',
               celsius=34,
               )


def test10_model_47804508():
    '''Test 10 model 47804508 (L1) from the Allen celltypes data base'''
    path = os.path.join(TESTDATA_PATH, 'Test_10_Allen_47804508/')
    run_reduce(morphology_file=os.path.join(path, "reconstruction.swc"),
               model_file=os.path.join(path, "AllenTemplate.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "origRandomSynapses-10000"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='allen',
               celsius=34,
               )


def test11_human_Eyal_2016():
    '''Test 11 model Human L2/3 Cell from Eyal et al 2016'''
    path = os.path.join(TESTDATA_PATH, 'Test_11_Human_L2_3_Eyal/')
    run_reduce(morphology_file=os.path.join(path, "2013_03_06_cell08_876_H41_05_Cell2.ASC"),
               model_file=os.path.join(path, "model_0603_cell08_cm045.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "origRandomSynapses-10000"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='human',
               celsius=37,
               )


def test12_TPC_Markram_2016():
    '''Test 12 Tufted Pyramidal Cell (L6) Markram et al. Cell (2015) ----'''
    path = os.path.join(TESTDATA_PATH, 'Test_12_TPC_L6_Markram/')
    run_reduce(morphology_file=os.path.join(path, "dend-tkb070125a3_ch1_cc2_b_hw_60x_1_axon-tkb060223b3_ch1_cc2_o_ps_60x_1_-_Clone_5.asc"),
               model_file=os.path.join(path, "cADpyr231_L6_TPC_L1_44f2206f70.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "synapses_location.txt"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type='bbpnew',
               celsius=34,
               )


def test13_dbc_Markram_2015():
    '''Test 13 Double Bouquet Cell (L4) Markram et al. Cell (2015)'''
    path = os.path.join(TESTDATA_PATH, 'Test_13_DBC_L4_Markram/')
    run_reduce(morphology_file=os.path.join(path, "C140600C-I1_-_Clone_2.asc"),
               model_file=os.path.join(path, "cNAC187_L4_DBC_23ffe29c8b.hoc"),
               frequency=0,
               synapse_file=os.path.join(path, "synapses_locations.txt"),
               voltage_file=os.path.join(path, "voltage_vectors_for_unit_test.txt"),
               create_type ='bbpnew',
               celsius=34,
               )


if __name__ == '__main__':
    test1()
    test2()
    test3_a71075_passive()
    test4_amsalem_2016()
    test5_Hay_2011_active_dendrite()
    test6_L4_LBC_cNAC187_5_for_run()
    test7_Almog_Korngreen_2014()
    test8_Marasco_Limongiello_Migliore_2012()
    test9_model_48310820()
    test10_model_47804508()
    test11_human_Eyal_2016()
    test12_TPC_Markram_2016()
    test13_dbc_Markram_2015()
