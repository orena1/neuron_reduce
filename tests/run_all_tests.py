import os
import os


WriteVectors = False
plot_voltages = True
def Createtxt():
    return("python test_script_helper.py " +  orig_morphology_file + ' ' + orig_model_file  + ' ' +  reduced_model_file  + 
           ' ' +  str(frequency)  + ' ' +  str(manual_total_nsegs)  + ' ' +  synapse_file  + ' ' +  voltage_file  + ' ' +  
           str(write_unit_test_vectors)  + ' ' +  str(plot_voltages) + ' ' + create_type + ' ' + str(celsius))

# ### Test 1 passive neuron
Path = 'TestsFiles/Test_1/'
orig_morphology_file = Path + "2013_03_06_cell08_876_H41_05_Cell2.ASC"
orig_model_file = Path + "model.hoc"
reduced_model_file = "model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test_0.txt"
write_unit_test_vectors = WriteVectors
create_type = 'basic'
celsius = 37

os.system(Createtxt())
print("Done Test 1 freq 0 \n")

frequency = 10
voltage_file = Path + "voltage_vectors_for_unit_test_" + `frequency` +".txt"


os.system(Createtxt())
print("Done Test 1 freq 10 \n")

frequency = 38
voltage_file = Path + "voltage_vectors_for_unit_test_" + `frequency` +".txt"

os.system(Createtxt())
print("Done Test 1 freq 38 \n")

frequency = 200
voltage_file = Path + "voltage_vectors_for_unit_test_" + `frequency` +".txt"

os.system(Createtxt())
print("Done Test 1 freq 200 \n")

### Test 2 passive neuron not deleting the axon ###
Path = 'TestsFiles/Test_2/'
orig_morphology_file = Path + "2013_03_06_cell08_876_H41_05_Cell2.ASC"
orig_model_file = Path + "model.hoc"
reduced_model_file ="model.hoc"
frequency = 38
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'basic'
celsius = 37
os.system(Createtxt())
print("Done Test 2\n")
   
    
### Test 3 a71075 passive ###
Path = 'TestsFiles/Test_3/'
orig_morphology_file = Path + "dend-C050800E2_cor_axon-C120398A-P2_-_Scale_x1.000_y1.050_z1.000_-_Clone_81.asc"
orig_model_file = Path + "cADpyr230_L4_SS_4_dend_C050800E2_cor_axon_C120398A_P2___Scale_x1_000_y1_050_z1_000___Clone_81.hoc"
reduced_model_file ="model.hoc"
frequency = 38
manual_total_nsegs = -1
synapse_file = Path + "synapse_fromh5a71075.txt"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'bbp'
celsius = 34
os.system(Createtxt())
print("Done Test 3 a71075 passive\n")


### Test 4 amsalem et al 2016  ### *
Path = 'TestsFiles/Test_4_LBC_amsalem/'
orig_morphology_file = Path + "C230300D1.asc"
orig_model_file = Path + "cNAC187_L23_LBC_3_C230300D1_new_new_fit.hoc"
reduced_model_file ="model.hoc"
frequency = 9
manual_total_nsegs = -1
synapse_file = Path + "synapse_fromh5a71075.txt"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'bbpactive'
celsius = 34
os.system(Createtxt())
print("Done Test 3 a71075 passive\n")



### Test 5 Hay et al 2011 active dendrite ###
Path = 'TestsFiles/Test_5_Hay_2011/'
orig_morphology_file = Path + "cell1.asc"
orig_model_file = Path + "L5PCtemplate.hoc"
reduced_model_file ="model.hoc"
frequency = 38
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'hay'
celsius = 37
os.system(Createtxt())
print("Done Test 5 Hay Model\n")

## Test 6 Path = L4_LBC_cNAC187_5_for_run/'
Path = 'TestsFiles/L4_LBC_cNAC187_5_for_run/'
orig_morphology_file = Path + "2013_03_06_cell08_876_H41_05_Cell2.ASC"
orig_model_file = Path + "cNAC187_L4_LBC_8e834c24cb.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "1487081844_732516.txt"
voltage_file = Path + "voltage_vectors_for_unit_test_0.txt"
write_unit_test_vectors = WriteVectors
create_type = 'bbpnew'
celsius = 34
os.system(Createtxt())
print("Done Test 6 L4_LBC\n")
#asda

## Test 7 Almog and Korngreen 2014'
Path = 'TestsFiles/Test_7_Almog/'
orig_morphology_file = Path + ""
orig_model_file = Path + "A140612_1.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test_0.txt"
write_unit_test_vectors = WriteVectors
create_type = 'almog'
celsius = 34
os.system(Createtxt())
print("Done Test 7, Almog and Korngreen 2014\n")

# Test 8 Marasco Limongiello Migliore 2012'
Path = 'TestsFiles/Test_8_C1_Marasco/'
orig_morphology_file = Path + ""
orig_model_file = Path + "geo5038801modMod.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "1487081844_732516.txt"
voltage_file = Path + "voltage_vectors_for_unit_test_0.txt"
write_unit_test_vectors = WriteVectors
create_type = 'almog'
celsius = 34
os.system(Createtxt())
print("Done Test 8 Marasco Limongiello Migliore 2012\n")

### Test 9 model 48310820 (L5PC) from the Allen celltypes data base
Path = 'TestsFiles/Test_9_Allen_483108201/'
orig_morphology_file = Path + "reconstruction.swc"
orig_model_file = Path + "AllenTemplate.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'allen'
celsius = 34
print Createtxt()
os.system(Createtxt())
print("Done Test 9, Allen V1 L5PC\n")

### Test 10 model 47804508 (L1) from the Allen celltypes data base
Path = 'TestsFiles/Test_10_Allen_47804508/'
orig_morphology_file = Path + "reconstruction.swc"
orig_model_file = Path + "AllenTemplate.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'allen'
celsius = 34
print Createtxt()
os.system(Createtxt())
print("Done Test 10, Allen V1 L1\n")

### Test 11 model Human L2/3 Cell from Eyal et al 2016
Path = 'TestsFiles/Test_11_Human_L2_3_Eyal/'
orig_morphology_file = Path + "2013_03_06_cell08_876_H41_05_Cell2.ASC"
orig_model_file = Path + "model_0603_cell08_cm045.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "origRandomSynapses-10000"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'human'
celsius = 37
print Createtxt()
os.system(Createtxt())
print("Done Test 11, Human L2/3\n")

### Test 12 Tufted Pyramidal Cell (L6) Markram et al. Cell (2015) ---- 
Path = 'TestsFiles/Test_12_TPC_L6_Markram/'
orig_morphology_file = Path + "dend-tkb070125a3_ch1_cc2_b_hw_60x_1_axon-tkb060223b3_ch1_cc2_o_ps_60x_1_-_Clone_5.asc"
orig_model_file = Path + "cADpyr231_L6_TPC_L1_44f2206f70.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "synapses_location.txt"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type = 'bbpnew'
celsius = 34
print Createtxt()
os.system(Createtxt())
print("Done Test 12, Tufted Pyramidal Cell (L6) Markram et al. Cell (2015)\n")

### Test 13 Double Bouquet Cell (L4) Markram et al. Cell (2015) ---- 
Path = 'TestsFiles/Test_13_DBC_L4_Markram/'
orig_morphology_file = Path + "C140600C-I1_-_Clone_2.asc"
orig_model_file = Path + "cNAC187_L4_DBC_23ffe29c8b.hoc"
reduced_model_file ="model.hoc"
frequency = 0
manual_total_nsegs = -1
synapse_file = Path + "synapses_locations.txt"
voltage_file = Path + "voltage_vectors_for_unit_test.txt"
write_unit_test_vectors = WriteVectors
create_type  = 'bbpnew'
celsius = 34
print Createtxt()
os.system(Createtxt())
print("Done Test 13, Double Bouquet Cell (L4) Markram et al. Cell (2015)\n")

