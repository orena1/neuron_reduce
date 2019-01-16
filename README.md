Introduction
===========

Neuron_Reduce provides an analytical method for reducing neuron model complexity. It enables the mapping of synapses and active ion channels to a computationally simpler model while accelerating simulation speed by up to 200-fold for inputs consisting of thousands of dendritic synapses. Full details are available in the accompanied paper - https://doi.org/10.1101/506485

Installation
===========

``` pip install --user neuron_reduce```

Quick Start
===========
The following code show the main function that is used to reduce a complex cell. 
```python
complex_cell  # The model cell
synapses_list # A list of all synapse on this cell
netcon_list   # A list of all netcons for the synapses on the cell

import neuron_reduce
reduced_cell, synapses_list, netcons_list =  neuron_reduce.subtree_reductor(complex_cell, synapses_list, netcons_list)
```

Detailed example
===========

Copy example folder from github
```bash
git clone https://github.com/orena1/neuron_reduce.git
```

Go to example folder
```bash
cd neuron_reduce
cd example
nrnivmodl mod #compile the mod files
```

Open python and run the following code


```python
from __future__ import division
from neuron import gui,h
import numpy as np
import neuron_reduce
import time
import matplotlib.pyplot as plt



#Create a L5_PC model
h.load_file('L5PCbiophys3.hoc')
h.load_file("import3d.hoc")
h.load_file('L5PCtemplate.hoc')
complex_cell = h.L5PCtemplate('cell1.asc')
h.celsius = 37
h.v_init = complex_cell.soma[0].e_pas


#Add synapses to the model
synapses_list, netstims_list, netcons_list, randoms_list = [], [], [] ,[]

all_segments = [i for j in map(list,list(complex_cell.apical)) for i in j] + [i for j in map(list,list(complex_cell.basal)) for i in j]
len_per_segment = np.array([seg.sec.L/seg.sec.nseg for seg in all_segments])
rnd = np.random.RandomState(10)
for i in range(10000):
    seg_for_synapse = rnd.choice(all_segments,   p=len_per_segment/sum(len_per_segment))
    synapses_list.append(h.Exp2Syn(seg_for_synapse))
    if rnd.uniform()<0.85:
        e_syn, tau1, tau2, spike_interval, syn_weight = 0, 0.3, 1.8,  1000/2.5, 0.0016
    else:
        e_syn, tau1, tau2, spike_interval, syn_weight = -86, 1,   8,   1000/15.0, 0.0008
    #set synaptic varibales
    synapses_list[i].e, synapses_list[i].tau1, synapses_list[i].tau2 = e_syn, tau1, tau2
    #set netstim variables
    netstims_list.append(h.NetStim())
    netstims_list[i].interval, netstims_list[i].number, netstims_list[i].start, netstims_list[i].noise = spike_interval, 9e9, 100, 1
    #set random
    randoms_list.append(h.Random())
    randoms_list[i].Random123(i)
    randoms_list[i].negexp(1)
    netstims_list[i].noiseFromRandom(randoms_list[i])       
    #set netcon varibales 
    netcons_list.append(h.NetCon(netstims_list[i], synapses_list[i] ))
    netcons_list[i].delay, netcons_list[i].weight[0] = 0, syn_weight

#Simulate the full neuron for 1 seconds
soma_v = h.Vector()
soma_v.record(complex_cell.soma[0](0.5)._ref_v)

time_v = h.Vector()
time_v.record(h._ref_t)

h.tstop = 1000
st = time.time()
h.run()
print('complex cell simulation time {:.4f}'.format(time.time()-st))
complex_cell_v = list(soma_v)



#apply Neuron_Reduce to simplify the cell
reduced_cell, synapses_list, netcons_list = neuron_reduce.subtree_reductor(complex_cell, synapses_list, netcons_list, reduction_frequency=0, total_segments_manual=-1)
for r in randoms_list:r.seq(1) #reset random


#Running the simulation again but now on the reduced cell
st = time.time()
h.run()
print('reduced cell simulation time {:.4f}'.format(time.time()-st))
reduced_celll_v = list(soma_v)

#plotting the results
plt.figure()

plt.plot(time_v, complex_cell_v, label='complex cell')
plt.plot(time_v, reduced_celll_v,  label='redcued cell')
plt.show()
```

Citation
===========
https://doi.org/10.1101/506485
