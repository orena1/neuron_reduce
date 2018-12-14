: $Id: netstim.mod 2212 2008-09-08 14:32:26Z hines $
: comments at end
:###### Maybe I need to add another erand()!!!! ######
NEURON	{ 
  ARTIFICIAL_CELL NetStimO
  RANGE interval, number, start, OrEvent, stimend
  RANGE noise
  RANGE mod_freq, fire_rate, vari, factor, phase
  THREADSAFE : only true if every instance has its own distinct Random
  POINTER donotuse
}

PARAMETER {
	interval	= 10 (ms) <1e-9,1e9>: time between spikes (msec)
	number	    = 10 <0,1e9>	: number of spikes (independent of noise)
	start		= 50 (ms)	: start of first spike
	noise		= 0 <0,1>	: amount of randomness (0.0 - 1.0)
	mod_freq    = 40 <1,1e8> : the frequency the netstim will track.
	fire_rate   =  7 <0,1e8> :  the netstim firerate.
	vari        =  1  <0.00001, 100>  : the pahse preference strength
	phase       =  0  <0,2>           : the pahse releated to t
        OrEvent     = 0  : This value will contain the first event that is sent from INITAL, and will be resent in SaveState
        stimend         = 9e9 (ms)      : time of stimulus end
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
	factor
}

PROCEDURE seed(x) {
	set_seed(x)
}

PROCEDURE factors() { LOCAL temps,temps1, precision : I need to think how can I change this code so that changing the vari.
        temps  = 0
        temps1 = 0
        precision =10000
    FROM i = 0 TO precision {

        temps  = temps +  (sin( (1.5+(i/precision)*2)  * acos(-1))+1)*(1/precision)
        
        temps1 = temps1 + (( (sin((1.5+(i/precision)*2)* acos(-1))+1))^vari)*(1/precision)
    }
    :printf("temps= %g , temps1=%g\n", temps,temps1)
 factor = temps/temps1   
}

INITIAL {
	on = 0 : off
	ispike = 0
	factors()
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
		on = 1
		: randomize the first spike so on average it occurs at
		: start + noise*interval
		event = start + invl(interval) - interval*(1. - noise)
		: but not earlier than 0
		if (event < 0) {
			event = 0
		}
                OrEvent = event
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = 0
		ispike = 0
	}
}

FUNCTION invl(mean (ms)) (ms) {
    LOCAL x,u,t1,q1, lamd
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	if (noise == 0) {
		invl = mean
	}else{
	
	lamd = fire_rate*(2^vari)
	q1 = 0
	t1= (t+phase)/1000
	while (	q1==0){ :this Algorithm is from http://freakonometrics.hypotheses.org/724  from this paper. http://filebox.vt.edu/users/pasupath/papers/nonhompoisson_streams.pdf
	    u = erand()
	    t1 = t1 - log(u)/lamd
	    if (erand()<= (lamda(t1)/lamd)){
	        invl=t1*1000-(t+phase)
	        q1 =1
        }
    }
	}
}
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
		_lerand = nrn_random_pick(_p_donotuse);
	}else{
		/* only can be used in main thread */
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		: the old standby. Cannot use if reproducible parallel sim
		: independent of nhost or which host this instance is on
		: is desired, since each instance on this cpu draws from
		: the same stream
		erand = exprand(1)
VERBATIM
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
	void** pv = (void**)(&_p_donotuse);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
ENDVERBATIM
}

PROCEDURE next_invl() {
	if (number > 0) {
		event = invl(interval)
	}
	if (ispike >= number) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on spike sequence
			: but not if a netsend is on the queue
			init_sequence(t)
			: randomize the first spike so on average it occurs at
			: noise*interval (most likely interval is always 0)
			next_invl()
			event = event - interval*(1. - noise)
			net_send(event, 1)
		}else if (w < 0) { : turn off spiking definitively
			on = 0
		}
	}
	if (flag == 3) { : from INITIAL
		if (on == 1) { : but ignore if turned off by external event
			init_sequence(t)
			net_send(0, 1)
		}
	}
	if (flag == 1 && on == 1) {
		if (t <= stimend)  {
			ispike = ispike + 1
			net_event(t)
			next_invl()
		}else{
			on = 0
		}
		if (on == 1) {
			net_send(event, 1)
		}
	}
}



FUNCTION lamda(t1) { LOCAL mod_freq1
    mod_freq1 = mod_freq * 2 
    
    lamda =  factor*fire_rate*(((sin(t1*mod_freq1* acos(-1))+1))^vari)


}


COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed)

Parameters;
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   interval: 	mean time between spikes (ms)
   number: 	number of spikes (independent of noise)

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode
The intrinsic bursting parameters have been removed since
generators can stimulate other generators to create complicated bursting
patterns with independent statistics (see below)

Modified by Michael Hines to use logical event style with NET_RECEIVE
This stimulator can also be triggered by an input event.
If the stimulator is in the on==0 state (no net_send events on queue)
 and receives a positive weight
event, then the stimulator changes to the on=1 state and goes through
its entire spike sequence before changing to the on=0 state. During
that time it ignores any positive weight events. If, in an on!=0 state,
the stimulator receives a negative weight event, the stimulator will
change to the on==0 state. In the on==0 state, it will ignore any ariving
net_send events. A change to the on==1 state immediately fires the first spike of
its sequence.

ENDCOMMENT

