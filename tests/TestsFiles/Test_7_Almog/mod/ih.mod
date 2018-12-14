TITLE hyperpolarization-activated current (H-current) 

COMMENT
Based on Williams and Stuart J. Neurophysiol 83:3177,2000
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iH
	USEION h READ eh WRITE ih VALENCE 1
	RANGE gbar, h_inf, tau, ih  
	GLOBAL t0,t1,t2,t3,v05, z
	GLOBAL q10, temp, tadj, vmin,vmax
	
}

UNITS {
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(molar)	= (1/liter)
	(mM) 	= (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)

}

PARAMETER {
	v		(mV)
	
	celsius		(degC)
	eh	   (mV)     
	gbar	= 0.0	(pS/um2)

	v05 = -91       (mV)   		: V1/2 of activation	
	z=6		(mV)	 	: slope of activation
	
	t0 = 0.0003933	(1/ms) 	 	: parameters for time constant of activation    
	t1 = -0.0249	(1/mV)     
	t2 = 0.0877	(1/ms)     
	t3 = 0.062	(1/mV)
			
	temp = 21	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity
	vmin = -120 (mV)
	vmax = 100 (mV)     
	
}


ASSIGNED {
	ih		(mA/cm2)
        h_inf
        tau        (ms)
	tadj
	
}

STATE { h }


INITIAL {
	rates(v)
      	h = h_inf
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
      	ih = (1e-4) * gbar * h * (v-eh)
}


DERIVATIVE states  { 

	rates(v) 
	h' = (h_inf-h)/tau  
}


PROCEDURE rates( v (mV)) {

	tadj= q10^((celsius-22)/10)
	h_inf = 1/(1+exp((v-v05)/z))	
        tau = 1/(tadj*(t0*exp(t1*v)+t2*exp(t3*v)))			
}	

