TITLE SINUOISY current

COMMENT
--------------------------------------------------------------------------------------------------------------------

    Sinusoidal + Fluctuating current model for temporally-modulated synaptic bombardment
    ====================================================================================

  The present version implements and generate a realization of an Ornstein-Uhlenbeck (OU) process
  (see Cox & Miller, 1969; see Tuckwell) to mimick the somatic impact of linearly adding EPSPs and
  IPSPs. Thus, it generates and injects in the specified neuronal compartment a fluctuating current
  waveform (a noise), characterized by a gauss-distributed amplitude, where neighboring amplitude 
  samples are by definition linearly correlated on a time scale set by the correlation time-length
  "tau" of the process.
 
  The numerical scheme for integration of OU processes takes advantage of the fact that these
  processes are gaussian, which led to an exact update rule independent of the time step dt..
  (see Gillespie DT, Am J Phys 64: 225, 1996):

  x(t+dt) = x(t) + (1. - exp(-dt/tau)) * (m - x) + sqrt(1.-exp(-2.*dt/tau)) * s * N(0,1)  
  where N(0,1) is a normal random number (avg=0, sigma=1)..

  Please note that only fixed integration time-step methods makes sense, since the stochastic current
  synthesized by the present mechanism is produced randomly and on-line. In other words, it is wrong to
  assume that neglecting the present integration step, reducing it and resynthesizing the current, lead
  to the same overall trajectory in the compartment output voltage.

  As opposed to the previously developed mechanisms Ifluct1.mod (see the ModelDB), here the STANDARD DEVIATION of the
  current is sinusoidally oscillating as indicated below, as a function of time: A * sin (2 pi f t ).

 IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a point process, mimicking a current-
  clamp stimulation protocol, injecting a sinusoidally oscillating waveform overlapped to a noisy component.
  
  I(t)  = m + x(t) * (A * sin (2 pi f t ) + s)

  Note: 
  Since this is an electrode current, positive values of i depolarize the cell and in the presence of the
  extracellular mechanism there will be a change in vext since i is not a transmembrane current but a current
  injected directly to the inside of the cell.
  
 REFERENCES

 Koendgen, H., Geisler, C., Fusi, S., Wang, X.-J., Luescher, H.-R., and Giugliano, M. (2007). 
 Arsiero, M., Luescher, H.-R., Lundstrom, B.N., and Giugliano, M. (2007). 
 La Camera, G., Rauch, A., Thurbon, D., Luescher, H.-R., Senn, W., and Fusi, S. (2006). 
 Giugliano, M., Darbon, P., Arsiero, M., Luescher, H.-R., and Streit, J. (2004). 
 Rauch, A., La Camera, G., Luescher, H.-R., Senn, W., and Fusi, S. (2003). 
 Boucsein, C. Tetzlaff, T., Meier, R., Aertsen, A., and B. Naundorf (2009)

 The present mechanism has been inspired by "Gfluct.mod", by A. Destexhe (1999), as taken from ModelDB.
 Destexhe, A., Rudolph, M., Fellous, J-M. and Sejnowski, T.J. (2001). 

 AUTHORS
 M. Giugliano, Theoretical Neurobiology, Department of Biomedical Sciences, University of Antwerp, Antwerp
 		and Brain Mind Institute, EPFL Lausanne
  
 V. Delattre, Brain Mind Institute, EPFL Lausanne


 PARAMETERS

  The mechanism takes as input the following parameters (reported with their default values):

    m   = 0. (nA)       : DC offset of the overall current
    s   = 0.5 (nA)       : square root of the steady-state variance of the (noisy) stimulus component
    tau = 2. (ms)       : steady-state correlation time-length of the (noisy) stimulus component
    amp = 0. (nA)       : amplitude of the (sinusoidal) stimulus component
    freq= 0. (Hz)       : steady-state correlation time-length of the (noisy) stimulus component

--------------------------------------------------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS Isinunoisy2
    RANGE amp, i, freq, m, s, tau, x, new_seed
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    dt   (ms)
    m   = 0. (nA)       : DC offset of the overall current
    s   = 0.5 (nA)       : square root of the steady-state variance of the (noisy) stimulus component
    tau = 0. (ms)       : steady-state correlation time-length of the (noisy) stimulus component
    amp = 0. (nA)       : amplitude of the (sinusoidal) stimulus component
    freq= 0. (Hz)       : frequency of the sinusoidal input current 
    phas= 0. (HZ)       : phase of the sinusoidal input current
    fr2 = 0.(Hz)        : steady-state correlation time-length of the (noisy) stimulus component
}

ASSIGNED { 
    i (nA)              : overall sinusoidal noisy current
    x                   : state variable
}

INITIAL {
    i = m
    x = 0               : to reduce the transient, the state is set to its (expected) steady-state    
}

BREAKPOINT {  
    SOLVE oup
    
    if (tau <= 0) {  x =  normrand(0,1)  }
     
    :if (amp>s) { printf("ERROR..! amp: %f s: %f \n",amp,s) }
     i = x * (s + amp * sin(0.0062831853071795866 * freq * t)) + m

}


PROCEDURE oup() {       : uses "Scop" function normrand(mean, std_dev)
if (tau > 0) {  x = x + (1. - exp(-dt/tau)) * ( - x) + sqrt(1.-exp(-2.*dt/tau)) * normrand(0,1)}
}

PROCEDURE new_seed(seed) {      : procedure to set the seed
    set_seed(seed)
    VERBATIM
      printf("Setting random generator with seed = %g\n", _lseed);
    ENDVERBATIM
}

