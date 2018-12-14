TITLE SIN current

COMMENT
-----------------------------------------------------------------------------

    Sinusoidal current model for membrane impedance analysis
    ========================================================

 IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a
  point process, mimicking a current-clamp stimulation protocol, injecting
  a sinusoidally oscillating waveform I(t).
  
  I(t) = A * sin (2 pi f ( t - ttstart) )

  Note: 
  Since this is an electrode current, positive values of i depolarize the cell and in the
  presence of the extracellular mechanism there will be a change in vext since i is not a
  transmembrane current but a current injected directly to the inside of the cell.

  Note:
  This point mechanism has been created by modifying Iclamp..
          
  Refer to: Cali' et al. (2007)

 PARAMETERS

  This mechanism takes the following parameters:

  off (nA) : initial current offset.
  del (ms) : initial stimulation latency
  dur (ms) : duration of the stimulation
  amp (nA) : amplitude of the stimulation
  freq (Hz) : frequency of the sinusoidal current.

 Written by M. Giugliano and C. Cali', Brain Mind Institute, EPFL, March 2006
 Modified by C. Cali', Brain Mind Institute, EPFL, April 2006
-----------------------------------------------------------------------------
ENDCOMMENT

NEURON {
    POINT_PROCESS Isin
    RANGE del, dur, amp, i, freq, off
    ELECTRODE_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    off (nA)
    del (ms)
    dur (ms)    <0,1e9>
    amp (nA)
    freq (Hz)
}
ASSIGNED { i (nA) }

INITIAL {
    i = 0
}

BREAKPOINT {
    at_time(del)
    at_time(del+dur)

    if (t < del + dur && t >= del) {
        i = off + amp * sin(0.0062831853071795866*freq*t)
    }else{
        i = off
    }
}
