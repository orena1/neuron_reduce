TITLE Four-State Model of chanrhod channel for Subthalamic Nucleus

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (mW) = (milliwatt)
    (photons)  = (1)
}

NEURON { :public interface of the mechanism
	SUFFIX chanrhod
	:Name of the Channel
    NONSPECIFIC_CURRENT	icat
    :Specific Channel for Na, Ca, K
    RANGE photons, flux, irradiance, U, U0, U1
    :Calculated optics values given source_values
    RANGE channel_density, gdens1, gdens2
    :Channel density and conductance. gdens1 = density conductance for o1. gdens2 = density conductance for o2
    RANGE x, y, z
    :Location of that segment (of nsegs)
	:RANGE Ka1, Ka2, Kd1, Kd2, e12, e21, e12dark, e21dark, Kr, delta1, delta2, o10, o20, Islow, Ifast, c_1, c_2, b, h, c, amp, photon_energy, wavelength, phi0, phio, gcat
	:Various Rate Constants, phi0:constant, phio: static value of phi
    GLOBAL source_photons, source_flux, source_irradiance
    :Optical Input
    RANGE Ka1, Ka2
    GLOBAL Kd1, Kd2
    GLOBAL Kr
    RANGE e12, e21, e12dark, e21dark
    RANGE delta1, delta2
    RANGE o10, o20
    RANGE Islow, Ifast
    RANGE c_1, c_2, b
    GLOBAL h, c
    RANGE amp, photon_energy
    RANGE wavelength
    RANGE phi0, phio, gcat    
    GLOBAL sigma_retinal, gcat1, gcat2, ecat, Imax, gamma, tChR
    :sigma_retinal = cross-sectional area of chanrhod, gcat1/2 = single channel conductance of chanrhod, ecat = nernst potential for chanrhod, Imax=maximum current if all channels are in o1, gamma=g2/g1
    GLOBAL epsilon1, epsilon2
    :epsilon1,2 = quantum efficiency for ka1, ka2
    RANGE Tx, phi
    :Tx = Transfer Resistance between two point locations assuming homogenous tissue (V=ITx), phi = photon flux / channel, phi = flux x Tx, flux = irradiance x cross section of retinal
    GLOBAL tstimon, tstimoff
}

:NEURON doesn't change ion concentrations automatically - need another mechanism that will write cai and cao

PARAMETER {
    :channel_density        = 0
    channel_density        = 1.3e10              (1/cm2) : variable : number of channels per cm2
    : 2e9 per oocyte : Nagel 2003
    : 100 ChR2s/um2 = 1e10 (Grossman 2011)

    :gcat1    = 40e-15    (mho)   40 fS
    :gcat2    = 40e-15    (mho)   40 fS
    gcat1    = 50e-15    (mho)   : (Grossman 2011) 50 fS
    gcat2    = 250e-17   (mho)   : Figure 8 Nickolic 2009
    :g=150e-15 (mho) 15fA at -100mV (Harz 1992, Feldbauer 2009, Lin 2009)
    : 100 fS : Channelrhodopsins: Molecular Properties and Applications Ernst Bamberg, PhD
    : single-channel conductance of ~=50 fS : Channelrhodopsin-2, a directly light-gated cation-selective membrane channel Georg Nagel
    : A value of 40 fS was obtained : Channelrhodopsin-2 is a leaky proton pump Katrin Feldbauer

    gamma=0.05 : Figure 8 Nickolic 2009

    sigma_retinal = 1.2e-16  (cm2)        : 1.2e-8 um2 which is cross-sectional area of retinal

    epsilon1 = 0.5             (1)         : quantum efficiency of ChR-2 system : A typical value would be e 0.5 for rhodopsin; nikolic 2009, Hegemann (1999), Constant Value (Nikolic 2009)
	:epsilon2 = 0.12                        : (Grossman 2011)
	epsilon2 = 0.1                      : Figure 8 at source_irradiance = 0.9 mW/mm2 (Nikolic 2009)
    ecat     = 0      (mV)     : Nagel 2003
                               : -42 mV: Kang, Y., Okada, T., Ohmori, H., 1998. A phenytoin-sensitive cationic current... 1998.
                               : -48 mV: Christian Alzheimer, A novelvoltage-dependentcationcurrentinratneocortical neurones. 1994.

    Tx      = 1       (1)      : Default light "transfer resistance" between optrode and compartment; geometry dependent
    vshift  = 0        (mV)     : Adjust the voltage for different resting potentials (resting potential of pyramidal cell is -65 and of fiber tract is -70)

	:tChR=0.2 (ms) : Volvox - Figure 8 (Nickolic 2009)
	tChR=1.3 (ms) : Axons (Nickolic 2009)

    x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)
	Ka1 = 0.5 :(Grossman 2011)
	Ka2 = 0.12 :(Grossman 2011)

	:Kd1 = 0.1 : Constant Value (Grossman 2011)
	:Kd2 = 0.05 : Constant Value (Grossman 2011)
	Kd1 = 0.13 : Figure 8 at source_irradiance = 0.9 mW/mm2 (Nikolic 2009) - 0.09W/cm2
	Kd2 = 0.025 : Figure 8 at source_irradiance = 0.9 mW/mm2 (Nikolic 2009)

	Kr = 0.0004 : Figure 8 (Nikolic 2009)

	e12 = 0.053 : Figure 8 at source_irradiance = 0.9 mW/mm2 (Nikolic 2009)
	e21 = 0.023 : Figure 8 at source_irradiance = 0.9 mW/mm2 (Nikolic 2009)
	e12dark = 0.022
	e21dark = 0.011

	delta1= 0.03 (ms): Figure 8 (Nikolic 2009)
	delta2= 0.15 (ms): Figure 8 (Nikolic 2009)

	h = 6.6260693e-34           (m2 kg/s)  : planck's constant
    c = 299792458.0             (m/s)      : speed of light
	wavelength = 4.45e-7
}

ASSIGNED {  :calculated by the mechanism (computed by NEURON)
    v           (mV)
    icat        (mA/cm2)
    gdens1        (mho/cm2)
    gdens2        (mho/cm2)
    source_irradiance  (W/cm2)      : Light irradiance (W/mm2) exiting optrode, from ostim.mod
    source_photons     (photons/ms) : number of photons exiting optrode per millisecond, from ostim.mod
    source_flux        (photons/ms cm2) : flux of photons exiting optrode per millisecond, from ostim.mod
    irradiance          (W/cm2)      : number of photons exiting optrode per millisecond, from ostim.mod
    flux               (photons/ms cm2) : number of photons exiting optrode per millisecond, from ostim.mod
    phi                (photons/ms) : number of photons hitting channel per millisecond
    U
    U0
    U1
    Imax
    Islow
    Ifast
    c_1
    c_2
    b
    amp
    photon_energy
    phi0
    phio
    gcat
    tstimon
    tstimoff
}

STATE { :state or independent variables
	o1 o2 c1 c2
}

INITIAL {
    irradiance = 0
    flux = 0
    phi = 0
    Islow=0
    Ifast=0
    :gdens1 = gcat1 * channel_density : (mho/cm2)
	:gdens2 = gcat2 * channel_density : (mho/cm2)
    tstimon = 0
    tstimoff = 0

    : STATES
	c1 = 1 :Amount of channels at initial time
	c2 = 0
	o1 = 0
	o2 = 0

	phio = 0
    o10=0
    o20=0
}

BREAKPOINT {
    irradiance = source_irradiance * Tx : (W/cm2)
    flux      = source_flux * Tx           : (photons/ms cm2)
    phi       = flux             * sigma_retinal
                : (photons/ms cm2) * (cm2)
                :  --> (photons/ms / channel)

	U=v-ecat-vshift-75
	U0=40
	U1=15
	Imax=(v-ecat-vshift)*gcat1*channel_density : mA/cm2

    b= (Kd1+Kd2+e12dark+e21dark)/2

    c_1 = 0.1029797709
    c_2 = 0.0398631371

    gcat=(o1+gamma*o2)*gcat1*channel_density

	if (phi>0) {
		Ka1 = epsilon1 * phi * (1 - exp( -(t - tstimon) / tChR)) :Need to add t-td for optical stimulation starting after t=0
		Ka2 = epsilon2 * phi * (1 - exp( -(t - tstimon) / tChR))
		:e12 = 0.011 + 0.005*log(phi/0.024) :(Grossman 2011) - have to use at suprathreshold phi values - need to have equation for phi at certain properities?
		:e21 = 0.008 + 0.004*log(phi/0.024) :(Grossman 2011)
		e12=0.053
		e21=0.023
		:e12 = e12dark + c_1*log10(1+phi/phi0) : e12=0.053 phi>0, e12=0.022 phi=0
   		:e21 = e21dark + c_2*log10(1+phi/phi0) : e21=0.023 phi>0, e21=0.011 phi=0
		icat = Imax * (o1 + gamma * o2) * (1-exp(-U/U0))/(U/U1) : (mA/cm2)
		o10=o1
		o20=o2
		phio=phi
	} else {
		Ka1 = epsilon1 * phio * (exp (- ((t-tstimoff) / tChR)) - exp( -(t - tstimon) / tChR) ) : 500ms duration light pulse - very rapidly drops off due to long duration pulse and short tChR
		Ka2 = epsilon2 * phio * (exp (- ((t-tstimoff) / tChR)) - exp( -(t - tstimon) / tChR) ) : 500ms duration light pulse
		e12 = 0.022
		e21 = 0.011
		:icat = Imax * (o1 + gamma * o2) * (1-exp(-U/U0))/(U/U1) : (mA/cm2)
    	Islow = Imax * ( ( (delta2 - (Kd1 + (1 - gamma) * e12dark)) * o10 + (( 1 - gamma) * e21dark + gamma*(delta2-Kd2)) * o20 ) / (delta2 - delta1) )
    	Ifast = Imax * ( ( (Kd1 + (1 - gamma) * e12dark - delta1) * o10 + (-(1 - gamma)*e21dark+gamma*(Kd2-delta1))*o20 ) / (delta2 - delta1) )
    	icat = Islow * exp(-delta1*(t-tstimoff)) + Ifast*exp(-delta2*(t-tstimoff)) : (mA/cm2)
	}

    SOLVE states METHOD cnexp
    if (o1>1){o1=1}
    if (o1<0){o1=0}
	if (o2>1){o2=1}
    if (o2<0){o2=0}
    if (c1>1){c1=1}
    if (c1<0){c1=0}
    if (c2>1){c2=1}
    if (c2<0){c2=0}
    c1    = 1 - o1 - o2 - c2
}

DERIVATIVE states {  :states the set of diffy qs
	o1' = Ka1*c1 - (Kd1 + e12)*o1 + e21*o2
	o2' = Ka2*c2 + e12*o1 - (Kd2+e21)*o2
	:c1' = c1* (Kd1 - Ka1) + c2 * Kr
	c2' = Kd2*o2 - (Ka2+Kr)*c2
}
