: Calcium channel taken from Dayan & Abbot and modified with shifts

NEURON	{
	SUFFIX cat
	USEION ca READ eca WRITE ica
	RANGE gcabar, gca, ical, catempfactor
    RANGE MSHFT, MSLOP, HSHFT, HSLOP
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gcabar = 0.0001 (S/cm2)
	MSHFT = -40
	MSLOP = 0
	HSHFT = -40
	HSLOP = 0
	catempfactor = 1
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gca	(S/cm2)
}

STATE	{
	mcal
	hcal
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gca = gcabar*mcal*mcal*hcal
	ica = gca*(v-eca)
}

INITIAL	{
	mcal = mcalinf(v)
	hcal = hcalinf(v)
}
FUNCTION gate(vv){LOCAL lm,lh
	lm = mcalinf(vv)
	lh = hcalinf(vv)
	gate = lm*lm*lh
}

DERIVATIVE states	{
    LOCAL inf, tau
    inf = mcalinf(v)
    tau = mcaltau(v)
    mcal' = ((inf-mcal)/tau)
    inf = hcalinf(v)
    tau = hcaltau(v)
    hcal' = ((inf-hcal)/tau)
}

FUNCTION mcalinf(Vm (mV))	{ :what units?
	UNITSOFF
	:mcalinf = 0.204 + 0.333 / ( exp( ( Vm + 15.8 + MSHFT) / 18.2 ) + exp( ( - Vm - 131 - MSHFT) / 16.7 ) )
	mcalinf = 1/(1+exp(- ((Vm+57+MSHFT)/(6.2+MSLOP))) )
	UNITSON
}

FUNCTION hcalinf(Vm (mV))	{ :what units?
	UNITSOFF
	hcalinf = 1/(1+exp(((Vm+81+HSHFT)/(4+HSLOP))))
	UNITSON
}

FUNCTION mcaltau(Vm	(mV))	{
	UNITSOFF
	mcaltau = (0.612+(1/(exp(-(Vm+132)/16.7)+exp((Vm+16.8)/18.2))))*catempfactor
	UNITSON
}

FUNCTION hcaltau(Vm	(mV))	{
	UNITSOFF
:	if (Vm < -80)	{
:		hcaltau = (exp((Vm+467)/66.6))*catempfactor
:	}else{
		hcaltau = (28+exp(-(Vm+22)/10.5))*catempfactor
:	}
	UNITSON
}

