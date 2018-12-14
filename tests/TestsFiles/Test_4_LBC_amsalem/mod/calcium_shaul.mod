: Dynamics that track inside calcium concentration

NEURON	{
	SUFFIX CaDynamics
	USEION ca READ ica WRITE cai
	RANGE decay, gamma, minCai
}

UNITS	{
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER	{
	gamma = 20 : Represents the diffusion
:		     and buffering of calcium within the cell.
:	             Should be corrected by two for Ca ion num
	decay = 2000 (ms)
	minCai = 1e-8 (mM)
}

ASSIGNED	{ica (mA/cm2)}

STATE	{
	cai (mM)
	}

BREAKPOINT	{ SOLVE states METHOD cnexp }

DERIVATIVE states	{
	cai' = -(ica*gamma/(2*FARADAY)) - (cai - minCai)/decay
}

INITIAL {
        cai = minCai
}
