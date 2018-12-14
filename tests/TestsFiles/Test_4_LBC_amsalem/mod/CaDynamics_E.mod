: Dynamics that track inside calcium concentration
:
: Etay (9.5.07)
: Gamma and Decay set based on multiple runs

NEURON	{
	SUFFIX CaDynamics_E
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
	gamma = 30 : Represents the diffusion
:		     and buffering of calcium within the cell.
:	             Should be corrected by two for Ca ion num
	decay = 800 (ms)
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
