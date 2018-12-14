: Dynamics that track inside and outside calcium concentrations

NEURON	{
	SUFFIX calcium
	USEION ca READ ica WRITE cai,cao
	RANGE decay, surftovol
	RANGE trans, fhspace
}

UNITS	{
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER	{
	surftovol = 10
	decay = 1500 (ms)
	cabath = 2 (mM)
	fhspace = 300 (angstrom)
	trans = 50 (ms)
}

ASSIGNED	{ica (mA/cm2)}

STATE	{
	cai (mM)
	cao (mM)
	}

BREAKPOINT	{ SOLVE states METHOD cnexp }

DERIVATIVE states	{
	cai' = -(ica*surftovol/FARADAY) - (cai-(5*1e-5))/decay
	cao' = ica*(1e6)/(fhspace*FARADAY) + (cabath-cao)/trans :unclear
}