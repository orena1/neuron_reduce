NEURON	{
	SUFFIX caint
	USEION ca READ ica WRITE cai
	RANGE decay, surftovol
	RANGE init_cai
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
	init_cai = 0
}

ASSIGNED	{ica (mA/cm2)}

STATE	{cai (mM) }

INITIAL{
	cai = init_cai
}

BREAKPOINT	{ SOLVE state METHOD cnexp }

DERIVATIVE state	{
	cai' = -(ica*surftovol/FARADAY) - cai/decay
}
