:Comment :
:Reference :Magistretti.J and Alonso.A, J. Gen. Physiol, 1999

NEURON	{
	SUFFIX Nap
	USEION na READ ena WRITE ina
	RANGE gNapbar, gNap, ina 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNapbar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNap	(S/cm2)
	mInf
	mTau
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNap = gNapbar*m*m*m
	ina = gNap*(v-ena)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
	UNITSOFF
		mInf = 1.0000/(1+ exp((v - -44.0000)/-4.8500))
		mTau = 1
	UNITSON
}
