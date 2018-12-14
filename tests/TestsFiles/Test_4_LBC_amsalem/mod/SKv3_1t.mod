:Comment :
:Reference : :		Characterization of a Shaw-related potassium channel family in rat brain, The EMBO Journal, vol.11, no.7,2473-2486 (1992)
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX SKv3_1t
	USEION k READ ek WRITE ik
	RANGE gSKv3_1tbar, gSKv3_1t, ik 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gSKv3_1tbar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gSKv3_1t	(S/cm2)
	mInf
	mTau
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gSKv3_1t = gSKv3_1tbar*m
	ik = gSKv3_1t*(v-ek)
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
  LOCAL qt
  qt = 2.3^((34-21)/10)
	UNITSOFF
		mInf =  1/(1+exp(((v -(18.700))/(-9.700))))
		mTau =  0.2*20.000/(1+exp(((v -(-46.560))/(-44.140))))/qt
	UNITSON
}
