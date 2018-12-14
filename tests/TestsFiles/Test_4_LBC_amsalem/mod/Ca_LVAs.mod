:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential

NEURON	{
	SUFFIX Ca_LVAs
	USEION ca READ eca WRITE ica
	RANGE gCa_LVAsbar, gCa_LVAs, ica
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_LVAsbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_LVAs	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_LVAs = gCa_LVAsbar*m*m*h
	ica = gCa_LVAs*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	UNITSOFF
		v = v + 10
		mInf = 1.0000/(1+ exp((v - -30.000)/-6))
		mTau = 5.0000 + 20.0000/(1+exp((v - -25.000)/5))
		hInf = 1.0000/(1+ exp((v - -80.000)/6.4))
		hTau = 20.0000 + 50.0000/(1+exp((v - -40.000)/7))
		v = v - 10
	UNITSON
}
