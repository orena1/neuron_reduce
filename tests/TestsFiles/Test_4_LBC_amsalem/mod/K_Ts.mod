:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential

NEURON	{
	SUFFIX K_Ts
	USEION k READ ek WRITE ik
	RANGE gK_Tsbar, gK_Ts, ik
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gK_Tsbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Ts	(S/cm2)
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
	gK_Ts = gK_Tsbar*(m^4)*h
	ik = gK_Ts*(v-ek)
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
		mInf =  1/(1 + exp(-(v+0)/19))
		mTau =  (0.34+0.92*exp(-((v+71)/59)^2))
		hInf =  1/(1 + exp(-(v+66)/-10))
		hTau =  (8+49*exp(-((v+73)/23)^2))
		v = v - 10
	UNITSON
}
