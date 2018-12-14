TITLE slowly activating potassium (M-) current from bullfrog sympathetic ganglion cells

COMMENT Equations from
   Yamada WM, Koch C, Adams P (1989) In: Methods in Neuronal Modeling. MIT
   Press. 

>< Temperature adjusts time constant measured at 23.5 degC.
>< Written by Arthur Houweling for MyFirstNEURON.
ENDCOMMENT

NEURON {
	SUFFIX iM
	USEION k READ ek WRITE ik 
    RANGE gk, minf, mtau, ik
    RANGE SLOP,SHFT
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
	ek		(mV)
	gk= 0:3.1e-4	(mho/cm2)
	SLOP=0
	SHFT=0

}

STATE { m }

ASSIGNED {
	ik	(mA/cm2)
	mtau	(ms)
	htau	(ms)
	minf
	tadj
	 
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik= gk* m* (v- ek)
}

FUNCTION gate(vv){
	gate = 1/ ( 1+ exp(-(vv+ 35 + SHFT)/ (10 + SLOP)) )
	 
}

DERIVATIVE states {
       rates()
       m'= (minf- m)/ mtau 
}
  
INITIAL { UNITSOFF
	tadj= 3^ ((celsius- 23.5)/ 10)
	rates()
	m= minf
}

PROCEDURE rates() {  LOCAL a,b
	mtau= 1000/( 3.3* (exp((v+ 35+ SHFT)/ 20)+ exp(-(v+ 35+ SHFT)/ 20)))/ tadj
	minf= 1/ (1+ exp(-(v+ 35+ SHFT)/ (10+ SLOP)))
}
