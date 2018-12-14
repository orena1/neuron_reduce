: SK-type calcium-activated potassium current
: Reference : 

NEURON {
       SUFFIX SK_E
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gSK_Ebar, gSK_E, ik
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gSK_Ebar = .000001 (mho/cm2)
          zTau = 1              (ms)
          ek           (mV)
          cai          (mM)
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         gSK_E	       (S/cm2)
}

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           gSK_E  = gSK_Ebar * z * z 
           ik   =  gSK_E * (v - ek)
}

DERIVATIVE states {
        rates(cai)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(ca(mM)) {
          LOCAL v
          if(ca < 1e-7){
	              ca = ca + 1e-07
          }
          zInf = 1/(1 + 0.001 / ca)
}

INITIAL {
        rates(cai)
        z = zInf
}
