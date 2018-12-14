: SK-type calcium-activated potassium current


NEURON {
       SUFFIX SK_S
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gSK_Sbar, gSK_S, ik
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gSK_Sbar = 0 (mho/cm2)
          zTau = 1              (ms)
          ek           (mV)
          cai          (mM)
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         gSK_S	       (S/cm2)
}

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           gSK_S  = gSK_Sbar * z * z 
           ik   =  gSK_S * (v - ek)
}

DERIVATIVE states {
        rates(v,cai)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(Vm (mV), ca(mM)) {
          LOCAL v
          v = Vm + 5
          if(ca < 1e-7){
	              ca = ca + 1e-07
          }
          zInf = 1/(1 + 0.001 / ca)
}

INITIAL {
        rates(v,cai)
        z = zInf
}

