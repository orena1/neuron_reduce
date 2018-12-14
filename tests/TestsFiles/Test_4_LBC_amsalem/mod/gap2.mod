NEURON {
POINT_PROCESS Gap2
RANGE vgap
RANGE r, i
NONSPECIFIC_CURRENT i
}
PARAMETER { r = 1e10 (megohm) }
ASSIGNED {
v (millivolt)
vgap (millivolt)
i (nanoamp)
}
BREAKPOINT { i = (v - vgap)/r }
