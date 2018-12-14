NEURON {
POINT_PROCESS GapMul
RANGE g, i, vgap
NONSPECIFIC_CURRENT i
THREADSAFE
}
PARAMETER { g = 1 (nanosiemens) }
ASSIGNED {
v (millivolt)
vgap (millivolt)
i (nanoamp)
}
BREAKPOINT { i = (v - vgap)*(g*1e-3) }
