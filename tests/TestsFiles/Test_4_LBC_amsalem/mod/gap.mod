NEURON {
POINT_PROCESS Gap
POINTER vgap
RANGE g, i
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
