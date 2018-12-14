: test units: modlunit ostim.mod

TITLE Optical Stimulation

NEURON {
    POINT_PROCESS ostim
    RANGE radius
    RANGE delay, amp, dur
    RANGE pulses, isi, pcount
    POINTER irradiance
    POINTER flux
    POINTER photons
    POINTER tstimon
    POINTER tstimoff
    RANGE wavelength
}

UNITS {
    (mW) = (milliwatt)
    (W)  = (watt)
    PI   = (pi)             (1)

}

PARAMETER {
    radius    = 100             (um)
    delay     = 1               (ms)
    dur       = 5               (ms)
    : Aravanis et al 2007
    :     Psrc = 20 mW
    :     Af = 0.6
    :     r = 0.1 mm
    :     Irradiance = Af * Psrc/(pi * r**2) = 380 mW/mm2 --> 0.38 W/mm2 --> 38 W/cm2
    amp = 38                    (W/cm2)
    wavelength = 4.73e-7        (m)
    h = 6.6260693e-34           (m2 kg/s)  : planck's constant
    c = 299792458.0             (m/s)      : speed of light

    : Number of pulses
    pulses = 1
    isi    = 1                 (ms)        : Inter-spike interval
}

ASSIGNED {
    on
    irradiance (W/cm2)         : W/cm2
    flux      (1/ms cm2)       : photons/ms cm2
    photons   (photons/ms)     : photons/ms
    pcount
    tstimon
    tstimoff
}

INITIAL {
    on          = 0
    pcount      = pulses
    irradiance   = 0
    flux        = 0
    photons     = 0
    net_send(delay,1)          : Sends a packet that should arrive after given initial_delay
    tstimon     = 0
    tstimoff    = 0
}

BEFORE BREAKPOINT {
    if (on==0) {
        irradiance = 0
        photons   = 0
        flux      = 0
    } else if (on==1) {
        calc_irradiance_photons_flux()
    }
}

NET_RECEIVE (w) {
    if (flag==1) {
        on = 1
        tstimon = t
        pcount = pcount - 1
        if (pcount>0) {
            net_send(dur,3)
        } else {
            net_send(dur,2)
        }
    }else if (flag==2) { : stimulus complete
        on = 0
        tstimoff = t
    }else if (flag==3) { : spike train
        on = 0
        net_send(isi,1)
        tstimoff = t
    }

}

FUNCTION calc_irradiance_photons_flux() {
    LOCAL photon_energy, area
    irradiance = amp
    photon_energy = h * c / wavelength : J / photon
    flux = (1 / 1000) * amp / photon_energy : (1 s / 1000 ms) * (W/cm2) * (photons/W s)
                                              :     --> photons/ms cm2
    area = PI * (radius ^ 2) : um2
    photons = 1e-8 * flux * area : (1 cm2/1e8 um2) * (um2) * (photons/ms cm2)
}

