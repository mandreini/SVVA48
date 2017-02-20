# blah blah heading
# this file determines the forces that will be acting on the beam

import matplotlib.pyplot as plt
import numpy
import parameters

step_size = 0.001
z_disc = numpy.arange(0, parameters.L, step_size)

def getValAtLoc(force, z_loc):
    # used to check the value of the force at a given location
    # determine which index of z_disc is closest to given z_loc
    z_ind = numpy.where(min(abs(z_disc-z_loc)) == abs(z_disc-z_loc))
    return force[z_ind]


def normalInternal():
    return 0


def shearyInternal():
    # determines the internal shear force along the fuselage in the y-direction
    sheary = - parameters.q * z_disc
    sheary = sheary + numpy.where(z_disc>parameters.Lf02, parameters.Sy2, 0.)
    sheary = sheary + numpy.where(z_disc>parameters.Lf01, parameters.Sy1, 0.)

    assert abs(sheary[0]) < 0.001, 'Internal shear_y at aft of fuselage is not 0.'
    assert abs(sheary[-1] - parameters.q * step_size) < 0.001, 'Internal shear_y at nose of fuselage is not 0.'

    return sheary
    # L-Lf1-Lf2 is about 33.8 m, which is index 33801. shear[33800] is similar to shear[33801] but shear[33802] is much different


def shearyInternal_alt(forces, locations):
    # the high-brow thing i wanted to try to make
    sheary_val = numpy.array([0.] * len(z_disc))
    for ind, f in enumerate(forces):
        if isinstance(locations[ind], tuple) or isinstance(locations[ind], list):
            # distributed load
            dist_start, dist_end = locations[ind]
            sheary_val += numpy.where((z_disc>=dist_start) * (z_disc<=dist_end), f*z_disc, 0.)
        elif isinstance(locations[ind], float):
            # point force
            sheary_val += numpy.where(z_disc>=locations[ind], f, 0.)

    assert abs(sheary_val[0]) < 0.001, 'Internal shear_y at aft of fuselage is not 0.'
    assert abs(sheary_val[-1] - parameters.q * step_size) < 0.001, 'Internal shear_y at nose of fuselage is not 0.'

    return sheary_val


def shearxInternal():
    # determines the internal shear force along the fuselage in the x-direction
    shearx = parameters.Sx
    shearx = shearx + numpy.where(z_disc>parameters.Lf02, parameters.Sx2, 0)
    shearx = shearx + numpy.where(z_disc>parameters.Lf01, parameters.Sx1, 0)

    # not sure how to check that the shear calculations work out, as the shear starts at non-zero.
    # also, I don't know how Sx1, Sx2 are to be determined. Will need to talk to Phillip
    return shearx


def momentxInternal():
    # determines the internal moment in x along the fuselage
    momentx = - parameters.q * z_disc * z_disc/2
    momentx = momentx + numpy.where(z_disc>parameters.Lf02, parameters.Sy2 * (z_disc - parameters.Lf02), 0.)
    momentx = momentx + numpy.where(z_disc>parameters.Lf01, parameters.Sy1 * (z_disc - parameters.Lf01), 0.)

    assert abs(momentx[0] < 0.001), 'Internal moment_x at aft of fuselage is not 0.'
    assert abs(momentx[-1] + parameters.q * step_size * step_size / 2) < 0.001, 'Internal moment_x at nose of fuselage is not 0.'

    return momentx


def momentxInternal_alt(forces, locations):
    # more pretentious stuff
    momentx_val = numpy.array([0.] * len(z_disc))
    for ind, f in enumerate(forces):
        if isinstance(locations[ind], tuple) or isinstance(locations[ind], list):
            dist_start, dist_end = locations[ind]
            momentx_val += numpy.where((z_disc>=dist_start) * (z_disc<=dist_end), f*z_disc*z_disc/2, 0.)
        elif isinstance(locations[ind], float):
            f_loc = locations[ind]
            momentx_val += numpy.where(z_disc>f_loc, f*(z_disc-f_loc), 0.)

    assert abs(momentx_val[0] < 0.001), 'Internal moment_x at aft of fuselage is not 0.'
    assert abs(momentx_val[-1] + parameters.q * step_size * step_size / 2) < 0.001, 'Internal moment_x at nose of fuselage is not 0.'

    return momentx_val


def momentzInternal():
    # determines the internal moment in z along the fuselage
    momentz = parameters.M3 + parameters.Sx*parameters.R  # this is taken from equation 2.4; momentz = Mz,A-Mz,B
    return momentz


sy = shearyInternal_alt(parameters.external_forces_y, parameters.external_locations_y)
pass
