# blah blah heading
# this file determines the forces that will be acting on the beam

import matplotlib.pyplot as plt
import numpy
import parameters

step_size = 1.
z_disc = numpy.arange(0, parameters.L, step_size)

def getValAtLoc(force, z_loc):
    """
    Determine the internal (given) at a location (given)
    :param force: numpy.ndarray(force)
    :param z_loc: float(z)
    :return: float(force)
    """
    z_ind = numpy.where(min(abs(z_disc-z_loc)) == abs(z_disc-z_loc))
    return force[z_ind]


def shearyInternal():
    """
    This method determines the shear in y along the beam in z. Takes into account the distributed load (q) and both
    reaction forces (Sy1, Sy2) due to the landing gear
    :return: numpy.ndarray(sheary)
    """
    sheary = - parameters.q * z_disc
    sheary = sheary + numpy.where(z_disc>parameters.Lf02, parameters.Sy2, 0.)
    sheary = sheary + numpy.where(z_disc>parameters.Lf01, parameters.Sy1, 0.)

    assert abs(sheary[0]) < 0.001, 'Internal shear_y at aft of fuselage is not 0.'
    assert abs(sheary[-1] - parameters.q * step_size) < 0.001, 'Internal shear_y at nose of fuselage is not 0.'

    return sheary


def shearyInternal_alt(forces, locations):
    """
    This is a more general way to determine the shear in y along the fuselage in z. forces is a list of the applied
    forces, point or distributed. locations is a list of the locations of said forces; distributed forces are to be
    given as a tuple/list of [start, end] locations. Can be used for any beam under shear, really.
    :param forces: iterable(force)
    :param locations: iterable(start[, stop] location)
    :return: numpy.ndarray(sheary)
    """
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
    """
    Same thing as shearxInternal, but for shear in x along fuselage in z. Is not (currently) complete due to missing
    reaction forces.
    :return: numpy.ndarray(shearx)
    """
    shearx = parameters.Sx
    shearx = shearx + numpy.where(z_disc>parameters.Lf02, parameters.Sx2, 0)
    shearx = shearx + numpy.where(z_disc>parameters.Lf01, parameters.Sx1, 0)

    return shearx


def momentxInternal():
    """
    This method will determine the internal moment force x along the beam in z. Takes into account the distributed load
    and the 2 reaction forces.
    :return: numpy.ndarray(Mx)
    """
    momentx = - parameters.q * z_disc * z_disc/2
    momentx = momentx + numpy.where(z_disc>parameters.Lf02, parameters.Sy2 * (z_disc - parameters.Lf02), 0.)
    momentx = momentx + numpy.where(z_disc>parameters.Lf01, parameters.Sy1 * (z_disc - parameters.Lf01), 0.)

    assert abs(momentx[0] < 0.001), 'Internal moment_x at aft of fuselage is not 0.'
    assert abs(momentx[-1] + parameters.q * step_size * step_size / 2) < 0.001, 'Internal moment_x at nose of fuselage is not 0.'

    return momentx


def momentxInternal_alt(forces, locations):
    """
    Same thing as shearyInternal_alt, but for momentxInternal: a general method to determine the internal moment.
    :param forces: iterable(forces)
    :param locations: locations: iterable(start[, stop] location)
    :return:
    """
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
    # Determines moment in z, but (currently) incomplete due to not having all reaction forces. Currently is a constant.
    momentz = z_disc * 0.
    momentz += parameters.M3 + parameters.Sx*parameters.R  # this is taken from equation 2.4; momentz = Mz,A - Mz,B
    return momentz


sy = shearyInternal_alt(parameters.external_forces_y, parameters.external_locations_y)
pass
