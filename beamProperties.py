# blah blah heading

import numpy
import parameters

def _determineBooms(debug=parameters.debug):
    # used to determine booms properties

    # skin booms
    delta_theta = 2*numpy.pi / parameters.ns
    skin_boom_spacing = delta_theta * parameters.R
    stringer_area = parameters.tst * (parameters.wst + parameters.hst)
    skin_boom_area = skin_boom_spacing*parameters.ts +  stringer_area

    skin_boom_x = parameters.R * numpy.cos(delta_theta * numpy.arange(0, parameters.ns, 1))
    skin_boom_y = parameters.R * numpy.sin(delta_theta * numpy.arange(0, parameters.ns, 1))

    # floor booms
    floor_width = 2 * numpy.sqrt((2*parameters.hf*parameters.R - parameters.hf**2))
    n_booms_floor = int(floor_width / skin_boom_spacing)
    floor_boom_spacing = floor_width / n_booms_floor
    area_floor = floor_width * parameters.tf
    floor_boom_area = area_floor / n_booms_floor

    offset = - floor_width/2 + floor_boom_spacing/2
    floor_boom_x = numpy.array([i*floor_boom_spacing for i in range(n_booms_floor)]) + offset
    floor_boom_y = numpy.array([-parameters.R+parameters.hf]*int(n_booms_floor))

    assert isinstance(n_booms_floor, int), 'Amount of floor booms is not an integer.'
    assert n_booms_floor * floor_boom_spacing < floor_width, 'Spacing of floor booms is too large.'
    assert (skin_boom_x**2 + skin_boom_y**2 - parameters.R**2 < 0.001).all(), 'Skin boom lies (roughly) on the fuselage'
    assert sum(floor_boom_x) < 0.001, 'Floor booms are symmetrical'

    if debug: print("skin_boom_spacing: %f \nfloor_boom_spacing: %f" % (skin_boom_spacing, floor_boom_spacing))

    # return skin_boom_area, skin_boom_spacing, n_booms_floor, floor_boom_spacing, floor_boom_area
    return (skin_boom_area, skin_boom_spacing, skin_boom_x, skin_boom_y), \
           (floor_boom_area, floor_boom_spacing, floor_boom_x, floor_boom_y, n_booms_floor, )


def centroidCalculation(boom_vals, debug=parameters.debug):
    # determines location of the centroid of the structure.
    # Also used to check the accuracy of the structural idealisation.

    assert len(boom_vals) == 2 and len(boom_vals[0]) == 4 and len(boom_vals[1]) == 5, 'boom_vals argument is not correct format.'
    assert isinstance(boom_vals[0][0], float) and isinstance(boom_vals[0][1], float) and isinstance(boom_vals[1][0], float) and isinstance(boom_vals[1][1], float) and isinstance(boom_vals[1][4], int), 'boom areas or spacing not properly given'
    assert isinstance(boom_vals[0][2], numpy.array) and isinstance(boom_vals[0][3], numpy.array) and isinstance(boom_vals[1][2], numpy.array) and isinstance(boom_vals[1][3], numpy.array), 'boom locations not properly given'

    skin_booms, floor_booms = boom_vals

    x_bar = 0  # due to symmetry
    y_bar = (skin_booms[0]*sum(skin_booms[3]) + floor_booms[0]*sum(floor_booms[3])) / \
            (skin_booms[0]*len(skin_booms[3]) + floor_booms[0]*len(floor_booms[3]))

    assert abs(y_bar) < parameters.R, 'y_bar is not inside the cross section.'
    assert y_bar < 0, 'y_bar is negative. This is only valid for floor at negative y'
    assert (numpy.abs(skin_booms[2]) <= parameters.R).all(), 'skin booms are not location in the cross section.'

    if debug: print("skin_boom_y_locs: %f", ', '.join([str(i) for i in skin_booms[3]]))

    return x_bar, y_bar


def inertiaCalculations(boom_vals, centroid_val, debug=parameters.debug):
    # this will determine Ixx and Iyy, and takes Ixy as 0 due to symmetry in y-axis
    # only accounts for steiner terms in the booms!

    assert len(boom_vals) == 2 and len(boom_vals[0]) == 4 and len(boom_vals[1]) == 5, 'boom_vals argument is not correct format ([sk_area], [sk_spacing], [sk_x_vals], [sk_y_vals], [fl_area], [fl_spacing], [fl_x_vals], [fl_y_vals], [fl_number]).'
    assert isinstance(boom_vals[0][0], float) and isinstance(boom_vals[0][1], float) and isinstance(boom_vals[1][0],float) and isinstance(boom_vals[1][1], float), 'boom areas or spacing not properly given as floats.'
    assert isinstance(boom_vals[1][4], int), 'Number of booms in the floor is not given as an integer.'
    assert isinstance(boom_vals[0][2], numpy.array) and isinstance(boom_vals[0][3], numpy.array) and isinstance(boom_vals[1][2], numpy.array) and isinstance(boom_vals[1][3], numpy.array), 'boom locations not properly given as numpy.array objects.'
    assert len(centroid_val) == 2, 'centroid_val argument is not correct format (x_bar, y_bar).'
    assert isinstance(centroid_val[0], float) and isinstance(centroid_val[1], float), 'centroid values are not given as float datatypes.'

    skin_booms, floor_booms = boom_vals

    Ixx = sum(skin_booms[0]*(skin_booms[3]-centroid_val[1])**2) + sum(floor_booms[0]*(floor_booms[3]-centroid_val[1])**2)
    Iyy = sum(skin_booms[0]*(skin_booms[2]-centroid_val[0])**2) + sum(floor_booms[0]*(floor_booms[2]-centroid_val[0])**2)
    Ixy = 0.  # I am adding this in case it comes up later

    assert Ixx > 0, 'Error in Ixx, value is negative.'
    assert Iyy > 0, 'Error in Iyy, value is negative.'

    if debug: print("Ixx: %f \nIyy: %f" % (Ixx, Iyy))

    return Ixx, Iyy, Ixy

