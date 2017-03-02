__author__ = 'Matt'
__date__ = '18-Feb-17'

import numpy
import parameters


def _determineBooms1(debug=parameters.debug):
    """
    This method will determine the initial area of the booms by principle of area equivalence. The floor area component
    will be added to the (2) booms closest to it. (This is done because y-const in floor for Mx. Floor is symmetrical in
    x to allow for idealisation in that particular case. Only valid for this type of fuselage).
    :param debug: bool - prints out relevant values if True
    :return: numpy.ndarray([(B1area, B1x, B1y), (B2area, B2x, B2y), ...])
    """

    # skin booms
    delta_theta = 2*numpy.pi / parameters.ns
    skin_boom_spacing = delta_theta * parameters.R
    stringer_area = parameters.tst * (parameters.wst + parameters.hst)
    skin_boom_area = numpy.array([skin_boom_spacing*parameters.ts + stringer_area]*parameters.ns)

    skin_boom_x = parameters.R * numpy.cos(delta_theta * numpy.arange(0, parameters.ns, 1))
    skin_boom_y = parameters.R * numpy.sin(delta_theta * numpy.arange(0, parameters.ns, 1))

    # add floor component
    floor_y_loc = -parameters.R+parameters.hf
    floor_width = 2 * numpy.sqrt((2 * parameters.hf * parameters.R - parameters.hf ** 2))
    area_floor = floor_width * parameters.tf

    floor_boom_loc_1 = numpy.where(min(abs(skin_boom_y-floor_y_loc)) == abs(skin_boom_y-floor_y_loc))[0]  # locate closest boom to floor
    floor_boom_loc_2 = parameters.ns/2 + (parameters.ns-floor_boom_loc_1)  # find other boom through symmetry
    # add line to determine error between floor and floor boom
    skin_boom_area[floor_boom_loc_1] += area_floor/2
    skin_boom_area[floor_boom_loc_2] += area_floor/2

    skin_boom_data = numpy.column_stack((skin_boom_area, skin_boom_x, skin_boom_y))

    assert (skin_boom_x**2 + skin_boom_y**2 - parameters.R**2 < 0.001).all(), 'Skin booms lie (roughly) on the fuselage'
    assert abs(sum(skin_boom_data[:,[0]]) - (area_floor + stringer_area*parameters.ns + 2 * numpy.pi * parameters.R * parameters.ts)) < 0.001, 'Area is not equivalent.'

    return skin_boom_data

    # old code when fuselage was discretised differently with multiple booms in floor
    # floor booms
    # floor_width = 2 * numpy.sqrt((2*parameters.hf*parameters.R - parameters.hf**2))
    # n_booms_floor = int(floor_width / skin_boom_spacing)
    # floor_boom_spacing = floor_width / n_booms_floor
    # area_floor = floor_width * parameters.tf
    # floor_boom_area = area_floor / n_booms_floor
    #
    # offset = - floor_width/2 + floor_boom_spacing/2
    # floor_boom_x = numpy.array([i*floor_boom_spacing for i in range(n_booms_floor)]) + offset
    # floor_boom_y = numpy.array([-parameters.R+parameters.hf]*int(n_booms_floor))
    # return (skin_boom_area, skin_boom_spacing, skin_boom_x, skin_boom_y), \
    #        (floor_boom_area, floor_boom_spacing, floor_boom_x, floor_boom_y, n_booms_floor,)

    # assert isinstance(n_booms_floor, int), 'Amount of floor booms is not an integer.'
    # assert n_booms_floor * floor_boom_spacing < floor_width, 'Spacing of floor booms is too large.'
    # assert sum(floor_boom_x) < 0.001, 'Floor booms are symmetrical'
    # assert function to check area equivalence

    # if debug: print("skin_boom_spacing: %f \nfloor_boom_spacing: %f" % (skin_boom_spacing, floor_boom_spacing))


def centroidCalculation(boom_vals, debug=parameters.debug):
    """
    This method will determine the centroid of the structural idealisation. It uses the formula y_bar = sum(Ay)/sum(A)
    :param boom_vals: numpy.ndarray([(B1area, B1x, B1y), (B2area, B2x, B2y), ...])
    :param debug: bool - prints out relevant values if True
    :return: tuple(x_bar, y_bar)
    """

    assert isinstance(boom_vals, numpy.ndarray) and len(boom_vals[0]) == 3, 'boom_vals is not correct format'

    x_bar = 0.  # due to symmetry
    y_bar = sum(boom_vals[:,0] * boom_vals[:,2]) / sum(boom_vals[:,0])

    assert abs(y_bar) < parameters.R, 'y_bar is not inside the cross section.'
    assert y_bar < 0, 'y_bar is negative. This is only valid for floor at negative y'
    assert (numpy.abs(boom_vals[2]) <= parameters.R).all(), 'skin booms are not location in the cross section.'

    if debug: print("skin_boom_y_locs: %f", ', '.join([str(i) for i in boom_vals[3]]))

    return x_bar, y_bar


def inertiaCalculations(boom_vals, centroid_val, debug=parameters.debug):
    """
    This method will determine the area moment of inertia in both xx, yy and xy. It (currently) assumes Ixy to be 0 due
    to symmetry about the y-axis. Only the steiner terms (Adx^2) of the booms are used in calculation
    :param boom_vals: numpy.ndarray([(B1area, B1x, B1y), (B2area, B2x, B2y), ...])
    :param centroid_val: tuple(x_bar, y_bar)
    :param debug: bool - prints out relevant values if True
    :return: tuple(Ixx, Ixy, Iyy)
    """

    assert isinstance(boom_vals, numpy.ndarray) and len(boom_vals[0]) == 3, 'boom_vals is not correct format'
    assert len(centroid_val) == 2, 'centroid_val argument is not correct format (x_bar, y_bar).'
    assert isinstance(centroid_val[0], float) and isinstance(centroid_val[1], float), 'centroid values are not given as float datatypes.'

    Ixx = sum(boom_vals[:,0]*(boom_vals[:,2]-centroid_val[1])**2)
    Iyy = sum(boom_vals[:,0]*(boom_vals[:,1]-centroid_val[0])**2)
    Ixy = 0.  # I am adding this in case it comes up later

    assert Ixx > 0, 'Error in Ixx, value is negative.'
    assert Iyy > 0, 'Error in Iyy, value is negative.'

    if debug: print("Ixx: %f \nIyy: %f" % (Ixx, Iyy))

    return Ixx, Iyy, Ixy

booms = _determineBooms1()
centroid = centroidCalculation(booms)
inertia = inertiaCalculations(booms, centroid)
pass
