# This file will combine the bending/shear internal forces along with the cross-section shear/normal stress calculations

import time
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import parameters
import beamProperties
import iteration_current  # rename iteration_current to like stressDistribution.py
import internalForces


def main():
    """
    The process to be used is to first calculate the internal forces (bending and shear) along the fuselage. Then these
    values will be substituted for Mx, My, Sy, Sx in the cross-section stress distributions _for each location in z_.
    This will yield an n-dimensional array that can be used to determine all the idealisation and everything everywhere!
    :return: numpy.ndarray(results)
    """

    # beam0 = beamProperties.beam_properties  # the fuselage idealised using area-equivalence
    # booms = iteration_current.all_booms
    # beam0 = (iteration_current.centroid, iteration_current.Ixx_booms,
    #          iteration_current.Iyy_booms, iteration_current.Ixy_booms)
    # The 2 should be equal, but just in case, I don't want to have to restructure one to fit the other just yet...
    internal_forces = internalForces.all_forces  # all forces acting on the fuselage as a function of z
    z = internalForces.z_disc

    sigma_z = []
    tau_z = []
    boom_z = []

    for ind in range(len(z)):
        sheary_z = internal_forces['sheary'][ind]
        shearx_z = internal_forces['shearx'][ind]
        momentx_z = internal_forces['momentx'][ind]
        momenty_z = internal_forces['momenty'][ind]
        sigma_i, boom_i = iteration_current.normal_stress_calculator(momentx_z, momenty_z, (2,2))
        # sigma_i = iteration_current.normal_stress_calculation(booms, beam0[1], beam0[2], beam0[3], momentx_z, momenty_z)
        # boom_i = iteration_current.recalculateBoomAreaFromStresses(sigma_i, booms[:parameters.ns], booms[parameters.ns:])
        tau_i = iteration_current.shear_stress_calculator(boom_i, shearx_z, sheary_z, parameters.M3)  # for shear stress calculations

        sigma_z.append(sigma_i)
        boom_z.append(boom_i)
        tau_z.append(tau_i)

    # rearrange the sigma, areas and shear vals into a more readable format: val in x,y along z
    sigma = []
    boom = []
    tau = []

    for ind in range(len(sigma_z)):
        sigma.append([i[0] for i in sigma_z[ind]])
        boom.append([i[0] for i in boom_z[ind]])
        # tau.append([i[0] for i in boom_z[ind]])  # placeholder, sue me
        tau.append([i[3] for i in tau_z[ind]])

    x_locs = [i[1] for i in sigma_z[0]]
    y_locs = [i[2] for i in sigma_z[0]]

    a = lambda x: numpy.array(x)

    return a(sigma), a(boom), a(tau), a(x_locs), a(y_locs)


def plot_vals(z, locs, val, show_plot=True):
    """
    This function will use matplotlib to show the 3-D model for the fuselage in order to show the values calculated
    in a practical format. Also will look snazzy on the final report.
    val - 2-D array. axis=0 will be the values for area/normal-stress/shear-stress and have the length of the number of booms
    axis=1 will have the length of z_disc and have the length of the z discretisation
    :param z: (constant) z-locations
    :param val: 2-D array of value for the booms (len_booms, len_z)
    :return: something, idk
    """

    v2 = numpy.array(val, dtype=float)
    # v2 = v2/numpy.max(v2)
    v2 = v2.swapaxes(0, 1)
    v3 = [v2 for v in v2]  # fuck matplotlib for the 20th time i swear
    points = numpy.meshgrid(locs[0], locs[1], z)
    x_v, y_v, z_v = points
    y_v = numpy.swapaxes(y_v, 0, 1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('Von Mises Stress Distribution\n')
    ax.set_xlabel('\nx')
    ax.set_ylabel('\nz')
    ax.set_zlabel('\ny')
    # ax.set_xlim3d(-4, 4)
    # ax.set_ylim3d(-4, 4)
    # ax.set_zlim3d(-1, 4)
    ax.scatter(x_v, z_v, zs=y_v, c=v3)
    # yes im switching the y and z axes around because i cant rotate the figure since it lags to high heaven and I'm out of patience for matplotlib
    plt.savefig('results.png')
    # testing a circle
    # t = numpy.linspace(0, 2*numpy.pi)
    # x2 = 3.1*numpy.cos(t)
    # y2 = 3.1*numpy.sin(t)
    # ax.scatter(x2, y2)

    if show_plot:
        plt.show()

teset = False
if teset:
    sample_vals = [[50, 75, 100], [150, 175, 200], [250, 275, 300], [350, 375, 400]]
    z = [0, 1, 2, 3]
    x = [3.1, -3.1*numpy.cos(numpy.pi/3), -3.1*numpy.cos(numpy.pi/3)]
    y = [0., 3.1*numpy.sin(numpy.pi/3), -3.1*numpy.sin(numpy.pi/3)]

    plot_vals(z, (x, y), sample_vals)


def detImportantVals(z, loc, var):
    """
    This function will determine the maximum stress (shear and normal) and its location from the everything
    :param z: z-loc array (1-D)
    :param loc: x,y coordinates (2-D)
    :param var: value of the variable at the x,y locations (2-D)
    :return:
    """
    location_of_max = numpy.where(var==numpy.max(var))
    max_var = var[location_of_max]
    z_loc = z[location_of_max[0]]
    max_info = max_var[0], z_loc[0], location_of_max[1][0]

    location_of_min = numpy.where(var==numpy.min(var))
    min_var = var[location_of_min]
    z_loc = z[location_of_min[0]]
    min_info = min_var[0], z_loc[0], location_of_min[1][0]

    return max_info, min_info

von_mises = lambda n, s: numpy.sqrt(n*n + 1000./3 * s*s)
t0 = time.time()
vals = main()
dt = time.time() - t0
print 'time taken: %f' % round(dt, 3)
normal, area, shear, x, y = vals
total_stress = von_mises(normal, shear)
if True: plot_vals(internalForces.z_disc, (x, y), total_stress, show_plot=False)
maximum, minimum = detImportantVals(internalForces.z_disc, (x, y), total_stress)
# normal: 2-D array; len(normal) = 70 (z discretisation); len(normal[0]) = 41 (amount of booms). same for area and shear

if False:
    color = total_stress[34]
    plt.scatter(x, y, s=50, c=color)
    plt.title('von Mises stress distribution at point of maximum stress (z=34)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(b=True)
    plt.savefig('z34stress.png')
    plt.show()

pass
