# This will compute a more precise analysis of the shear through the open section frames.

import parameters

import matplotlib.pyplot as plt
import numpy


def plot_frame(f1, f2):
    """
    Do you seriously need a description for this?
    :param f1:
    :param f2:
    :return:
    """
    theta_test = numpy.linspace(0, 2 * numpy.pi, len(frame_1))

    plt.figure()
    plt.subplot(211)
    plt.plot(theta_test, frame_1, 'r')
    plt.title('frame1')
    plt.xlabel('theta')
    plt.ylabel('shear')

    plt.subplot(212)
    plt.plot(theta_test, frame_2, 'b')
    plt.title('frame2')
    plt.xlabel('theta')
    plt.ylabel('shear')

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()


def plot_circle(shear_flow, frame_no):
    """
    This will plot in a circle
    :param shear_flow: the q values to be plotted
    :return: a figure if you really want it
    """

    r = parameters.R
    theta_vals = numpy.linspace(0, 2.*numpy.pi, len(shear_flow))
    x = r * numpy.cos(theta_vals)
    y = r * numpy.sin(theta_vals)

    # set up arrays for plotting
    s2 = numpy.meshgrid(shear_flow, shear_flow)[0]
    locs = numpy.array([[x, y] for z in shear_flow])

    # set colormap to show the absolute value about the average, i.e. the symmetry of the shear distribution
    avg_shear = numpy.average(s2)
    s3 = s2 - avg_shear
    s3 = numpy.abs(s3)
    s3 = s3 + avg_shear

    # plot the shear distribution
    plt.scatter(locs[:,0], locs[:,1], c=s3)
    plt.title("Shear in frame %s" % frame_no)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(b=True)
    plt.xlim([-4, 4])
    plt.ylim([-4, 4])
    plt.savefig('shear_for_frame%s.png' % frame_no)
    # plt.show()

def getQ(theta):
    # blach blach uses elliptoids.
    alpha = numpy.pi/2 - theta
    r2 = parameters.R
    r1 = parameters.R - parameters.ts

    a2 = r2 * numpy.sin(alpha)
    b2 = r2 * (1 - numpy.cos(alpha))
    a1 = r1 * numpy.sin(alpha)
    b1 = r1 * (1 - numpy.cos(alpha))

    A1 = numpy.pi/2 * a1 * b1
    A2 = numpy.pi/2 * a2 * b2
    A  = A2 - A1

    y1 = 4./3*b1/numpy.pi
    y2 = 4./3*b2/numpy.pi

    ybar = (y2*A2 - y1*A1) / (A2 - A1)
    # note that ybar will be higher than y1 or y2 because math, this is inline with what 100% accuracy random internet sites tell me

    ybar = numpy.nan_to_num(ybar)
    A = numpy.nan_to_num(A)

    return ybar * A


def frame_analysis(Vy, Vx, M=0, debug=False):
    """
    This function will determine the shear in an open ring frame. The frame will assumed as a dz section of the skin that
    introduces the point loads (and torsion) into the fuselage. They will have the same radius and thickness as the skin.
    Approach is q=VyQyIxx - VxQxIyy. Where Qy, Qx is the first moment of area in y and x respectively. y_c, x_c and A
    will be determined with semi-elliptical areas
    :param V: float - shear force.
    :param M: floar - couple moment.
    :return: numpy.ndarray(shear flow values)
    """

    # section properties
    A = numpy.pi*parameters.R**2
    Ixx = numpy.pi*parameters.R**3*parameters.ts
    Iyy = Ixx  # symmetry
    Ixy = 0.  # symmetry

    # theta = numpy.linspace(0., numpy.pi/2.)
    # note that for some reason I seem to be starting with theta=0 at the bottom of the circle, not the x-axis.
    # I presume this is due to an inaccurate reference for the Q calc, but since it's an axis switch, I will do that.
    # Qquartervals = getQ(theta)
    # sheary = Vy*Qquartervals/Ixx
    # shearx = Vx*Qquartervals/Iyy
    #
    # circle = numpy.array([
    #     sheary - shearx[::-1],
    #     sheary[::-1] - shearx,
    #     sheary - shearx[::-1],
    #     sheary[::-1] - shearx
    # ]).flatten()
    #
    # circle += M / (2.*1/2 * numpy.pi*parameters.R**2)

    # determine shear with -Sx/Iy *int(t*x*ds); shear0 is closed section shear int(p*q*ds) + 2Aq_0 = 0; shearz is torsion
    theta = numpy.linspace(0, 2*numpy.pi, 200)

    sheary = -(Vy/Ixx)*parameters.ts*parameters.R*(-numpy.cos(theta)+1)
    shearx = -(Vx/Iyy)*parameters.ts*parameters.R*numpy.sin(theta)
    shear0 = -numpy.sum((shearx + sheary) * (2.*numpy.pi/200.)) / (2 * A)
    shearz = M / (2.*1/2 * numpy.pi*parameters.R**2)

    circle = shearx + sheary + shear0 + shearz

    if debug:
        t2 = numpy.linspace(0, 2*numpy.pi, len(circle))
        # plt.plot(t2, sheary)
        # plt.plot(t2, shearx[::-1])
        plt.plot(t2, circle)
        plt.show()

    return circle

frame_1 = frame_analysis(parameters.Sy1, parameters.Sx1, M=parameters.Mlg1)#, debug=True)
frame_2 = frame_analysis(parameters.Sy2, parameters.Sx2, M=parameters.Mlg2)#, debug=True)  # check if there is an additional torque at rear landing gear
# plot_frame(frame_1, frame_2)
plot_circle(frame_1, '1')
plot_circle(frame_2, '2')
'''
brief analysis:
both plots have the smooth-transition and rough-transition at the same point. However, the shape and magnitude of the
shear flow in the frames is pretty different, which shouldnt be surprising given the reaction forces.
'''