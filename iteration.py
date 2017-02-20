# header stuff

import parameters
import numpy

def recalculateBoomAreaFromStresses(stress_vals):
    """
    B1 = tD*b*(2+sigma2/sigma1)
    B36 = f(sigma1/sigma36)
    :param stress_vals:
    :return:
    """
    sigma1, x_loc, y_loc = stress_vals[:,0], stress_vals[:,1], stress_vals[:,2]
    b = numpy.sqrt((numpy.roll(x_loc, -1) - x_loc)**2 + (numpy.roll(y_loc, -1)-y_loc)**2)  # refactor to calculate arc length, not absolute distance
    sigma2 = numpy.roll(sigma1, -1)
    sigma0 = numpy.roll(sigma1, 1)
    Bi = parameters.ts*b/6 * ((2+sigma2/sigma1)+(2+sigma0/sigma1))

    return numpy.column_stack((Bi, x_loc, y_loc))

#filip's stress calculations
R = 3.1
A = 0.0022
n = 0
total_output = []
for n in range(36):
    boom_output = numpy.array([A,R*numpy.cos(n*2*numpy.pi/36),R*numpy.sin(n*2*(numpy.pi/36))])
    total_output.append(boom_output)

total_output = numpy.array(total_output)

def bending_stress(Moment,y,I):
    return Moment*y/I

def normal_stress_calculation(boom_array,Ixx,Iyy,Mx,My):
    stresses = []
    for n in range(36):
        sigma = bending_stress(Mx,boom_array[n][2],Iyy) + bending_stress(My,boom_array[n][1],Ixx)
        stresses.append([sigma,boom_array[n][1],boom_array[n][2]])

    return numpy.array(stresses)

s1 = normal_stress_calculation(total_output,5,1,2,0)
a1 = recalculateBoomAreaFromStresses(s1)
s2 = normal_stress_calculation(a1, 5, 1, 2, 0)
a2 = recalculateBoomAreaFromStresses(s1)

pass
