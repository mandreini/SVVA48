# header stuff

import parameters
import numpy

def recalculateBoomAreaFromStresses(stress_vals):
    """
    B1 = tD*b21*(2+sigma2/sigma1) + tD*b01*(2+sigma0/sigma1)
    tD - thickness, constant. b21, y-distance from boom2-boom1; b01, y-distance from boom0-boom1.
    sigma0/1/2, stresses in boom 0/1/2 respectively.
    :param stress_vals: numpy.ndarry([numpy.ndarray(stress in booms), x_location, y_location])
    :return: numpy.ndarray([numpy.ndarray(area of booms), x_location, y_location])
    """
    sigma1, x_loc, y_loc = stress_vals[:,0], stress_vals[:,1], stress_vals[:,2]

    b = numpy.sqrt((x_loc-numpy.roll(x_loc, 1))**2 + (y_loc-numpy.roll(y_loc, 1))**2)
    sigma2 = numpy.roll(sigma1, -1)
    sigma0 = numpy.roll(sigma1, 1)
    Bi = parameters.ts/6.*b*(2+sigma0/sigma1) + parameters.ts/6.*b*(2+sigma2/sigma1)

    if not numpy.isfinite(Bi).all(): Bi[numpy.where(not Bi)] = 0. # replace anything that's not a number with 0
    return numpy.column_stack((Bi, x_loc, y_loc))
    # I believe I have fixed the issue. The problem was that there wasn't one, the first iteration with a stress of 5
    # Simply didn't give results that changed anything. 0.00215 becomes 0.002148. If it doesn't show a noticible change
    # With numbers that aren't completely basic, I will be rather annoyed...


#filip's stress calculations - should be refactored to numpy calculations at some point, probably
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
        sigma = bending_stress(Mx,boom_array[n][2],Iyy) # + bending_stress(My,boom_array[n][1],Ixx)
        sigma += 0.1  # adding flat normal stress as assuming that stress will never be 0 at any point
        stresses.append([sigma,boom_array[n][1],boom_array[n][2]])

    return numpy.array(stresses)

s1 = normal_stress_calculation(total_output,500,1,2,0)
a1 = recalculateBoomAreaFromStresses(s1)

# TODO: add the iteration to optimize structure
# TODO: include all reaction forces to determine all internal stresses to correctly compute boom values
# TODO: put internal shear/normal stresses in presentable format
# TODO: determine shear in frames
# TODO: get a passing grade

pass
