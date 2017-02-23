# header stuff

import parameters
import numpy


def recalculateBoomAreaFromStresses(stress_vals):

    sigma1, x_loc, y_loc = stress_vals[:,0], stress_vals[:,1], stress_vals[:,2]
    stringer_area = parameters.tst * (parameters.wst + parameters.hst)
    constant_factor = (parameters.ts*parameters.R*(2*(numpy.pi)/parameters.ns)/6)
    Bi = []
    for i in range(parameters.ns):

        if i + 1 == len(sigma1):
            i = i - len(sigma1)


        forward_stress_ratio = constant_factor*(2+sigma1[i+1]/sigma1[i])
        backward_stress_ratio = constant_factor*(2+sigma1[i-1]/sigma1[i])


        Bi.append(stringer_area+forward_stress_ratio+backward_stress_ratio)

    circle_boom_amount = len(Bi)

    for i in range(len(sigma1)-parameters.ns):

        if i == 0:
            forward_stress_ratio = constant_factor*(2+sigma1[circle_boom_amount+i+1]/sigma1[circle_boom_amount+i])
            backward_stress_ratio = constant_factor*(2+sigma1[20]/sigma1[circle_boom_amount+i])

            Bi.append(stringer_area+forward_stress_ratio+backward_stress_ratio)

        elif i == (len(sigma1)-parameters.ns-1):

            forward_stress_ratio = constant_factor*(2+sigma1[34]/sigma1[circle_boom_amount+i])
            backward_stress_ratio = constant_factor*(2+sigma1[circle_boom_amount+i-1]/sigma1[circle_boom_amount+i])

            Bi.append(stringer_area+forward_stress_ratio+backward_stress_ratio)

        else:
            forward_stress_ratio = constant_factor*(2+sigma1[circle_boom_amount+i+1]/sigma1[circle_boom_amount+i])
            backward_stress_ratio = constant_factor*(2+sigma1[circle_boom_amount+i-1]/sigma1[i])

            Bi.append(stringer_area+forward_stress_ratio+backward_stress_ratio)


    stress_vals[:,0] = numpy.array(Bi)

    return stress_vals



def boom_creator(n_booms,n_floor_booms):
    total_output = []
    A = parameters.tst * (parameters.wst + parameters.hst)

    for n in range(n_booms):
        boom_output = numpy.array([A,parameters.R*numpy.cos(n*2*numpy.pi/n_booms),parameters.R*numpy.sin(n*2*(numpy.pi/n_booms))])
        total_output.append(boom_output)


    total_output[0] = [A,parameters.R,0]
    total_output[8] = [A,0,parameters.R]
    total_output[17] = [A,-parameters.R,0]
    total_output[26] = [A,0,-parameters.R]

    floor_width = 2*numpy.sqrt(parameters.R**2-(parameters.R-parameters.hf)**2)

    for n in range(n_floor_booms):
        floor_booms = numpy.array([A,-floor_width/2+floor_width*(n+1)/(n_floor_booms+1),parameters.R-parameters.hf])
        total_output.append(floor_booms)

    return numpy.array(total_output)

total_output = boom_creator(36,2)


def bending_stress(Moment,y,I):
    return Moment*y/I


def normal_stress_calculation(boom_array,Ixx,Iyy,Mx,My):
    stresses = []
    for n in range(len(boom_array)):
        sigma = bending_stress(Mx,boom_array[n][2],Iyy)  + bending_stress(My,boom_array[n][1],Ixx)
        sigma += 0.0  # adding flat normal stress as assuming that stress will never be 0 at any point
        stresses.append([sigma,boom_array[n][1],boom_array[n][2]])

    return numpy.array(stresses)


def idealized_structure_moment_of_inertia(booms):
    areas, x_loc, y_loc = booms[:,0], booms[:,1], booms[:,2]
    Ixx = 0
    Iyy = 0

    for i in range(len(booms)):
        Ixx = Ixx + areas[i]*y_loc[i]**2
        Iyy = Iyy + areas[i]*x_loc[i]**2

    return Ixx,Iyy






s1 = normal_stress_calculation(total_output,1,1,200,200)

a1 = recalculateBoomAreaFromStresses(s1)
print idealized_structure_moment_of_inertia(a1)

