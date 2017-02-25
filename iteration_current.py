# header stuff
import internalForces
import parameters
import numpy
import pygame
import beamProperties
import matplotlib.pyplot as plt


def boom_creator(n_booms,n_floor_booms):
    fuselage_booms = []
    floor_booms = []

    #Calculating fuselage booms initial area and location
    A_fuselage_booms = parameters.tst * (parameters.wst + parameters.hst) + parameters.R*(2*numpy.pi/36)*parameters.ts

    for n in range(n_booms):
        boom_output = numpy.array([A_fuselage_booms,parameters.R*numpy.cos(n*2*numpy.pi/n_booms),parameters.R*numpy.sin(n*2*(numpy.pi/n_booms))])
        fuselage_booms.append(boom_output)


    fuselage_booms[0] = [A_fuselage_booms,parameters.R,0]
    fuselage_booms[9] = [A_fuselage_booms,0,parameters.R]
    fuselage_booms[18] = [A_fuselage_booms,-parameters.R,0]
    fuselage_booms[27] = [A_fuselage_booms,0,-parameters.R]

    #Calculating floor booms initial area and location
    floor_width = 2*numpy.sqrt(parameters.R**2-(parameters.R-parameters.hf)**2)
    A_floor_booms = parameters.tf*floor_width/n_floor_booms
    for n in range(n_floor_booms):
        floor_boom = numpy.array([A_floor_booms,(-floor_width/2)+n*(floor_width)/(n_floor_booms-1),-(parameters.R-parameters.hf)])
        floor_booms.append(floor_boom)

    return numpy.array(fuselage_booms),numpy.array(floor_booms)


def bending_stress_precise(Mx,My,x,y,Ixx,Iyy,Ixy):
    return (Ixx*My-Ixy*Mx)*x/(Ixx*Iyy-Ixy**2)+(Iyy*Mx-Ixy*My)*y/(Ixx*Iyy-Ixy**2)

def neutral_axis_orientation(Mx,My,Ixx,Iyy,Ixy):
    return numpy.atan((Ixx*My-Ixy*Mx)/(Iyy*Mx-Ixy*My))*57.3

def bending_stress(Moment,y,I):
    return Moment*y/I

def normal_stress_calculation(boom_array,Ixx,Iyy,Ixy,Mx,My):
    stresses = []
   # na_orientation = neutral_axis_orientation(Mx,My,Ixx,Ixy,Ixy)
    for n in range(len(boom_array)):
        sigma = bending_stress_precise(Mx,My,boom_array[n][1],boom_array[n][2],Ixx,Iyy,Ixy)
        sigma += 0.0  # adding flat normal stress as assuming that stress will never be 0 at any point
        stresses.append([sigma,boom_array[n][1],boom_array[n][2]])

    return numpy.array(stresses)

def idealized_structure_moment_of_inertia(booms,centroid):
    areas, x_loc, y_loc = booms[:,0], booms[:,1], booms[:,2]
    Ixx = 0
    Iyy = 0
    Ixy = 0
    for i in range(len(booms)):
        Ixx = Ixx + areas[i]*(y_loc[i]-centroid[1])**2
        Iyy = Iyy + areas[i]*(x_loc[i]-centroid[0])**2
        Ixy = Ixy + areas[i]*(y_loc[i]-centroid[1])*(x_loc[i]-centroid[0])

    return Ixx,Iyy,Ixy


def recalculateBoomAreaFromStresses(stress_vals,fuselage_booms,floor_booms):

    sigma_fuselage, x_loc, y_loc = stress_vals[0:len(fuselage_booms),0], stress_vals[0:len(fuselage_booms),1], stress_vals[0:len(fuselage_booms),2]
    sigma_floor = stress_vals[-len(floor_booms):,0]
    floor_stringer_area = floor_booms[:,0]
    fuselage_stringer_area = fuselage_booms[:,0]
    constant_factor = (parameters.ts*parameters.R*(2*(numpy.pi)/parameters.ns)/6)
    Bi = []


    #RECALCULATING AREA OF FUSELAGE BOOMS

    forward_stress_ratio_fuselage = constant_factor*(2+numpy.roll(sigma_fuselage,-1)/sigma_fuselage)
    backward_stress_ratio_fueslage = constant_factor*(2+numpy.roll(sigma_fuselage,1)/sigma_fuselage)
    Bi_fuselage = forward_stress_ratio_fuselage + backward_stress_ratio_fueslage + fuselage_stringer_area

    #RECALCULATING AREA OF FLOOR BOOMS(FIRST AND LAST BOOMS DONE MANUALLY)

    #FIRST
    forward_stress_ratio_first = constant_factor*(2+sigma_fuselage[21]/sigma_floor[0])
    backward_stress_ratio_first = constant_factor*(2+sigma_fuselage[20]/sigma_floor[0])
    #MIDDLE
    forward_stress_ratio_floor = constant_factor*(2+numpy.roll(sigma_floor,-1)/sigma_floor)
    backward_stress_ratio_floor = constant_factor*(2+numpy.roll(sigma_floor,1)/sigma_floor)
    #LAST
    forward_stress_ratio_last = constant_factor*(2+sigma_fuselage[34]/sigma_floor[-1])
    backward_stress_ratio_last = constant_factor*(2+sigma_fuselage[33]/sigma_floor[-1])


    forward_stress_ratio_floor[0] = forward_stress_ratio_first
    forward_stress_ratio_floor[-1] = forward_stress_ratio_last

    backward_stress_ratio_floor[0] = backward_stress_ratio_first
    backward_stress_ratio_floor[-1] = backward_stress_ratio_last

    Bi_floor = forward_stress_ratio_floor + backward_stress_ratio_floor + floor_stringer_area

    return Bi_fuselage,Bi_floor


def divide_into_multicell(all_booms):
    areas, x_loc, y_loc = all_booms[:,0], all_booms[:,1], all_booms[:,2]
    cell_I = []
    cell_II = []

    for i in range(len(all_booms)):
        if all_booms[i][2] >= -(parameters.R -parameters.hf):
            cell_I.append(all_booms[i])

        if all_booms[i][2] <= -(parameters.R -parameters.hf):
            cell_II.append(all_booms[i])

    return numpy.array(cell_I),numpy.array(cell_II)

def create_force_matrix():
    shear_x = numpy.array(internalForces.shearxInternal())
    shear_y = numpy.array(internalForces.shearyInternal())
    moment_x = numpy.array(internalForces.momentxInternal())
    moment_z = numpy.array(internalForces.momentzInternal())
    moment_y = numpy.array(internalForces.momentyInternal())
    return numpy.column_stack((shear_x,shear_y,moment_x,moment_y,moment_z))



#Create booms
fuselage_booms,floor_booms = boom_creator(36,5)
all_booms = numpy.append(fuselage_booms,floor_booms,axis=0)
epsilon = 2
epsilon2 = 2





#ITERATION
while epsilon2>0.01 or epsilon>0.01:
    #Calculate idealized structure centroid
    centroid = beamProperties.centroidCalculation(all_booms)

    #Calculate idealized structure moments of inertia
    Ixx_booms,Iyy_booms,Ixy_booms = idealized_structure_moment_of_inertia(all_booms,centroid)

    #Calculate stresses in new idealized structure
    stresses = normal_stress_calculation(all_booms,Ixx_booms,Iyy_booms,Ixy_booms,70000000,20000000)

    #Recalculate areas for booms
    fuselage_booms[:,0],floor_booms[:,0] = recalculateBoomAreaFromStresses(stresses,fuselage_booms,floor_booms)

    epsilon = numpy.mean(fuselage_booms[:,0]/all_booms[0:len(fuselage_booms),0])-1
    epsilon2 = numpy.mean(floor_booms[:,0]/all_booms[-len(floor_booms):,0])-1

    all_booms = numpy.append(fuselage_booms,floor_booms,axis=0)


#FUSELAGE PLOTTER
plt.scatter(fuselage_booms[:,1],fuselage_booms[:,2])
plt.scatter(floor_booms[:,1],floor_booms[:,2])
plt.show()

cell1, cell2 = divide_into_multicell(all_booms)


def J_calculator(booms,fuselage_thickness,floor_thickness):

    #Getting booms for both fuselage and floor
    booms_fuselage = booms[booms[:,2] != -parameters.R+parameters.hf]
    booms_floor = booms[booms[:,2] == -parameters.R+parameters.hf]

    #Initalizing Js
    J_fuselage= 0
    J_floor = 0

    #Counting J for fuselage
    for i in range(len(booms_fuselage)-1):
        J_fuselage = J_fuselage + (fuselage_thickness**3)*numpy.sqrt((booms_fuselage[i,1]-booms_fuselage[i+1,1])**2 + (booms_fuselage[i,2]-booms_fuselage[i,2])**2 )

    #Counting J for floor
    for i in range(len(booms_floor)-1):
        J_floor = J_floor + (floor_thickness**3)*numpy.sqrt((booms_floor[i,1]-booms_floor[i+1,1])**2 + (booms_floor[i,2]-booms_floor[i,2])**2 )

    return J_fuselage+J_floor


def low_area_calculator(f_booms):

    #Getting area of triangle to calculate section area
    A_triangle = -f_booms[-1][2]*f_booms[-1][1]

    #Getting theta corresponding to the section that includes triangle and cell area
    theta = 2*numpy.arcsin(f_booms[-1][1]/parameters.R)

    #Getting cell area
    A_section = theta/(2*numpy.pi)*(parameters.R**2*numpy.pi)-A_triangle

    return A_section



def constant_shear_flow_calculator(cell_I,cell_II,torque):

    #Calculating Js for both cells
    J_1 = J_calculator(cell_I,parameters.ts,parameters.tf)
    J_2 = J_calculator(cell_II,parameters.ts,parameters.tf)

    #Ratio of T1/T2 to solve for T1 and T1
    ratio_T1_T2 = J_1/J_2
    T2 = torque/(ratio_T1_T2+1)
    T1 = torque - T2

    #Getting an array with only floor booms
    floor_booms = cell_I[cell_I[:,2]==-parameters.R+parameters.hf]

    #Calculating areas of both sections
    area_2 = low_area_calculator(floor_booms)
    area_1 = (parameters.R**2)*numpy.pi - area_2

    #Finding constant shear flows using T=2Aq
    constant_shear1 = T1/(2*area_1)
    constant_shear2 = T2/(2*area_2)

    return constant_shear1,constant_shear2


print constant_shear_flow_calculator(cell1,cell2,10)

'''
#FUSELAGE PLOTTER
forces = create_force_matrix()

#plot section
plt.subplot(221)
plt.title("Vy")
plt.plot(numpy.arange(parameters.L),forces[:,1])

plt.subplot(222)
plt.title("Mx")
plt.plot(numpy.arange(parameters.L),forces[:,2])

plt.subplot(223)
plt.title("Vx")
plt.plot(numpy.arange(parameters.L),forces[:,0])

plt.subplot(224)
plt.title("My")
plt.plot(numpy.arange(parameters.L),forces[:,3])

plt.show()




#TESTER FOR ALL QUANTITIES
x_axis = []
for i in range(len(stresses[:,0])):
    x_axis.append(i)

plt.scatter(x_axis,stresses[:,0])
plt.show()
#TESTER FOR ALL QUANTITITES



def calculate_shear_flows(booms,Sx,Sy,Ixx,Iyy):
    q_unit = []

    for i in range(parameters.ns):
        q_unit.append(-Sx/Iyy*booms[i][1]*booms[i][0]-Sy/Ixx*booms[i][2]*booms[i][0])


    return q_unit

'''

