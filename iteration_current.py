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
        #sigma = +bending_stress(Mx,boom_array[n][2],Ixx) +bending_stress(My,boom_array[n][1],Iyy)
        sigma = bending_stress_precise(Mx,My,boom_array[n][1],boom_array[n][2],Ixx,Iyy,Ixy)
        #sigma += 0.0  # adding flat normal stress as assuming that stress will never be 0 at any point
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
    floor_width = 2*numpy.sqrt(parameters.R**2-(parameters.R-parameters.hf)**2)
    floor_stringer_area = parameters.tf*floor_width/len(floor_booms)
    fuselage_stringer_area = parameters.tst * (parameters.wst + parameters.hst)
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
    cell_I_fuselage = []
    cell_II_fuselage = []
    cell_floor = []

    for i in range(len(all_booms)):
        if all_booms[i][2] > -(parameters.R -parameters.hf):
            cell_I_fuselage.append(all_booms[i])

        if all_booms[i][2] < -(parameters.R -parameters.hf):
            cell_II_fuselage.append(all_booms[i])

        if all_booms[i][2] == -(parameters.R -parameters.hf):
            cell_floor.append(all_booms[i])
    last_two = cell_I_fuselage[-2:]
    del cell_I_fuselage[-2:]
    cell_I_fuselage = numpy.append(cell_I_fuselage,cell_floor,axis=0)
    cell_I_fuselage = numpy.append(cell_I_fuselage,last_two,axis=0)

    return numpy.array(cell_I_fuselage),numpy.append(cell_II_fuselage,numpy.flipud(cell_floor),axis=0)

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

centroid = beamProperties.centroidCalculation(all_booms)
Ixx_booms,Iyy_booms,Ixy_booms = idealized_structure_moment_of_inertia(all_booms,centroid)

#ITERATION
for i in range(20):

    #Calculate idealized structure centroid
    centroid = beamProperties.centroidCalculation(all_booms)

    #Calculate idealized structure moments of inertia
    Ixx_booms,Iyy_booms,Ixy_booms = idealized_structure_moment_of_inertia(all_booms,centroid)

    #Calculate stresses in new idealized structure
    stresses = normal_stress_calculation(all_booms,Ixx_booms,Iyy_booms,Ixy_booms,1,-3.3e6)

    #Recalculate areas for booms
    fuselage_booms[:,0],floor_booms[:,0] = recalculateBoomAreaFromStresses(stresses,fuselage_booms,floor_booms)

    epsilon = numpy.mean(fuselage_booms[:,0]/all_booms[0:len(fuselage_booms),0])-1
    epsilon2 = numpy.mean(floor_booms[:,0]/all_booms[-len(floor_booms):,0])-1

    all_booms = numpy.append(fuselage_booms,floor_booms,axis=0)

cell1, cell2 = divide_into_multicell(all_booms)

#Properties cell1
centroid_cell1 = beamProperties.centroidCalculation(cell1)
Ixx_cell1, Iyy_cell1, Ixy_cell1 = idealized_structure_moment_of_inertia(cell1,centroid_cell1)

#Properties cell2
centroid_cell2 = beamProperties.centroidCalculation(cell2)
Ixx_cell2, Iyy_cell2, Ixy_cell2 = idealized_structure_moment_of_inertia(cell2,centroid_cell2)



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


def angle_between(p1, p2):
    ang1 = numpy.arctan2(*p1[::-1])
    ang2 = numpy.arctan2(*p2[::-1])
    return numpy.rad2deg((ang1 - ang2) % (2 * numpy.pi))

def Aqb_term_calculator(cell,q_fuselage):
    Aqb = []

    for i in range(len(cell)):
        if i + 1 == len(cell):
            p1=(cell[i][1],cell[i][2])
            p2=(cell[0][1],cell[0][2])
        else:
            p1=(cell[i][1],cell[i][2])
            p2=(cell[i+1][1],cell[i+1][2])

        A = angle_between(p1,p2)*parameters.R**2/2
        Aqb.append(A*q_fuselage[i])


        return sum(Aqb)

def lengths_calculator(booms):
    upper_booms = booms[booms[:,2]> -parameters.R+parameters.hf]
    lower_booms = booms[booms[:,2]< -parameters.R+parameters.hf]
    floor_booms = booms[booms[:,2] == -parameters.R+parameters.hf]
    Lu = 0
    Lf = 0
    Ld = 0

    for i in range(len(upper_booms)-1):
        Lu = Lu + numpy.sqrt((upper_booms[i+1][1]-upper_booms[i][1])**2+(upper_booms[i+1][2]-upper_booms[i][2])**2)

    for i in range(len(lower_booms)-1):
        Ld = Ld + numpy.sqrt((lower_booms[i+1][1]-lower_booms[i][1])**2+(lower_booms[i+1][2]-lower_booms[i][2])**2)

    for i in range(len(floor_booms)-1):
        Lf = Lf + numpy.sqrt((floor_booms[i+1][1]-floor_booms[i][1])**2+(floor_booms[i+1][2]-floor_booms[i][2])**2)

    return Lu,Ld,Lf

length_upper,length_down,length_floor = lengths_calculator(all_booms)




'''
fuselage_booms = numpy.concatenate((fuselage_booms,[fuselage_booms[0]]),axis=0)
plt.ylabel('y (m)')
plt.xlabel('x (m)')

plt.scatter(fuselage_booms[:,1],fuselage_booms[:,2])
plt.plot(fuselage_booms[:,1],fuselage_booms[:,2])

plt.scatter(floor_booms[:,1],floor_booms[:,2])
plt.plot(floor_booms[:,1],floor_booms[:,2])
plt.grid(True)
plt.show()


def calculate_shear_flow_cell1(booms,Ixx,Iyy,centroid,Sx,Sy,cut):
    booms = numpy.roll(numpy.flipud(booms),-cut,axis=0)
    Aiyi_sum = []
    Aixi_sum = []
    i_list = []
    q_unit = []
    f_count_x = []
    f_count_y = []

    for i in range(len(booms)):

        Aixi = (booms[i][1]-centroid[0])*booms[i][0]
        Aiyi = (booms[i][2]-centroid[1])*booms[i][0]

        Aixi_sum.append(Aixi)
        Aiyi_sum.append(Aiyi)

        q = -Sx/Iyy*sum(Aixi_sum)-Sy/Ixx*sum(Aiyi_sum)

        if i == (len(booms)):
            q = 0

        if i + 1 == (len(booms)):
            y_dist = numpy.absolute(booms[i][2]-booms[0][2])
            x_dist = numpy.absolute(booms[i][1]-booms[0][1])
        else:
            y_dist = numpy.absolute(booms[i][2]-booms[i+1][2])
            x_dist = numpy.absolute(booms[i][1]-booms[i+1][1])

        if ((booms[i][1]>=0) & (booms[i][2]>=0)):
            f_count_y.append(q*y_dist)
            f_count_x.append(-q*x_dist)
        elif ((booms[i][1]<0) & (booms[i][2]>=0)):
            f_count_y.append(-q*y_dist)
            f_count_x.append(-q*x_dist)
        elif ((booms[i][1]<0) & (booms[i][2]<0)):
            f_count_y.append(-q*y_dist)
            f_count_x.append(q*x_dist)
        elif ((booms[i][1]>=0) & (booms[i][2]<0)):
            f_count_y.append(q*y_dist)
            f_count_x.append(q*x_dist)

        q_unit.append(q)

    return q_unit,f_count_x,f_count_y,booms

def calculate_shear_flow_cell2(booms,Ixx,Iyy,centroid,Sx,Sy,cut):
    booms = numpy.roll(numpy.flipud(booms),-cut,axis=0)

    Aiyi_sum = []
    Aixi_sum = []
    i_list = []
    q_unit = []
    f_count_x = []
    f_count_y = []

    for i in range(len(booms)):

        Aixi = (booms[i][1]-centroid[0])*booms[i][0]
        Aiyi = (booms[i][2]-centroid[1])*booms[i][0]

        Aixi_sum.append(Aixi)
        Aiyi_sum.append(Aiyi)

        q = -Sx/Iyy*sum(Aixi_sum)-Sy/Ixx*sum(Aiyi_sum)

        if i == (len(booms)):
            q = 0

        if i + 1 == (len(booms)):

            y_dist = numpy.absolute(booms[i][2]-booms[0][2])
            x_dist = numpy.absolute(booms[i][1]-booms[0][1])
        else:

            y_dist = numpy.absolute(booms[i][2]-booms[i+1][2])
            x_dist = numpy.absolute(booms[i][1]-booms[i+1][1])

        if ((booms[i][2]) == -(parameters.R -parameters.hf)):
            if ((booms[i+1][2]) != -(parameters.R -parameters.hf)):
                f_count_y.append(-q*y_dist)
                f_count_x.append(-q*x_dist)
            else:
                f_count_y.append(q*y_dist)
                f_count_x.append(q*x_dist)
        elif ((booms[i][1]<=0) & (booms[i][2]<0)):
            f_count_y.append(q*y_dist)
            f_count_x.append(-q*x_dist)
        elif ((booms[i][1]==0) & (booms[i][2]<0)):
            f_count_y.append(-q*y_dist)
            f_count_x.append(-q*x_dist)

        q_unit.append(q)

    return q_unit,f_count_x,f_count_y,booms

'''
def get_sum_floor_q(cell1,q1,cell2,q2):


    fuselage1 = boom_flows1[boom_flows1[:,2]!=-parameters.R+parameters.hf]
    fuselage2 = boom_flows2[boom_flows2[:,2]!=-parameters.R+parameters.hf]
    fuselage2 = np.roll(np.fuselage2,)
    #del boom_flows1[-2:]

    fuselage = numpy.append(boom_flows1,boom_flows2)

    boom_flows1 = numpy.column_stack((cell1,q1))
    boom_flows2 = numpy.column_stack((cell2,q2))

    q_floor_1 = boom_flows1[boom_flows1[:,2]==-parameters.R+parameters.hf][3]
    q_floor_2 = boom_flows2[boom_flows2[:,2]==-parameters.R+parameters.hf][3]
    q_floor = q_floor_1 + q_floor_2
    floor = boom_flows1[boom_flows1[:,2]==-parameters.R+parameters.hf]
    floor = numpy.delete(floor,1,0)
    floor[:,3]=q_floor

    return floor,boom_flows1,boom_flows2,q_floor_2


def go_around_sum_lengths(reordered_booms,q):

    booms_fuselage = reordered_booms[reordered_booms[:,2] != -parameters.R+parameters.hf]
    booms_floor = reordered_booms[reordered_booms[:,2] == -parameters.R+parameters.hf]

    sum_L_fuselage = []
    sum_L_floor = []

    for i in range(len(booms_fuselage)-1):
        sum_L_fuselage.append(numpy.sqrt((booms_fuselage[i][1]-booms_fuselage[i+1][1])**2) + (booms_fuselage[i][2]-booms_fuselage[i+1][2])**2)
    for i in range(len(booms_floor)-1):
        sum_L_floor.append(numpy.sqrt((booms_floor[i][1]-booms_floor[i+1][1])**2) + (booms_floor[i][2]-booms_floor[i+1][2])**2)

    for i in len(reordered_booms):
        qb_L.append()
    return sum(sum_L_fuselage),sum(sum_L_floor)


def go_around_cell(reordered_booms,centroid,Sx,Sy,Ixx,Iyy,cell,cut):
    q_unit_fueselage = []
    new_reordered_booms = numpy.empty([len(reordered_booms),len(reordered_booms[0])])
    f_count_x = []
    f_count_y = []

    AiBi_x_sum = []
    AiBi_y_sum = []

    for i in range(len(reordered_booms)):
        i = i + cut + 1
        AiBi_x = (reordered_booms[i][1]-centroid[0])*reordered_booms[i][0]
        AiBi_y = (reordered_booms[i][2]-centroid[1])*reordered_booms[i][0]

        AiBi_x_sum.append(AiBi_x)
        AiBi_y_sum.append(AiBi_y)

        q = -Sx/Iyy*sum(AiBi_x_sum)-Sy/Ixx*sum(AiBi_y_sum)

        if i == (len(reordered_booms)+cut):

            q = 0

        q_unit_fueselage.append(q)

        if i+1 == len(reordered_booms):
            y_dist = numpy.absolute(reordered_booms[i][2]-reordered_booms[0][2])
            x_dist = numpy.absolute(reordered_booms[i][1]-reordered_booms[0][1])
        else:
            y_dist = numpy.absolute(reordered_booms[i][2]-reordered_booms[i+1][2])
            x_dist = numpy.absolute(reordered_booms[i][1]-reordered_booms[i+1][1])

        if cell == 1:
            if ((reordered_booms[i][1]>=0) & (reordered_booms[i][2]>=0)):
                f_count_y.append(q*y_dist)
                f_count_x.append(-q*x_dist)
            elif ((reordered_booms[i][1]<0) & (reordered_booms[i][2]>=0)):
                f_count_y.append(-q*y_dist)
                f_count_x.append(-q*x_dist)
            elif ((reordered_booms[i][1]<0) & (reordered_booms[i][2]<0)):
                f_count_y.append(-q*y_dist)
                f_count_x.append(q*x_dist)
            elif ((reordered_booms[i][1]>=0) & (reordered_booms[i][2]<0)):
                f_count_y.append(q*y_dist)
                f_count_x.append(q*x_dist)

        if cell == 2:
            if ((reordered_booms[i][2]) == -(parameters.R -parameters.hf)):

                if ((reordered_booms[i+1][2]) == -(parameters.R -parameters.hf)):
                    f_count_y.append(q*y_dist)
                    f_count_x.append(-q*x_dist)
                else:
                    f_count_y.append(-q*y_dist)
                    f_count_x.append(q*x_dist)
            elif ((reordered_booms[i][1]<0) & (reordered_booms[i][2]<0)):
                f_count_y.append(-q*y_dist)
                f_count_x.append(q*x_dist)
            elif ((reordered_booms[i][1]>=0) & (reordered_booms[i][2]<0)):
                f_count_y.append(q*y_dist)
                f_count_x.append(+q*x_dist)
        new_reordered_booms[i] = reordered_booms[i+cut+1]
    return q_unit_fueselage,f_count_x,f_count_y,new_reordered_booms


def calculate_opensection_shear_flows(booms,Sx,Sy,Ixx,Iyy,centroid,cut):

    #Getting booms for both fuselage and floor
    booms_fuselage = booms[booms[:,2] != -parameters.R+parameters.hf]
    booms_floor = booms[booms[:,2] == -parameters.R+parameters.hf]

    if (min(booms[:,2]) < (-parameters.R+parameters.hf)):

        turn_around_floor = booms_floor
        reordered_booms = numpy.append(turn_around_floor,booms_fuselage,axis=0)
        q_unit_fueselage,f_count_x,f_count_y,new_reordered_booms = go_around_cell(reordered_booms,centroid,Sx,Sy,Ixx,Iyy,2,cut)


    else:

        reordered_booms = numpy.append(booms_fuselage[0:(len(booms_fuselage)-2)],booms_floor,axis=0)
        reordered_booms = numpy.append(reordered_booms,booms_fuselage[len(booms_fuselage)-2:len(booms)],axis=0)

        q_unit_fueselage,f_count_x,f_count_y,new_reordered_booms = go_around_cell(reordered_booms,centroid,Sx,Sy,Ixx,Iyy,1,cut)


    return q_unit_fueselage,f_count_x,f_count_y,new_reordered_booms


q1,fx1,fy1,cell1 = calculate_opensection_shear_flows(cell1,1,1,Ixx_cell1,Iyy_cell1,centroid_cell1,-2)
q2,fx2,fy2,cell2 = calculate_opensection_shear_flows(cell2,1,1,Ixx_cell2,Iyy_cell2,centroid_cell2,-2)
print fx2,cell2


#print sum(q_floor[:,3]),cell2,cell1
print sum(fx1),sum(fx2),sum(fy1),sum(fy2)
'''
calculate_shear_flow_cell1(cell1,Ixx_booms,Iyy_booms,centroid,0,1,-2)
'''

#print go_around_sum_lengths(reordered_cell,q_fuselage)
'''





test_box = numpy.empty([8, 3])
test_box[:,0] = [1,1,1,1,1,1,1,1]
test_box[:,1] = [0.5,0.5,0,-0.5,-0.5,-0.5,0,0.5]
test_box[:,2] = [0,0.5,0.5,0.5,0,-0.5,-0.5,-0.5]

test_box_1,test_box_2 = test_box[0:5],test_box[[4,5,6,7,0]]
centroid_test = beamProperties.centroidCalculation(test_box)
Ixx_test, Iyy_test, Ixy_test = idealized_structure_moment_of_inertia(test_box,centroid_test)
q1,fx1,fy1,nvm1 = calculate_opensection_shear_flows(test_box_1,0,1,Ixx_test,Iyy_test,centroid_test,0)
q2,fx2,fy2,nvm2 = calculate_opensection_shear_flows(test_box_2,0,1,Ixx_test,Iyy_test,centroid_test,0)
print sum(fy1)+sum(fy2),sum(fx1)+sum(fx2)
print q1,q2




'''

'''


x_axis = []
for i in range(len(fy)):
    x_axis.append(i)

plt.scatter(x_axis,fy)
plt.grid(True)
plt.show()
'''
'''



#FUSELAGE PLOTTER

#fuselage_booms = numpy.concatenate((fuselage_booms,[fuselage_booms[0]]),axis=0)
#fuselage_booms[len(fuselage_booms)+1] = [1,2,3]
print fuselage_booms
plt.ylabel('y (m)')
plt.xlabel('x (m)')

plt.scatter(fuselage_booms[:,1],fuselage_booms[:,2])
plt.plot(fuselage_booms[:,1],fuselage_booms[:,2])

plt.scatter(floor_booms[:,1],floor_booms[:,2])
plt.plot(floor_booms[:,1],floor_booms[:,2])
plt.grid(True)
plt.show()






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





'''

