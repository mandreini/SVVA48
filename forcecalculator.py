import parameters

#MOMENTS AROUND BACK LANDING GEAR LOOKING FROM THE SIDE
Fy = (3*parameters.W/parameters.Lf2)*(parameters.Lf1 + parameters.Lf2 - parameters.L/2)
By = parameters.W - Fy

#MOMENTS AROUND RIGHT LANDING GEAR
Ry = -Fy/2 + parameters.W/2 + parameters.Sx*(parameters.dlgy+parameters.dtaily)/parameters.Lf3
Ly = By - Ry

#MOMENTS AROUND BACK LANDING GEAR LOOKING FROM TOP
Fx = -parameters.Sx*(parameters.L-parameters.Lf1-parameters.Lf2+parameters.dtailz)/parameters.Lf2

#MOMENTS AROUND FRONT LANDING GEAR LOOKING FROM TOP
Bx = parameters.Sx*(parameters.L - parameters.Lf1 + parameters.dtailz)/parameters.Lf2


print Fx,Bx,parameters.Sx