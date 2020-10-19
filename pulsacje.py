import sys
import scipy.special as sci
import pandas as pd
import numpy as np
import math
import plotly.express as px

def welcome_msg():
    print("This programme generates a *.gif animation showing how star's surface pulsates in time")
    print()
    print("Usage: python3 pulsacje.py <f_in=> <f_out=> <t_total=> <delta_t=> <n_theta=> <n_phi=>")
    print("f_in     --> name of input text file containing modes parameters")
    print("f_out    --> name of the output file")
    print("t_total  --> total animation time expressed in days")
    print("delta_t  --> time step in the animation expressed in days")
    print("n_theta  --> number of discrete points on the star's astrophysical longitude")
    print("n_phi    --> number of discrete points on the star's astrophysical latitude")
    print()
    print("Content of the input text file:")
    print("frequency	amplitude	phase	lrank	mrank")
    print()

if len(sys.argv) != 7:
    welcome_msg()
    sys.exit(1)

for arg in sys.argv:
    if "f_in=" in arg:
        try:
            f_in = str(arg).replace("f_in=","")
        except ValueError:
            print("Wrong value of argument!")
        try:
            f = open(f_in)
            f.close()
        except FileNotFoundError:
            print("Your input file does not exist, try again!")
            sys.exit(1)
    if "f_out=" in arg:
        try:
            f_out = str(arg).replace("f_out=","")
        except ValueError:
            print("Wrong value of argument!")
    if "t_total=" in arg:
        try:
            t_total = float(str(arg).replace("t_total=",""))
        except ValueError:
            print("Wrong value of argument, float needed!")
    if "delta_t=" in arg:
        try:
            delta_t = float(str(arg).replace("delta_t=",""))
        except ValueError:
            print("Wrong value of argument, float needed!")
    if "n_theta=" in arg:
        try:
            n_theta = int(str(arg).replace("n_theta=",""))
        except ValueError:
            print("Wrong value of argument, integer needed!")
    if "n_phi=" in arg:
        try:
            n_phi = int(str(arg).replace("n_phi=",""))
        except ValueError:
            print("Wrong value of argument, integer needed!")

#function to calculate position vector
def position(f_in, t_total, delta_t, n_theta, n_phi):
    parameters = np.loadtxt(f_in)
    parameters = pd.DataFrame(parameters)
    num_lines = len(open(f_in).readlines())
   
    if num_lines == 1:
        parameters = parameters.T

    #quotient of t_total and delta_t will be the number of images in animation
    #the number needs to be integer 
    n_img = int(t_total / delta_t)
    
    #theta angle runs from 0 to pi
    #phi angle runs from 0 to 2pi 
    #by dividing those ranges by number of discrete points n_theta and n_phi we get delta_theta and delta_phi 
    delta_theta = np.pi / n_theta
    delta_phi = 2 * np.pi / n_phi

    cos_theta = []
    polynominal = []
    for i in range (1, n_theta + 2):
        cos1 = np.cos((i - 1) * delta_theta)
        cos_theta.append(cos1)

        for j in range (0, num_lines):
            y = sci.lpmv(parameters.iloc[j, 4], parameters.iloc[j, 3], cos_theta[i - 1])
            polynominal.append(y)

    cos_m_phi = []
    for i in range (1, n_phi + 2):
        phi_i = (i - 1) * delta_phi
        
        for j in range (0, num_lines):
            cos2 = np.cos(phi_i * parameters.iloc[j, 4])
            cos_m_phi.append(cos2)

    cos_ti = []
    for i in range (0, n_img):
        t_i = (i - 1) * delta_t
        arg = np.cos(2 * np.pi * parameters.iloc[0, 0] * t_i)
        cos_ti.append(arg)

    amplitude = []
    for i in range (0, num_lines):
        amplitude.append(parameters.iloc[0, 1])

    phi_i = []
    theta_i = []
    r = []
    for n in range (0, 3):
        for i in range(1, n_theta + 1):
            for j in range (1, n_phi + 1):
                r_i = 1
               
                for k in range (0, num_lines):
                    theta = (i - 1) * delta_theta
                    phi = (j - 1) * delta_phi

                    r_i = 1 + parameters.iloc[k, 1] * math.sqrt(((2 * parameters.iloc[k, 3] + 1) / (4 * np.pi)) *
                    ((math.factorial(int(parameters.iloc[k, 3]) - int(parameters.iloc[k, 4]))) / 
                    (math.factorial(int(parameters.iloc[k, 3]) + int(parameters.iloc[k, 4]))))) * polynominal[i] * cos_m_phi[j] * cos_ti[n]

                    r.append(r_i)
                    phi_i.append(phi)
                    theta_i.append(theta)

        r = pd.DataFrame(r)
        phi_i = pd.DataFrame(phi_i)
        theta_i = pd.DataFrame(theta_i)
        
        #coordinates are now spherical and need to be changed into cartesian
        position_vector = pd.concat([r * np.cos(phi_i) * np.sin(theta_i), r * np.sin(phi_i) * np.sin(theta_i), r * np.cos(theta_i)], axis = 1)
        position_vector.columns = ['x', 'y', 'z']

        fig = px.scatter_3d(position_vector, x = 'x' , y = 'y', z = 'z', color =  r,
                            range_color = (1 - 0.7 * sum(amplitude), 1 + 0.7 * sum(amplitude)))
        fig.update_layout(scene_aspectmode='cube', scene = dict(
                          xaxis = dict(range = [- 1 - sum(amplitude), 1 + sum(amplitude)],),
                          yaxis = dict(range = [- 1 - sum(amplitude), 1 + sum(amplitude)],),
                          zaxis = dict(range = [- 1 - sum(amplitude), 1 + sum(amplitude)])))

        fig.show()
        r = []
        phi_i = []
        theta_i = []
        position_vector = pd.DataFrame()
x = position(f_in, t_total, delta_t, n_theta, n_phi)
