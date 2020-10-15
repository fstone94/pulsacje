import sys 
import scipy.special as sci
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def welcome_msg():
    print("This programme generates a *.gif animation showing how star's surface pulsates in time")
    print()
    print("Usage: python3 pulsacje.py <f_in=> <f_out=> <t_total=> <delta_t=> <n_theta=> <n_phi=>")
    print("f_in     --> name of input file containing modes parameters")
    print("f_out    --> name of the output file")
    print("t_total  --> total animation time expressed in days")
    print("delta_t  --> time step in the animation expressed in days")
    print("n_theta  --> number of discrete points on the star's astrophysical longitude")
    print("n_phi    --> number of discrete points on the star's astrophysical latitude")
    print()
    print("Content of the input file:")
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
    num_lines = len(open(f_in).readlines())

    #quotient of t_total and delta_t will be the number of images in animation
    n_img = int(t_total / delta_t)
    
    #theta angle runs from 0 to pi
    #phi angle runs from 0 to 2pi 
    #division of those ranges by number of discrete points n_theta and n_phi - delta_theta/phi 
    delta_theta = np.pi / n_theta
    delta_phi = 2 * np.pi / n_phi

    cos_theta = []
    polynominal = []
    for i in range (1, n_theta):
        cos1 = np.cos((i - 1) * delta_theta)
        cos_theta.append(cos1)
        
        for j in range (0, num_lines-1):
            y = sci.lpmv(parameters[j, 4], parameters[j, 3], cos_theta[i - 1])
            polynominal.append(y)
 
    cos_m_phi = []
    for i in range (1, n_phi):
        phi_i = (i - 1) * delta_phi
        
        for j in range (0, num_lines - 1 ):
            cos2 = np.cos(phi_i * parameters[j, 4])
            cos_m_phi.append(cos2)

    cos_ti = []
    for i in range (0, n_img):
        t_i = (i - 1) * delta_t
        arg = np.cos(2 * np.pi * parameters[0, 0] * t_i)
        cos_ti.append(arg)

    amplitude = []
    for i in range (0, num_lines - 1 ):
        amplitude.append(parameters[0, 1])
    
    phi_i = []
    theta_i = []
    r = []
    for n in range (1, n_img):
        for i in range(1, n_theta):
            for j in range (1, n_phi):

                r_i = 1
               
                for k in range (0, num_lines - 1 ):
                    theta = (i - 1) * delta_theta
                    phi = (j - 1) * delta_phi

                    r_i = 1 + parameters[k, 1] * math.sqrt(((2 * parameters[k, 3] + 1) / (4 * np.pi)) * 
                    (math.factorial(int(parameters[k, 3]) - int(parameters[k, 4])) / 
                    (math.factorial(int(parameters[k, 3]) + int(parameters[k, 4]))))) * polynominal[i] * cos_m_phi[j] * cos_ti[n]
    
                    r.append(r_i)
                    phi_i.append(phi)
                    theta_i.append(theta)


x = position(f_in, t_total, delta_t, n_theta, n_phi)

