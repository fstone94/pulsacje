import sys 

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
