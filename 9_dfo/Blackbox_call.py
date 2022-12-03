import numpy as np
import os
import subprocess
import json
from subprocess import PIPE,STDOUT
from SML.src.SML import Kriging

def log_output(file, x, g, f):
    x_str = " ".join(map(str, x))
    g_str = " ".join(map(str, g))
    with open(file, "a") as fid:
        fid.write("%s %s %f\n" %(x_str,g_str,f))

def read_output(file):
    with open(file, "r") as f:
        outs = [[float(num) for num in line.rstrip().split()] for line in f]
    return outs

def try_remove(filename):
    # Delete files
    if os.path.isfile(filename):
        os.remove(filename)

def system_command(command):
    #launches a string command in the windows terminal
    p = subprocess.Popen(command,shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT,
                            ) # disable windows errors

    for line in iter(p.stdout.readline, b''):
        line = line.decode('utf-8')
        print(line.rstrip()) # print line by line
        # rstrip() to reomove \n separator

    p.communicate()[0]
    rc = p.returncode
    # print('returncode : %i' %rc)
    if rc != 0:
        raise Exception('The following command failed : %s' %(command))
    return rc

def Blackbox_call(x, *argv):

    # # load parameters
    # with open("bb_parameters.json","r") as file:
    #     bb_params = json.load(file)

    bb_params = argv[0] # retrieve the blackbox parameters vector p, and other blackbox options

    # %% Run Blackbox

    # Delete output files (from previous runs)
    try_remove("output.txt")

    command_str = "-x "+" ".join(map(str, x))+" -p "+" ".join(map(str, bb_params["parameters"]))
    if bb_params["use_gradients"]:
        command_str += " -G -v"
    command = "./bb "+command_str # <--------------------------------for MacOS/Linux
    # command = "bb "+command_str # <--------------------------------for windows

    if not bb_params["use_surrogate"]:
        #####################
        # Real model
        #####################
        status = system_command(command) # this runs the blackbox
        outs = read_output("output.txt")
        f = outs[0][0]
        g = outs[1]
        if bb_params["use_gradients"]:
            df = outs[2]
            dg = np.array([outs[3], outs[4], outs[5]])
        else:
            df = []
            dg = []

        #####################
    else:
        #####################
        # Surrogate model
        #####################
        sm = bb_params["sur_model"] # use this to load a trained model
        sm.x.points = np.array(x).reshape(-1,2)
        sm.predict()
        f = float(sm.yp.points[0,0])
        g = list(sm.yp.points[0,1:])
        pass
        #####################

    if bb_params["save_data"]:
        log_output("log.txt",x,g,f)

    if bb_params["use_gradients"]:
        return [f,g,df,dg]
    else:
        return [f,g]


if __name__ == "__main__":
    save_data = True
    sur = False
    use_gradients = True
    p = [8, -4, -3, -3, -3, 3, 2, 0.8]
    sur_model = None

    P = [save_data,sur,use_gradients,p,sur_model]
    x = np.array([1.0,2.0])

    with open('log.txt', 'w'):
        os.utime('log.txt', None)

    os.chdir('./9_dfo')
    f,df,g,dg = Blackbox_call(x,P)

    print(f)
    print(df)
    print(g)
    print(dg)