import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
from datap_template import DATAP_TEMPLATE
from scipy.optimize import minimize_scalar



def make_run_folder(heading: str, base_dir: str = None) -> str:
    # Makes the folder to put all the run info
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"{heading}_{ts}"
    if base_dir is None:
        base_dir = os.getcwd()
    run_dir = os.path.join(base_dir, folder_name)
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(run_dir+"/pChan", exist_ok=True)
    os.makedirs(run_dir+"/pChan/tfnr", exist_ok=True)
    os.makedirs(run_dir+"/hormat", exist_ok=True)
    os.makedirs(run_dir+"/hormat/tfnrNL", exist_ok=True)

    return run_dir, ts


print("This code tries to find the best harmonic oscilator radius for the hormat code.")
print("It compares hormat against pChan (Greg's code) and tries to optimise the (p,t) xsec to be identical.\n")
heading=str(input("Heading file name (max 20 chars)\n"))
a=92
z=40
j=0
n=0
E_E=0
s2n=15.829
print(f"A={a}, Z={z}, J={j}, n={n}, E_E={E_E}, s2n={s2n}")

elab=float(input("Lab energy (MeV)\n"))
lmax = 30

#--------------------------------------------------------------
# Set up directories
#--------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

cwd = os.getcwd()  # directory of wherever the user ran the command
# create directory to write output
run_dir, ts = make_run_folder(heading=heading)
print(f"Creating and writting calculations to:\n{run_dir}\n")

# pChan input:
pChan_in = (
    f"1\n1\n{heading}\n"  # Use non-local, use TPM pot, 
)
# pChan datap input:
pChan_dp = DATAP_TEMPLATE.format(
    a=a,
    z=z,
    elab=elab,
    qval=0.0,
    nstep=200,
    lmax=lmax
)
# make datap
with open(run_dir+"/pChan/datap", "w") as f:
    f.write(pChan_dp)

pChan_result = subprocess.run(
    [script_dir+"/pChan_DOM"], 
    input=pChan_in, 
    text=True, 
    capture_output=True,
    cwd=run_dir+"/pChan"
)

# make front21 input 
f21_in = (
    f"{heading}\n"           # set identifier
    f"{heading}\n"           # title information
    "8\n1\n"              # which reaction, read p wfns
    f"0\n{elab}\n"         # calc t wfns, Elab
    f"{a} {z}\n"            # A, Z targ
    f"2\n19.9 0.1\n"      # Integration ranges 0 - 19.9fm in 0.1 fm step
    f"2\n{lmax}\n"        # number of partial waves
    "0 0 0\n"             # default angles
    f"{j} {j}\n{n}\n"     # quantum numbers L and J and nodes
    f"1\n{E_E+s2n}\n"     # use s2n, excitation E +s2n
    f"1\n"                # no nonlocality in p channel
    "0\n1\n1\n"           # spin in incident, p channel pot (ignore)
    "1\n"                 # no nonlocality in t channel
    f"{j}\n1\n2\n"        # spin in outgoing channel, use Li potential
    "1\n1\n"              # D0, zero-range
    "1.25 0.65\n"         # di-neutron binding potential (radius and diffuseness)
    "0\n0\n"              # di-neutron spin orbit, di-neutron nonlocality
)

# Run front
f21_result = subprocess.run(
    [script_dir+"/front21.out"], 
    input=f21_in, 
    text=True, 
    capture_output=True,
    cwd=run_dir+"/pChan"+"/tfnr"
)

shutil.copy(run_dir+"/pChan/wfns_proton."+heading, run_dir+"/pChan/tfnr/fort.16")
twofnr20_result = subprocess.run(
    [script_dir+"/twofnr20.out"], 
    input=f"tran.{heading}", 
    text=True, 
    capture_output=True,
    cwd=run_dir+"/pChan/tfnr"
)

data = np.loadtxt(run_dir+f"/pChan/tfnr/21.{heading}")
pChan_xsec = data[:,1]
angles = data[:,0]
# Run front
f21_result = subprocess.run(
    [script_dir+"/front21.out"], 
    input=f21_in, 
    text=True, 
    capture_output=True,
    cwd=run_dir+"/hormat"+"/tfnrNL"
)

# make function that returns the 
def func(ho_rad):
    hormat_xsec = run_hormat_tfnr(ho_rad)
    return np.sum((np.log(pChan_xsec) - np.log(hormat_xsec))**2)

def run_hormat_tfnr(ho_rad):
    # hormat input:
    hormat_in = (
        f"{a}  {z} 1 1\n"  # Targ A, Z, proj A, Z
        f"{elab}\n"      # Lab energy of proj
        "-70.95 1.29 0.58 1.34 0.88\n"
        "-9.03 1.24 0.50  -15.74  1.20 0.45\n"
        "-0.00001   1.20 0.45 0.59\n"
        f"00 {lmax}\n"            # number of partial waves, Lmin, Lmax, 00 30
        "  13\n"             # 13
        f"  {ho_rad}\n"   # harmonic oscilator radius 3.0 2.71472------2.68 for 65 MeV
        "     16\n"          # nmax 16
    )
    # make file
    with open(run_dir+"/hormat/data_hormat", "w") as f:
        f.write(hormat_in)
    # Run hormat
    hormat_result = subprocess.run(
        [script_dir+"/hormat-so"], 
        input=hormat_in, 
        text=True, 
        capture_output=True,
        cwd=run_dir+"/hormat"
    )
    shutil.copy(run_dir+"/hormat/fort.213", run_dir+"/hormat/tfnrNL/fort.16")
    # run twofnr
    twofnr20_result = subprocess.run(
        [script_dir+"/twofnr20.out"], 
        input=f"tran.{heading}", 
        text=True, 
        capture_output=True,
        cwd=run_dir+"/hormat/tfnrNL"
    )
    data = np.loadtxt(run_dir+f"/hormat/tfnrNL/21.{heading}")
    hormat_xsec = data[:,1]
    return hormat_xsec

ho_result = minimize_scalar(func, bounds=(2.0, 3.0), method='bounded')
print(f"Optimised harmonic oscilator result is: {ho_result}")