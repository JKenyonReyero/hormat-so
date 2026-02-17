import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
from datap_template import DATAP_TEMPLATE


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
    os.makedirs(run_dir+"/hormat/tfnrLELO", exist_ok=True)
    os.makedirs(run_dir+"/hormat/tfnrLENLO", exist_ok=True)
    os.makedirs(run_dir+"/hormat/tfnrNL", exist_ok=True)

    return run_dir, ts

def check_process(result, process_name):
    # Check if the process succeeded
    if result.returncode == 0:
        print(f"{process_name} process finished successfully")
    else:
        print(f"{process_name} process failed with return code {result.returncode}")
        print("Error output:", result.stderr)
        exit()
    pass

def read_wavefunctions(filename):
    data = {}  # dictionary: L -> array([[r, Re, Im], ...])
    current_L = None
    block = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect new block
            if line.startswith("# L="):
                if current_L is not None and block:
                    data[current_L] = np.array(block, dtype=float)
                    block = []
                current_L = int(re.findall(r"\d+", line)[0])

            # Parse data lines
            else:
                # Convert Fortran 'D' exponents to Python 'E'
                line = line.replace("D", "E")
                try:
                    vals = [float(x) for x in line.split()]
                    if len(vals) == 3:
                        block.append(vals)
                except ValueError:
                    continue

        # Save last block
        if current_L is not None and block:
            data[current_L] = np.array(block, dtype=float)

    return data

J = [0,2,4,6,3,0,5,1,3,1,3,2,4,6,8,5,1,7,3,5,3,3,5,5,0,0,2,2,2,4,1,4,0,2,3,2,4,5,2,7,2]
N = [4,3,2,1,2,4,1,3,2,3,2,3,2,1,0,1,3,0,2,1,2,2,1,1,3,3,2,2,2,1,3,1,3,2,2,2,1,1,2,0,2]
E = [0.000,3.661,4.267,4.450,7.010,7.559,7.559,7.710,7.840,7.935,8.281,8.822,9.016,9.102,9.159,11.623,11.635,12.116,12.397,12.419,12.488,12.518,12.605,12.921,14.020,14.984,15.002,15.520,15.755,15.889,15.910,16.121,16.353,16.388,16.872,16.887,17.076,17.111,17.242,17.246,17.343]
#s2n = 15.829  # MeV
#Es2n = E+s2n

#--------------------------------------------------------------
# Take input
#--------------------------------------------------------------
# Heading of files
heading=str(input("Heading file name (max 11 chars)\n"))
use_row=int(input("Do you want to use a row of the table for inputs?\n" \
"0: no, let me enter my own\n" \
"1-41: that row\n"))
if use_row != 0:
    a=92
    z=40
    j=J[use_row-1]
    n=N[use_row-1]
    E_E=E[use_row-1]
    s2n=15.829
    print(f"A={a}, Z={z}, J={j}, n={n}, E_E={E_E}, s2n={s2n}")
else:
    a=int(input("Mass of target\n"))
    z=int(input("Charge of target\n"))
    j=int(input("J (L=J)\n"))
    n=int(input("nodes\n"))
    E_E=float(input("Excitation energy\n"))
    if a==90 and z==40:
        s2n=21.28
        print(s2n)
    elif a==92 and z==40:
        s2n=15.829
        print("s2n = "+str(s2n))
    else:
        s2n=float(input("2 neutron separation energy (MeV)\n"))

elab=float(input("Lab energy (MeV)\n"))


#--------------------------------------------------------------
# Set up directories
#--------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

cwd = os.getcwd()  # directory of wherever the user ran the command
# create directory to write output
run_dir, ts = make_run_folder(heading=heading)
print(f"Creating and writting calculations to:\n{run_dir}\n")


# hormat input:
hormat_in = (
    f"{a}  {z} 1 1\n"  # Targ A, Z, proj A, Z
    f"{elab}\n"      # Lab energy of proj
    "-70.95 1.29 0.58 1.34 0.88\n"
    "-9.03 1.24 0.50  -15.74  1.20 0.45\n"
    "-0.00001   1.20 0.45 0.59\n"
    "00 30\n"            # number of partial waves, Lmin, Lmax
    "  12\n"
    "  3.0  2.71472\n"
    "     16\n"
)
# make file
with open(run_dir+"/hormat/data_hormat", "w") as f:
    f.write(hormat_in)

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
    nstep=200
)
# make datap
with open(run_dir+"/pChan/datap", "w") as f:
    f.write(pChan_dp)

print("Done making inputs")

names = ["hormat", "pChan"]

for i in range(0,len(names)):
    # make front21 input 
    f21_in = (
        f"{names[i]}\n"           # set identifier
        f"{names[i]}\n"           # title information
        "8\n1\n"              # which reaction, read p wfns
        f"0\n{elab}\n"         # calc t wfns, Elab
        f"{a} {z}\n"            # A, Z targ
        f"2\n19.9 0.1\n"      # Integration ranges 0 - 19.9fm in 0.1 fm step
        "2\n30\n"             # number of partial waves
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

    if i == 0:  # hormat
        # Run hormat
        hormat_result = subprocess.run(
            [script_dir+"/hormat-so"], 
            input=hormat_in, 
            text=True, 
            capture_output=True,
            cwd=run_dir+"/"+names[i]
        )
        print("Done running hormat")
        for p in range(0,3):
            # make front for hormat
            end = ["LELO", "LENLO", "NL"]
            files = ["fort.211", "fort.212", "fort.213"]
            # Copy hormat LELO/LENLO/NL wfns to tfnr dirs and rename to fort.16
            shutil.copy(run_dir+"/"+names[i]+"/"+files[p], run_dir+"/"+names[i]+"/tfnr"+end[p]+"/fort.16")
            # Run front
            f21_result = subprocess.run(
                [script_dir+"/front21.out"], 
                input=f21_in, 
                text=True, 
                capture_output=True,
                cwd=run_dir+"/"+names[i]+"/tfnr"+end[p]
            )
            # run twofnr
            twofnr20_result = subprocess.run(
                [script_dir+"/twofnr20.out"], 
                input=f"tran.{names[i]}", 
                text=True, 
                capture_output=True,
                cwd=run_dir+"/"+names[i]+"/tfnr"+end[p]
            )
            print("Done running twoFNR for hormat "+end[p])            
    elif i == 1:  # pChan
        # todo run pChan
        pChan_result = subprocess.run(
            [script_dir+"/pChan_DOM"], 
            input=pChan_in, 
            text=True, 
            capture_output=True,
            cwd=run_dir+"/"+names[i]
        )
        print("Done running pChan")
        shutil.copy(run_dir+"/"+names[i]+"/wfns_proton."+heading, run_dir+"/"+names[i]+"/tfnr/fort.16")
        # Run front
        f21_result = subprocess.run(
            [script_dir+"/front21.out"], 
            input=f21_in, 
            text=True, 
            capture_output=True,
            cwd=run_dir+"/"+names[i]+"/tfnr"
        )
        # run twofnr
        twofnr20_result = subprocess.run(
            [script_dir+"/twofnr20.out"], 
            input=f"tran.{names[i]}", 
            text=True, 
            capture_output=True,
            cwd=run_dir+"/"+names[i]+"/tfnr"
        )
        print("Done running twoFNR for pChan")
# todo plot result of codes
subprocess.run(
    ["xmgrace"] + [run_dir+"/hormat/tfnrNL/21.hormat"] + [run_dir+"/pChan/tfnr/21.pChan"]
)

subprocess.run(
    ["xmgrace"] + [run_dir+"/hormat/tfnrNL/21.hormat"] + [run_dir+"/hormat/tfnrLENLO/21.hormat"] + [run_dir+"/hormat/tfnrLELO/21.hormat"]
)