import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
import json


def make_run_folder(heading: str, base_dir: str = None) -> str:
    # Makes the folder to put all the run info
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"{heading}_{ts}"
    if base_dir is None:
        base_dir = os.getcwd()
    run_dir = os.path.join(base_dir, folder_name)
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(run_dir+"/hormat", exist_ok=True)
    os.makedirs(run_dir+"/NL", exist_ok=True)
    os.makedirs(run_dir+"/LENLO", exist_ok=True)
    os.makedirs(run_dir+"/LELO", exist_ok=True)
    

    return run_dir, ts

def find_int_xsec(twofnr20_result):
    # find integrated xsec and return
    match = re.search(
        r"integrated cross section=\s*([0-9.+-Ee]+)",
        twofnr20_result.stdout
    )
    if match:
        return float(match.group(1))
    else:
        print("-----------------------\nNO OUTPUT INTEGRATED XSEC FOUND FROM TWOFNR\n-----------------------")
        return 0 
    
def find_exit_ecm(twofnr20_result):
    match = re.search(
        r"exit channel\s*.*?\becm=\s*([0-9.+-Ee]+)",
        twofnr20_result.stdout,
        re.DOTALL
    )
    if match:
        return float(match.group(1))
    else:
        print("-----------------------\nNO EXIT CHANNEL ECM FOUND\n-----------------------")
        return 0

J = [0,2,4,6,3,0,5,1,3,1,3,2,4,6,8,5,1,7,3,5,3,3,5,5,0,0,2,2,2,4,1,4,0,2,3,2,4,5,2,7,2]
N = [4,3,2,1,2,4,1,3,2,3,2,3,2,1,0,1,3,0,2,1,2,2,1,1,3,3,2,2,2,1,3,1,3,2,2,2,1,1,2,0,2]
E = [0.000,3.661,4.267,4.450,7.010,7.559,7.559,7.710,7.840,7.935,8.281,8.822,9.016,9.102,9.159,11.623,11.635,12.116,12.397,12.419,12.488,12.518,12.605,12.921,14.020,14.984,15.002,15.520,15.755,15.889,15.910,16.121,16.353,16.388,16.872,16.887,17.076,17.111,17.242,17.246,17.343]

a=92
z=40
s2n=15.829
heading=str(input("Heading file name (max 11 chars)\n"))
elab=float(input("Lab energy (MeV)\n"))
nl_pot = int(input("Which Non-local potential?\n[1] Tian Pang Ma potential\n[2] Perey-Buck potential\n"))
li_or_pang = int(input("[1] Li [2] Pang\n"))
if li_or_pang == 1:
    li_or_pang = 2
elif li_or_pang == 2:
    li_or_pang = 3
else:
    print("---Wrong---")
if nl_pot==1:  # Tian Pang Ma potential
    Vr, rr, ar = -70.95, 1.29, 0.58
    Wv, rv, av = -9.03,  1.24, 0.50
    Wd, rd, ad = -15.74, 1.20, 0.45
    Vso, rso, aso = -8.13, 1.02, 0.59  # -8.13, 1.02, 0.59
    rc = 1.34
    beta = 0.88
elif nl_pot==2:  # Perey-Buck potential
    Vr, rr, ar = -71, 1.22, 0.65
    Wv, rv, av = 0.0, 0.0, 0.0           # not included in perey-buck pot so zero
    Wd, rd, ad = -15, 1.22, 0.47       
    Vso, rso, aso = -14.35, 1.22, 0.65 # -14.35, 1.22, 0.65
    rc = 1.25                          # Generally 1.25
    beta = 0.85

nl_cor = int(input("Would you like non-local corrections in the triton channel?\n[1] No\n[2] Yes, default 0.20fm\n[3] Yes, I want to specify\n"))
if nl_cor == 1:
    nl_range = ""
elif nl_cor == 2:
    nl_range = "0.20\n"
elif nl_cor == 3:
    nl_range = str(float(input("What non-locality range do you want?")))+"\n"
    nl_cor = 2
elif nl_cor < 1 or nl_cor > 3:
    print("---Wrong nl_cor---")

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
    f"{Vr} {rr} {ar} {rc} {beta}\n"
    f"{Wv} {rv} {av}  {Wd}  {rd} {ad}\n"
    f"{Vso}   {rso} {aso}\n"
    "00 30\n"            # number of partial waves, Lmin, Lmax
    "  12\n"
    "  2.5\n"
    "     16\n"
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
print("Finished running hormat")
found_exit_ecm = False
exit_ecm = []

loc_type = ["LELO","LENLO","NL"]
for p in range(0,3):
    # write hormat wfns to row dir
    shutil.copy(run_dir+f"/hormat/fort.{211+p}",run_dir+f"/{loc_type[p]}/fort.16")

    int_xsec = []

    for i in range(1,42):
        # make front21 input 
        f21_in = (
            f"row_{i}\n"           # set identifier
            f"row_{i}\n"           # title information
            "8\n1\n"              # which reaction, read p wfns
            f"0\n{elab}\n"         # calc t wfns, Elab
            f"{a} {z}\n"            # A, Z targ
            f"2\n19.9 0.1\n"      # Integration ranges 0 - 19.9fm in 0.1 fm step
            "2\n30\n"             # number of partial waves
            "0 0 0\n"             # default angles
            f"{J[i-1]} {J[i-1]}\n{N[i-1]}\n"     # quantum numbers L and J and nodes
            f"1\n{E[i-1]+s2n}\n"     # use s2n, excitation E +s2n
            f"1\n"                # no nonlocality in p channel
            "0\n1\n1\n"           # spin in incident, p channel pot (ignore)
            f"{nl_cor}\n{nl_range}"                 # no nonlocality in t channel
            f"{J[i-1]}\n1\n{li_or_pang}\n"        # spin in outgoing channel, use 2:Li potential 3:Pang
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
            cwd=run_dir+f"/{loc_type[p]}"
        )
        # run twofnr
        twofnr20_result = subprocess.run(
            [script_dir+"/twofnr20nlmod.out"], 
            input=f"tran.row_{i}", 
            text=True, 
            capture_output=True,
            cwd=run_dir+f"/{loc_type[p]}"
        )
        # save twofnr output
        with open(run_dir+f"/{loc_type[p]}/out.row_{i}", "w") as f:
            f.write(twofnr20_result.stdout)
        # save tot xsec to array

        int_xsec.append(find_int_xsec(twofnr20_result))
        if found_exit_ecm == False:
            exit_ecm.append(find_exit_ecm(twofnr20_result))
    found_exit_ecm = True

    with open(run_dir+f"/{loc_type[p]}/int_xsecs", "w") as f:
        for i in range(1,42):
            f.write(f"{E[i-1]:12.6f} {int_xsec[i-1]:12.6e}\n")
    print(f"Finished running twoFNR for hormat {loc_type[p]}")

with open(run_dir+"/s2n.txt", "w") as f:
    json.dump((np.array(E)+15.829).tolist(), f)

with open(run_dir+"/exit_ecm.txt", "w") as f:
    json.dump(exit_ecm, f)