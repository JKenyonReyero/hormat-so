import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
import json


def make_run_folder(loc_types, nodes, ang_moms, heading: str, base_dir: str = None) -> str:
    pmin, pmax = loc_types
    nmin, nmax = nodes
    lmin, lmax = ang_moms
    # Makes the folder to put all the run info
    start_time = datetime.now()
    ts = start_time.strftime("%Y%m%d_%H%M%S")
    folder_name = f"{heading}_{ts}"
    if base_dir is None:
        base_dir = os.getcwd()
    run_dir = os.path.join(base_dir, folder_name)
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(run_dir+"/hormat", exist_ok=True)
    # make all NL/LE dirs
    type_loc = ["NL", "LENLO", "LELO"]
    for i in range(pmin, pmax+1):
        # make all NL/LE dirs
        os.makedirs(run_dir+f"/{type_loc[i]}", exist_ok=True)
        # make all nodes/ang mom dirs
        for p in range(nmin, nmax+1):
            os.makedirs(run_dir+f"/{type_loc[i]}/nodes_{p}", exist_ok=True)  # make nodes dirs
            for j in range(lmin, lmax+1):
                os.makedirs(run_dir+f"/{type_loc[i]}/nodes_{p}/l_{j}", exist_ok=True)  # make ang mom dirs

    

    return run_dir, ts, start_time

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

def save_int_xsec_n_l(my_dict, my_dir):
    np.savez(
        my_dir + "/integrated_xsecs_n_l.npz",
        **{f"{p}_{n}_{l}": arr
        for (p, n, l), arr in my_dict.items()}
    )
    return

def read_int_xsec_n_l(file_path):
    # Currently unused but saving for later if needed
    data = np.load(file_path)

    xsec_data = {
        tuple(k.split("_")[0:3]): data[k]
        for k in data.files
    }
    return xsec_data

a=92
z=40
s2n=15.829
heading=str(input("Heading file name (max 11 chars)\n"))
elab=float(input("Lab energy (MeV)\n"))
nl_pot = int(input("Which Non-local potential?\n[1] Tian Pang Ma potential\n[2] Perey-Buck potential\n"))
if nl_pot==1:  # Tian Pang Ma potential
    Vr, rr, ar = -70.95, 1.29, 0.58
    Wv, rv, av = -9.03,  1.24, 0.50
    Wd, rd, ad = -15.74, 1.20, 0.45
    Vso, rso, aso = -8.13, 1.02, 0.59
    rc = 1.34
    beta = 0.88
elif nl_pot==2:  # Perey-Buck potential
    Vr, rr, ar = -71, 1.22, 0.65
    Wv, rv, av = 0.0, 0.0, 0.0           # not included in perey-buck pot so zero
    Wd, rd, ad = -15, 1.22, 0.47       
    Vso, rso, aso = -14.35, 1.22, 0.65  
    rc = 1.25                          # Generally 1.25
    beta = 0.85

# Set up loop parameters and arrays
start = 0.5
end = 15.829
step = 0.5

# Create the values with arange
s2n = np.arange(start, end, step)

# Ensure the exact end is included
if s2n[-1] != end:
    s2n = np.append(s2n, end)

int_xsec = np.zeros(len(s2n))
int_xsec_n_l = {}

type_loc = ["LELO", "LENLO", "NL"]

pmin = 0  # which potential type to calculate
pmax = 2  # which potential type to calculate
nmin = 0  # min nodes: 0
nmax = 4  # max nodes: 4
lmin = 0  # min ang mom: 0
lmax = 5  # max ang mom: 5


#--------------------------------------------------------------
# Set up directories
#--------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

cwd = os.getcwd()  # directory of wherever the user ran the command
# create directory to write output
run_dir, ts, start_time = make_run_folder((pmin,pmax), (nmin, nmax), (lmin, lmax), heading=heading)
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
    "  3.0  2.71472\n"
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


# loop NL/LE
for p in range(pmin,pmax+1):                     # loop over NL/LE
    p_dir = run_dir+f"/{type_loc[p]}"
    # loop nodes
    for n in range(nmin, nmax+1):        # loop over nodes
        n_dir = p_dir+f"/nodes_{n}"
        # loop ang mom
        for l in range(lmin, lmax+1):    # loop over ang mom
            current_dir = n_dir+f"/l_{l}"  # find current dir
            # write hormat wfns to pot type/nodes/angmom dir
            shutil.copy(run_dir+f"/hormat/fort.{211+p}",current_dir+"/fort.16")
            int_xsec = np.zeros(len(s2n))  # reset integrated xsec
            # loop over s2n
            for i in range(0,len(s2n)):  # loop over s2n
                # make front21 input 
                identifier = f"s2n_{str(s2n[i]).replace('.', ',')}"
                f21_in = (
                    f"{identifier}\n" # set identifier    names are s2n followed by the s2n with a comma
                    f"{identifier}\n" # title information names are s2n followed by the s2n with a comma
                    "8\n1\n"              # which reaction, read p wfns
                    f"0\n{elab}\n"         # calc t wfns, Elab
                    f"{a} {z}\n"            # A, Z targ
                    f"2\n19.9 0.1\n"      # Integration ranges 0 - 19.9fm in 0.1 fm step
                    "2\n30\n"             # number of partial waves
                    "0 0 0\n"             # default angles
                    f"{l} {l}\n{n}\n"     # quantum numbers L and J and nodes
                    f"1\n{s2n[i]}\n"     # use s2n, excitation E +s2n
                    f"1\n"                # no nonlocality in p channel
                    "0\n1\n1\n"           # spin in incident, p channel pot (ignore)
                    "1\n"                 # no nonlocality in t channel
                    f"{l}\n1\n2\n"        # spin in outgoing channel, use Li potential
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
                    cwd=current_dir
                )
                # run twofnr
                twofnr20_result = subprocess.run(
                    [script_dir+"/twofnr20.out"], 
                    input=f"tran.{identifier}", 
                    text=True, 
                    capture_output=True,
                    cwd=current_dir
                )
                # save twofnr output
                with open(current_dir+f"/out.{identifier}", "w") as f:
                    f.write(twofnr20_result.stdout)
                
                int_xsec[i] = find_int_xsec(twofnr20_result)
            
            with open(current_dir+"/int_xsecs", "w") as f:
                for i in range(0,len(s2n)):
                    f.write(f"{s2n[i]:12.6f} {int_xsec[i]:12.6e}\n")
            int_xsec_n_l[(p,n,l)] = int_xsec
            print(f"Finished running twoFNR for {type_loc[p]} and n={n}, l={l}")

print("Finished running, plotting")
# save int_xsec_n_l
save_int_xsec_n_l(int_xsec_n_l, run_dir)

print(f"Integrated cross sections save to: {run_dir+"/integrated_xsecs_n_l.npz"}")
pnl_dict = {
    "pmin": pmin,
    "pmax": pmax,
    "nmin": nmin,
    "nmax": nmax,
    "lmin": lmin,
    "lmax": lmax
}
with open(run_dir+"/pnl.txt", "w") as f:
    json.dump(pnl_dict, f, indent=4)

print(f"pmin, pmax, nmin, nmax, lmin, lmax saved to {run_dir+"/pnl.txt"} in JSON format")

with open(run_dir+"/s2n.txt", "w") as f:
    json.dump(s2n.tolist(), f)

print(f"s2n array saved to {run_dir+"/s2n.txt"} in JSON format")

end_time = datetime.now()
elapsed = end_time - start_time

print(f"Total runtime: {elapsed}")
