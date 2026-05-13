import re
from tran_template import TRAN_TEMPLATE, TFOLD_SO_TEMPLATE


def search_tran_file(tran_string, search_string):
    """
    Search string should be literal (has r in front i.e r"Looking for this")
    """
    # find integrated xsec and return
    match = re.search(
        search_string,
        tran_string
    )
    if match:
        return float(match.group(1))
    else:
        print("-----------------------\nString not found in tran file\n-----------------------")
        return 0  

def split_by_form(no_header):
    # print(no_header)
    pattern = r'(?=^\(5e14\.7\)\s+.*$)'

    sections = re.split(pattern, no_header, flags=re.MULTILINE)

    # remove empty entries
    sections = [s for s in sections if s.strip()]
    print(sections)

    return sections

def split_by_form_foot_start(last_pot_foot):
    parts = re.split(r'(?m)\n(?=\s*10\.0000)', last_pot_foot, maxsplit=1)
    return parts


tran_file = "tran.testKDtli"
tfold_re_file = "tfold_relo.dat"
tfold_im_file = "tfold_imlo.dat"
which_channel_tfold = "ex"  # "en" or "ex"
read_pots = 0
pots_dict = {"pot_1":"",
             "pot_2":"",
             "pot_3":"",
             "pot_4":"",
             "pot_5":"",
             "pot_6":"",
             }
headinfo = {}

with open(tran_file, "r") as file:
    tran_string = file.read()



# Counting the number of (5e14.7) in the tran file to double check it lines up with the previous number of read in files
count_form = tran_string.count("(5e14.7)")
print(f"Number of potentials read in: {count_form}")
if count_form == 0:
    print("""-------------------
          I haven't handled this case yet (no potentials read in).
          Please either add functionality or use different tran file.
          -------------------""")
    
if count_form >= 1 and count_form < 7:
    sections = split_by_form(tran_string) # section[0] is header, then some potentials then the potential with footer
    # print(sections[1:])
    header = sections[0]
    last_pot, footer = split_by_form_foot_start(sections[-1])
    print(last_pot)
    print("---------")
    print(footer)
    print(sections[1:-1]+[last_pot])
    pots_ar = sections[1:-1]+[last_pot]
    for i, v in enumerate(pots_ar):
        pots_dict[f"pot_{i+1}"] = v
    
    # print(pots_dict)
    # I now have: header, the potentials in order and the footer.

split_header = header.split("\n")
print(split_header)
for i in range(0,len(split_header)):  # iterate over each line in the split header
    if (5 <= i <= 7):
        fields = split_header[i].split()
        # look for entrance channel potentials
        if fields[0] == "5.1000":
            headinfo["enVD"] = fields[1]     #! V   ; Real well depth of the optical potential.
            headinfo["enWD"] = fields[2]     #! W   ; Imaginary well depth.
            headinfo["enVSOD"] = fields[3]   # Vso ; Real well depth of spin-orbit term.
            headinfo["enWSOD"] = fields[4]   # Wso ; Imaginary well depth of spin-orbit term.
            headinfo["enRRD"] = fields[5]    # r0  ; Real well radius parameter.
            headinfo["enARD"] = fields[6]    # ar  ;Real well diffuseness parameter.
            headinfo["enRCD"] = fields[7]    # rc  ; Coulomb charge radius parameter.
        elif fields[0] == "6.1000":
            headinfo["enRSORD"] = fields[1]  # rsr ; Real well radius parameter of spin-orbit term.
            headinfo["enASORD"] = fields[2]  # asr ; Real well diffuseness parameter of spin-orbit term.
            headinfo["enRSOID"] = fields[3]  # rsi ; Imaginary well radius parameter of spin-orbit term.
            headinfo["enASOID"] = fields[4]  # asi ; Imaginary well diffuseness parameter of spin-orbit term.
        elif fields[0] == "7.1000":
            headinfo["enCSDGD"] = fields[1]  # Csd ; Mixing factor of volume and surface imaginary well.
            headinfo["enRID"] = fields[2]    #! ri  ; Imaginary well radius parameter.
            headinfo["enAID"] = fields[3]    # ai  ; Imaginary well diffuseness parameter.
            headinfo["enRGD"] = fields[4]    # rg  ; Gaussian type imaginary well radius parameter.
            headinfo["enAGD"] = fields[5]    # ag  ;Gaussian type imaginary well range parameter.
    # elif i == 9:
    #     # read entrance channel proj and targ mass/charge
    #     headinfo["enPMAS"] = fields[1]
    #     headinfo["enTMAS"] = fields[2]
    #     headinfo["enPZ"] = fields[3]
    #     headinfo["enTZ"] = fields[4]
    #     headinfo["enPSPN"] = fields[5]
    #     headinfo["enTSPN"] = fields[6]
    #     headinfo["enQVLUE"] = fields[7]
    elif (10 <= i <= 12):
        # look for exit channel potentials
        fields = split_header[i].split()
        if fields[0] == "5.2000":
            headinfo["exVD"] = fields[1]     #! V   ; Real well depth of the optical potential.
            headinfo["exWD"] = fields[2]     #! W   ; Imaginary well depth.
            headinfo["exVSOD"] = fields[3]   # Vso ; Real well depth of spin-orbit term.
            headinfo["exWSOD"] = fields[4]   # Wso ; Imaginary well depth of spin-orbit term.
            headinfo["exRRD"] = fields[5]    # r0  ; Real well radius parameter.
            headinfo["exARD"] = fields[6]    # ar  ;Real well diffuseness parameter.
            headinfo["exRCD"] = fields[7]    # rc  ; Coulomb charge radius parameter.
        elif fields[0] == "6.2000":
            headinfo["exRSORD"] = fields[1]  # rsr ; Real well radius parameter of spin-orbit term.
            headinfo["exASORD"] = fields[2]  # asr ; Real well diffuseness parameter of spin-orbit term.
            headinfo["exRSOID"] = fields[3]  # rsi ; Imaginary well radius parameter of spin-orbit term.
            headinfo["exASOID"] = fields[4]  # asi ; Imaginary well diffuseness parameter of spin-orbit term.
        elif fields[0] == "7.2000":
            headinfo["exCSDGD"] = fields[1]  # Csd ; Mixing factor of volume and surface imaginary well.
            headinfo["exRID"] = fields[2]    #! ri  ; Imaginary well radius parameter.
            headinfo["exAID"] = fields[3]    # ai  ; Imaginary well diffuseness parameter.
            headinfo["exRGD"] = fields[4]    # rg  ; Gaussian type imaginary well radius parameter.
            headinfo["exAGD"] = fields[5]    # ag  ;Gaussian type imaginary well range parameter.
    else:
        headinfo[f"line_{i+1}"] = split_header[i]

# Check entrance potentials
print("Hi")
bool_read_pot = {
    "check_enrepot": headinfo["enVD"] == "1.0000" and headinfo["enRRD"] == "99.0000",
    "check_enimpot": headinfo["enWD"] == "1.0000" and headinfo["enRID"] == "99.0000",
    "check_ensopot": headinfo["enVSOD"] == "1.0000" and headinfo["enRSORD"] == "99.0000",
    "check_exrepot": headinfo["exVD"] == "1.0000" and headinfo["exRRD"] == "99.0000",
    "check_eximpot": headinfo["exWD"] == "1.0000" and headinfo["exRID"] == "99.0000",
    "check_exsopot": headinfo["exVSOD"] == "1.0000" and headinfo["exRSORD"] == "99.0000",
}

numread_pots = sum(bool_read_pot.values())


# Counting the number of (5e14.7) in the tran file to double check it lines up with the previous number of read in files
count_form = tran_string.count("(5e14.7)")

print(f"Total number of potentials to be read in according to header: {numread_pots}")
print(f"Actual number of potentials read in: {count_form}")
if numread_pots != count_form:
    raise ValueError(f"Number of potentials to be read in according to header ({numread_pots}) does not equal the number of detected potentials in the tran file ({count_form}).")

# Assign pot_x to whichever enrepot ... exsopot is read in
pot_iter = 1
for key, check in bool_read_pot.items():
    if check:
        print(f"First valid pot: {key[6:]} is assigned to pot_{pot_iter}")
        pots_dict[f"{key[6:]}"] = pots_dict[f"pot_{pot_iter}"]
        pot_iter += 1
    else:
        pots_dict[f"{key[6:]}"] = ""
# pots_dict now has 12 keys, 6 that are pot_x and 6 that match the string to which the read in
# print(pots_dict)

#* Read in tfold pot
if not (which_channel_tfold == "en" or which_channel_tfold == "ex"):
    raise ValueError("Don't know what channel to put tfold pot into. which_channel_tfold not clear.")

with open(tfold_re_file, "r") as file:
    tfold_re = file.read()

with open(tfold_im_file, "r") as file:
    tfold_im = file.read()

tfold_so = TFOLD_SO_TEMPLATE.format(channel=which_channel_tfold)  # This is just a big array of zeroes at the moment

# replace card 5.? pot params
headinfo[f"{which_channel_tfold}VD"] = "1.0000"
headinfo[f"{which_channel_tfold}WD"] = "1.0000"
headinfo[f"{which_channel_tfold}VSOD"] = "1.0000"
headinfo[f"{which_channel_tfold}WSOD"] = "0.0000"
headinfo[f"{which_channel_tfold}RRD"] = "99.0000"
headinfo[f"{which_channel_tfold}ARD"] = "1.0000"
# replace card 6.? pot params
headinfo[f"{which_channel_tfold}RSORD"] = "99.0000"
headinfo[f"{which_channel_tfold}ASORD"] = "1.0000"
headinfo[f"{which_channel_tfold}RSOID"] = "1.0000"
headinfo[f"{which_channel_tfold}ASOID"] = "1.0000"
# replace card 7.? pot params
headinfo[f"{which_channel_tfold}CSDGD"] = "0.0000"
headinfo[f"{which_channel_tfold}RID"] = "99.0000"
headinfo[f"{which_channel_tfold}AID"] = "1.0000"
headinfo[f"{which_channel_tfold}RGD"] = "0.0000"
headinfo[f"{which_channel_tfold}AGD"] = "0.0000"

pots_dict[f"{which_channel_tfold}repot"] = tfold_re
pots_dict[f"{which_channel_tfold}impot"] = tfold_im
pots_dict[f"{which_channel_tfold}sopot"] = tfold_so

for key, value in headinfo.items():
    if key[0:2] == "en" or key[0:2] == "ex":
        headinfo[key] = f"{value:>10}"

mod_tran = TRAN_TEMPLATE.format(
    head_line_1=headinfo["line_1"],
    head_line_2=headinfo["line_2"],
    head_line_3=headinfo["line_3"],
    head_line_4=headinfo["line_4"],
    head_line_5=headinfo["line_5"],
    enVD=headinfo["enVD"],
    enWD=headinfo["enWD"],
    enVSOD=headinfo["enVSOD"],
    enWSOD=headinfo["enWSOD"],
    enRRD=headinfo["enRRD"],
    enARD=headinfo["enARD"],
    enRCD=headinfo["enRCD"],
    enRSORD=headinfo["enRSORD"],
    enASORD=headinfo["enASORD"],
    enRSOID=headinfo["enRSOID"],
    enASOID=headinfo["enASOID"],
    enCSDGD=headinfo["enCSDGD"],
    enRID=headinfo["enRID"],
    enAID=headinfo["enAID"],
    enRGD=headinfo["enRGD"],
    enAGD=headinfo["enAGD"],
    head_line_9=headinfo["line_9"],
    head_line_10=headinfo["line_10"],
    exVD=headinfo["exVD"],
    exWD=headinfo["exWD"],
    exVSOD=headinfo["exVSOD"],
    exWSOD=headinfo["exWSOD"],
    exRRD=headinfo["exRRD"],
    exARD=headinfo["exARD"],
    exRCD=headinfo["exRCD"],
    exRSORD=headinfo["exRSORD"],
    exASORD=headinfo["exASORD"],
    exRSOID=headinfo["exRSOID"],
    exASOID=headinfo["exASOID"],
    exCSDGD=headinfo["exCSDGD"],
    exRID=headinfo["exRID"],
    exAID=headinfo["exAID"],
    exRGD=headinfo["exRGD"],
    exAGD=headinfo["exAGD"],
    head_line_14=headinfo["line_14"],
    enrepot=pots_dict["enrepot"],
    enimpot=pots_dict["enimpot"],
    ensopot=pots_dict["ensopot"],
    exrepot=pots_dict["exrepot"],
    eximpot=pots_dict["eximpot"],
    exsopot=pots_dict["exsopot"],
    footer=footer,
                                )
print(mod_tran)

# make datap
with open(f"./tfold_{tran_file}", "w") as f:
    f.write(mod_tran)