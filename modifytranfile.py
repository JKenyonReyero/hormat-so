import re

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


filename = "tran.testKDtli"
read_pots = 0
pots_dict = {"pot_1":None,
             "pot_2":None,
             "pot_3":None,
             "pot_4":None,
             "pot_5":None,
             "pot_6":None,
             }


with open(filename, "r") as file:
    tran_string = file.read()


no_header = tran_string[:]

# Counting the number of (5e14.7) in the tran file to double check it lines up with the previous number of read in files
count_form = tran_string.count("(5e14.7)")
print(f"Number of potentials read in: {count_form}")
if count_form == 0:
    print("""-------------------
          I haven't handled this case yet (no potentials read in).
          Please either add functionality or use different tran file.
          -------------------""")
    
if count_form >= 1 and count_form < 7:
    sections = split_by_form(no_header) # section[0] is header, then some potentials then the potential with footer
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



with open(filename, "r") as file:
    for i, line in enumerate(file):
        print(line)
        if i >= 14:
            break
        fields = line.split()
        # look for entrance channels
        if fields[0] == "5.1000":
            enVD = fields[1]     #! V   ; Real well depth of the optical potential.
            enWD = fields[2]     #! W   ; Imaginary well depth.
            enVSOD = fields[3]   # Vso ; Real well depth of spin-orbit term.
            enWSOD = fields[4]   # Wso ; Imaginary well depth of spin-orbit term.
            enRRD = fields[5]    # r0  ; Real well radius parameter.
            enARD = fields[6]    # ar  ;Real well diffuseness parameter.
            enRCD = fields[7]    # rc  ; Coulomb charge radius parameter.
        elif fields[0] == "6.1000":
            enRSORD = fields[1]  # rsr ; Real well radius parameter of spin-orbit term.
            enASORD = fields[2]  # asr ; Real well diffuseness parameter of spin-orbit term.
            enRSOID = fields[3]  # rsi ; Imaginary well radius parameter of spin-orbit term.
            enASOID = fields[4]  # asi ; Imaginary well diffuseness parameter of spin-orbit term.
            
        elif fields[0] == "7.1000":
            enCSDGD = fields[1]  # Csd ; Mixing factor of volume and surface imaginary well.
            enRID = fields[2]    #! ri  ; Imaginary well radius parameter.
            enAID = fields[3]    # ai  ; Imaginary well diffuseness parameter.
            enRGD = fields[4]    # rg  ; Gaussian type imaginary well radius parameter.
            enAGD = fields[5]    # ag  ;Gaussian type imaginary well range parameter.

        # look for exit channels
        if fields[0] == "5.2000":
            exVD = fields[1]     # V   ; Real well depth of the optical potential.
            exWD = fields[2]     # W   ; Imaginary well depth.
            exVSOD = fields[3]   # Vso ; Real well depth of spin-orbit term.
            exWSOD = fields[4]   # Wso ; Imaginary well depth of spin-orbit term.
            exRRD = fields[5]    # r0  ; Real well radius parameter.
            exARD = fields[6]    # ar  ;Real well diffuseness parameter.
            exRCD = fields[7]    # rc  ; Coulomb charge radius parameter.
        elif fields[0] == "6.2000":
            exRSORD = fields[1]  # rsr ; Real well radius parameter of spin-orbit term.
            exASORD = fields[2]  # asr ; Real well diffuseness parameter of spin-orbit term.
            exRSOID = fields[3]  # rsi ; Imaginary well radius parameter of spin-orbit term.
            exASOID = fields[4]  # asi ; Imaginary well diffuseness parameter of spin-orbit term.
            
        elif fields[0] == "7.2000":
            exCSDGD = fields[1]  # Csd ; Mixing factor of volume and surface imaginary well.
            exRID = fields[2]    #! ri  ; Imaginary well radius parameter.
            exAID = fields[3]    # ai  ; Imaginary well diffuseness parameter.
            exRGD = fields[4]    # rg  ; Gaussian type imaginary well radius parameter.
            exAGD = fields[5]    # ag  ;Gaussian type imaginary well range parameter.

# Check entrance potentials
if enVD == "1.0000" and enRRD == "99.0000":
    read_enrepot = True
    read_pots += 1
else:
    read_enrepot = False
if enWD == "1.0000" and enRID == "99.0000":
    read_enimpot = True
    read_pots += 1
else:
    read_enimpot = False
if enVSOD == "1.0000" and enRSORD == "99.0000":
    read_ensopot = True
    read_pots += 1
else:
    read_ensopot = False

# Check exit potentials
if exVD == "1.0000" and exRRD == "99.0000":
    read_exrepot = True
    read_pots += 1
else:
    read_exrepot = False
if exWD == "1.0000" and exRID == "99.0000":
    read_eximpot = True
    read_pots += 1
else:
    read_eximpot = False
if exVSOD == "1.0000" and exRSORD == "99.0000":
    read_exsopot = True
    read_pots += 1
else:
    read_exsopot = False

# Counting the number of (5e14.7) in the tran file to double check it lines up with the previous number of read in files
count_form = tran_string.count("(5e14.7)")

print(f"Total number of potentials read in in input tran file: {read_pots}")
print(f"Actual number of potentials read in: {count_form}")

# def parse_file(filename):
#     with open(filename, "r") as f:
#         lines = f.readlines()

#     header_lines = []
#     footer_lines = []
#     blocks = []

#     n = len(lines)
#     i = 0

#     # -----------------------
#     # 1. FIND FOOTER START
#     # -----------------------
#     footer_start = None
#     for idx, line in enumerate(lines):
#         if line.startswith("   10.0000"):
#             footer_start = idx
#             break

#     if footer_start is None:
#         raise ValueError("Footer not found (missing '10.0000')")

#     # -----------------------
#     # 2. SPLIT HEADER / BODY / FOOTER
#     # -----------------------
#     header_lines = lines[:footer_start]
#     body_lines = lines[len(header_lines):footer_start]
#     footer_lines = lines[footer_start:]

#     # -----------------------
#     # 3. PARSE BLOCKS (if any)
#     # -----------------------
#     i = 0
#     n_body = len(body_lines)

#     while i < n_body:
#         line = body_lines[i]

#         if "(5e14.7)" in line:
#             block_label = line.strip()
#             i += 1

#             block_data = []

#             while i < n_body and "(5e14.7)" not in body_lines[i]:
#                 block_data.append(body_lines[i])
#                 i += 1

#             blocks.append({
#                 "label": block_label,
#                 "data": "".join(block_data).strip()
#             })
#         else:
#             i += 1

#     # -----------------------
#     # 4. SAFE OUTPUT (0–6 blocks)
#     # -----------------------
#     result = {
#         "header_string": "".join(header_lines).strip(),
#         "potential_1": blocks[0] if len(blocks) > 0 else None,
#         "potential_2": blocks[1] if len(blocks) > 1 else None,
#         "potential_3": blocks[2] if len(blocks) > 2 else None,
#         "potential_4": blocks[3] if len(blocks) > 3 else None,
#         "potential_5": blocks[4] if len(blocks) > 4 else None,
#         "potential_6": blocks[5] if len(blocks) > 5 else None,
#         "footer_string": "".join(footer_lines).strip()
#     }

#     return result

# result = parse_file(filename)

# print(result["header_string"])
# print(result["potential_1"])
# print(result["potential_2"])
# print(result["potential_3"])
# print(result["potential_4"])
# print(result["potential_5"])
# print(result["potential_6"])
# print(result["footer_string"])