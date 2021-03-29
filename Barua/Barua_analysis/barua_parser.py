from natsort import natsorted
import copy
import re

#We first parse the meanings of bngl into a dictionary
def parse_barua_meanings():
    f = open("baruaNames.txt", "r")
    lines = f.readlines()
    names_dict = {}
    for l in lines:
        substrates = l.split(": ")[0].strip()
        substrates = substrates.replace("s","S")
        var1 = l.split(": ")[1].strip()
        #var1 = var1.replace("Lig(l!1,l!2).", "")
        names_dict[substrates] = var1
    return names_dict

#Now we parse the equations
def parse_eq(txt):
    names_dict = parse_barua_meanings()
    f = open(txt, "r")
    lines = f.readlines()
    eq_dict = {}
    for l in lines:
        name = l.split(" = ")[0].strip()
        var2 = l.split(" = ")[1].strip()
        var2 = var2.split(" + ")
        elements = []
        for var in var2:
            if "2S" in var:
                elements.append(names_dict[var[1:]] + " x2 ")
            else:
                elements.append(names_dict[var])
        eq_dict[name] = elements
    return eq_dict


#compare_outputs()

#exit()

with open('barua_meanings','w') as afile:
    names_dict = parse_barua_meanings()
    list_of_strings = [ f'{key} : {names_dict[key]}' for key in natsorted(names_dict)]
    [afile.write(f'{st}\n') for st in list_of_strings ]

with open('2m/barua_eq.txt','w') as afile:
    names_dict = parse_eq("2m/2m_nontrivial.txt")
    list_of_strings = [ f'{key} : {names_dict[key]}' for key in natsorted(names_dict)]
    [afile.write(f'{st}\n') for st in list_of_strings ]
