from natsort import natsorted
import copy
import re

#We first parse the meanings of bngl into a dictionary
def parse_fceri_meanings():
    f = open("fceri_jiNames.txt", "r")
    lines = f.readlines()
    names_dict = {}
    for l in lines:
        substrates = l.split(": ")[0].strip()
        substrates = substrates.replace("s","S")
        var1 = l.split(": ")[1].strip()
        var1 = var1.replace("Lig(l!1,l!2).", "")
        var1 = var1.replace("Lig(l!1,l).", "Lig(1)")
        names_dict[substrates] = var1
    return names_dict

def parse_1():
    names_dict = parse_fceri_meanings()
    f = open("3m_1", "r")
    lines = f.readlines()
    eq_dict = {}
    for l in lines:
        name = l.split(" = ")[0].strip()
        var2 = l.split(" = ")[1].strip()
        var2 = var2.split(" + ")
        elements = []
        for var in var2:
            
            elements.append(names_dict[var])
        eq_dict[name] = elements
    return eq_dict


def parse_2(txt):
    names_dict = parse_fceri_meanings()
    f = open(txt, "r")
    lines = f.readlines()
    eq_dict = {}
    for l in lines:
        name = l.split(" = ")[0].strip()
        var2 = l.split(" = ")[1].strip()
        var2 = var2.split(" + ")
        list_ = []

        for var_ in var2:
            times2 = False
            receptor1 = []
            receptor2 = []
            if "2S" in var_:
                var_ = var_[1:]
                times2 = True
            var = names_dict[var_]
            elements = var.split(".")
            # separateing the receptors
            non_receptors = []
            for ele in elements:
                if "Rec(a!1" in ele:
                    receptor1.append(ele)
                elif "Rec(a!2" in ele:
                    receptor2.append(ele)
                else:
                    non_receptors.append(ele)
            elements = non_receptors
            # splitting the stuff between the receptors
            connection_re = r'!\d'
            for ele in elements:
                c = re.search(connection_re, ele)
                if c is not None:
                    if receptor1[0].find(c[0]) != -1:
                        receptor1.append(ele)
                    elif receptor2[0].find(c[0]) != -1:
                        receptor2.append(ele)
                    else:
                        print(f"No owner found for {ele}")
                else:
                    print(f"Unconnected element {ele}")
            #print(receptor1)
            # sorting elements: making Syk going before Lyn
            def sort_rec(receptor):
                if len(receptor) > 2 and receptor[1] < receptor[2]:    
                    receptor[1], receptor[2] = receptor[2], receptor[1]
            sort_rec(receptor1)
            print(receptor1)
            sort_rec(receptor2)

            # translation/simplification
            def subs_rec(receptor):
                if "b~Y" in receptor[0]:
                    if "g~Y" in receptor[0]:
                        receptor[0] = "Rec(b,g)"
                if "b~pY" in receptor[0]:
                    if "g~Y" in receptor[0]:
                        receptor[0] = "Rec(bp,g)"
                if "b~Y" in receptor[0]:
                    if "g~pY" in receptor[0]:
                        receptor[0] = "Rec(b,gp)"
                if "b~pY" in receptor[0]:
                    if "g~pY" in receptor[0]:
                        receptor[0] = "Rec(bp,gp)"
            subs_rec(receptor1)
            subs_rec(receptor2)

            def subs_syk(receptor):
                if len(receptor) >= 2:
                    if "Syk(a~Y,l~Y" in receptor[1]:
                        receptor[1] = "Syk(a,l)"
                    if "Syk(a~pY,l~Y" in receptor[1]:
                        receptor[1] = "Syk(ap,l)"
                    if "Syk(a~Y,l~pY" in receptor[1]:
                        receptor[1] = "Syk(a,lp)"
                    if "Syk(a~pY,l~pY" in receptor[1]:
                        receptor[1] = "Syk(ap,lp)"
            subs_syk(receptor1)
            subs_syk(receptor2)

            def subs_lyn(receptor):
                if len(receptor) == 3:
                    if re.search(r'Lyn\(SH2,U!\d\)', receptor[2]) is not None:
                        receptor[2] = "Lyn(c)"
                    elif re.search(r'Lyn\(SH2!\d,U\)', receptor[2]) is not None:
                        receptor[2] = "Lyn(s)"
                    else:
                        print(f"Unrecognized Lyn: {receptor[2]}")
                if len(receptor) == 2:
                    if re.search(r'Lyn\(SH2,U!\d\)', receptor[1]) is not None:
                        receptor[1] = "Lyn(c)"
                    elif re.search(r'Lyn\(SH2!\d,U\)', receptor[1]) is not None:
                        receptor[1] = "Lyn(s)"
            subs_lyn(receptor1)
            subs_lyn(receptor2)
            
            if times2:
                receptor1.append("2x")
            list_.append((receptor1,receptor2))
        eq_dict[name] = list_
    return eq_dict

with open('fceri_3_1','w') as afile:
    eq = parse_1()
    list_of_strings = [ f'{key} : {eq[key]}' for key in natsorted(eq)]
    [afile.write(f'{st}\n') for st in list_of_strings ]

with open('fceri_3_2','w') as afile:
    eq = parse_2("3m_2")
    list_of_strings = [ f'{key} : {eq[key]}' for key in natsorted(eq)]
    [afile.write(f'{st}\n') for st in list_of_strings ]


with open('fceri_3_3','w') as afile:
    eq = parse_2("3m_3")
    list_of_strings = [ f'{key} : {eq[key]}' for key in natsorted(eq)]
    [afile.write(f'{st}\n') for st in list_of_strings ]


with open('fceri_3_4','w') as afile:
    eq = parse_2("3m_4")
    list_of_strings = [ f'{key} : {eq[key]}' for key in natsorted(eq)]
    [afile.write(f'{st}\n') for st in list_of_strings ]