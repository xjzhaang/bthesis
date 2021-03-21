from natsort import natsorted
import copy 

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
        names_dict[substrates] = var1
    return names_dict

#Now we parse the equations
def parse_eq():
    names_dict = parse_fceri_meanings()
    f = open("2m_macrovariables1.txt", "r")
    lines = f.readlines()
    eq_dict = {}
    for l in lines:
        receptor1 = []
        receptor2 = []
        name = l.split(" = ")[0].strip()
        var2 = l.split(" = ")[1].strip().split(" + ")[1].strip()
        var = names_dict[var2]
        elements = var.split(".")
        
        for ele in elements:
            if "Rec(a!1" in ele:
                receptor1.append(ele)
            if "Rec(a!2" in ele:
                receptor2.append(ele)
        for ele in elements:
            if  "tSH2!6" in ele:
                if "g~pY!6" in receptor1[0]:
                    receptor1.append(ele)
                if "g~pY!6" in receptor2[0]:
                    receptor2.append(ele)
            if  "tSH2!5" in ele:
                if "g~pY!5" in receptor1[0]:
                    receptor1.append(ele)
                if "g~pY!5" in receptor2[0]:
                    receptor2.append(ele)
            if  "tSH2!4" in ele:
                if "g~pY!4" in receptor1[0]:
                    receptor1.append(ele)
                if "g~pY!4" in receptor2[0]:
                    receptor2.append(ele)
            if  "tSH2!3" in ele:
                if "g~pY!3" in receptor1[0]:
                    receptor1.append(ele)
                if "g~pY!3" in receptor2[0]:
                    receptor2.append(ele)
        for ele in elements:
            if "Lyn(SH2!3" in ele:
                if "b~pY!3" in receptor1[0]:
                    receptor1.append(ele)
                if "b~pY!3" in receptor2[0]:
                    receptor2.append(ele)
            if "Lyn(SH2!4" in ele:
                if "b~pY!4" in receptor1[0]:
                    receptor1.append(ele)
                if "b~pY!4" in receptor2[0]:
                    receptor2.append(ele)
            if "U!3" in ele:
                if "b~Y!3" in receptor1[0]:
                    receptor1.append(ele)
                if "b~Y!3" in receptor2[0]:
                    receptor2.append(ele)
            if "U!4" in ele:
                if "b~Y!4" in receptor1[0]:
                    receptor1.append(ele)
                if "b~Y!4" in receptor2[0]:
                    receptor2.append(ele)
        if "b~Y" in receptor1[0]:
            receptor1[0] = "Rec(b)"
        if "b~pY" in receptor1[0]: 
            receptor1[0] = "Rec(bp)"
        if "b~Y" in receptor2[0]:
            receptor2[0] = "Rec(b)"
        if "b~pY" in receptor2[0]: 
            receptor2[0] = "Rec(bp)"

        if "Syk(a~Y,l~Y" in receptor1[1]:
            receptor1[1] = "Syk(a,l)"
        if "Syk(a~Y,l~Y" in receptor2[1]:
            receptor2[1] = "Syk(a,l)"

        if "Syk(a~pY,l~Y" in receptor1[1]:
            receptor1[1] = "Syk(ap,l)"
        if "Syk(a~pY,l~Y" in receptor2[1]:
            receptor2[1] = "Syk(ap,l)"

        if "Syk(a~Y,l~pY" in receptor1[1]:
            receptor1[1] = "Syk(a,lp)"
        if "Syk(a~Y,l~pY" in receptor2[1]:
            receptor2[1] = "Syk(a,lp)"

        if "Syk(a~pY,l~pY" in receptor1[1]:
            receptor1[1] = "Syk(ap,lp)"
        if "Syk(a~pY,l~pY" in receptor2[1]:
            receptor2[1] = "Syk(ap,lp)"
        
        if len(receptor1) == 3:
            if "Lyn(SH2,U!3)" in receptor1[2]:
                receptor1[2] = "Lyn(c)"
            if "Lyn(SH2,U!4)" in receptor1[2]:
                receptor1[2] = "Lyn(c)"
            if "Lyn(SH2!3,U)" in receptor1[2]:
                receptor1[2] = "Lyn(s)"
            if "Lyn(SH2!4,U)" in receptor1[2]:
                receptor1[2] = "Lyn(s)"
        if len(receptor2) == 3:
            if "Lyn(SH2,U!4)" in receptor2[2]:
                receptor2[2] = "Lyn(c)"
            if "Lyn(SH2,U!3)" in receptor2[2]:
                receptor2[2] = "Lyn(c)"
            if "Lyn(SH2!3,U)" in receptor2[2]:
                receptor2[2] = "Lyn(s)"
            if "Lyn(SH2!4,U)" in receptor2[2]:
                receptor2[2] = "Lyn(s)"
        eq_dict[name] = (receptor1,receptor2)
    return eq_dict

#Check if all 84 are unique
def check_diff():
    eq_dict = parse_eq()
    sett = []
    for values in eq_dict.values():
        if values not in sett:
            sett.append(values)
        else:
            print("bruh")

#Now we construct our own tuples
def construct_pairs():
    receptors = ["Rec(b)", "Rec(bp)"]
    Lyn = ["", "Lyn(c)", "Lyn(s)"]
    Syk = ["Syk(a,l)","Syk(a,lp)","Syk(ap,l)","Syk(ap,lp)"]
    lists = []
    for k in range(len(Lyn)):
        for i in range(len(receptors)):
            for j in range(len(Syk)):
                if k == 0:
                    lists.append([receptors[i],Syk[j]])
                elif k == 1 and i == 0:
                    lists.append([receptors[i],Syk[j],Lyn[k]])
                elif k == 2 and i == 1:
                    lists.append([receptors[i],Syk[j],Lyn[k]])
    lists2 = []
    for element in lists:
        if element[1] == "Syk(a,lp)" or element[1] == "Syk(ap,lp)":
            lists2.append(element)
    list_of_tuples = []
    for element1 in lists:
        for element2 in lists2:
            if element1 != element2:
                if (element2, element1) not in list_of_tuples:
                    list_of_tuples.append((element1, element2))
    return list_of_tuples

#we compare the outputs
def compare_outputs():
    eq = parse_eq()
    list_of_tuples = construct_pairs()
    list2 = copy.deepcopy(list_of_tuples)
    for values in eq.values():
        for item in list_of_tuples:
            if values == item:
                list2.remove(item)
            if values[0] == item[1] and values[1] == item[0]:
                list2.remove(item)
    print(list2)
compare_outputs()

with open('fceri_meanings','w') as afile:
    names_dict = parse_fceri_meanings()
    list_of_strings = [ f'{key} : {names_dict[key]}' for key in natsorted(names_dict)]
    [afile.write(f'{st}\n') for st in list_of_strings ]
with open('fceri_eq','w') as afile:
    eq = parse_eq()
    list_of_strings = [ f'{key} : {eq[key]}' for key in natsorted(eq)]
    [afile.write(f'{st}\n') for st in list_of_strings ]

