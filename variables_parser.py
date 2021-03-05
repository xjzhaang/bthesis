import copy
def read_variables():
    #to 1479
    n = 0
    f = open("fceri_ode.txt", "r")
    lines = f.readlines()
    translations = {"S0": "IgE", "S1": "Lyn", "S2": "Syk", "S3": "r"}
    removal_dict = {"Lyn": ["cLyn", "SH2Lyn"], "cgammaP": "SH2gammaP", "SH2gammaP": "cgammaP"}
    for l in lines:
        n += 1
        #print(n)
        substrates = l.split("->")[0].strip()
        if " + " in substrates:
            sub1 = substrates.split(" + ")[0]
            sub2 = substrates.split(" + ")[1]
            two_sub = True   
        else:
            sub1 = substrates
            two_sub = False
        res_sub = l.split("->")[1].split(" ,")[0].strip()
        if " + " in res_sub:
            res_sub1 = res_sub.split(" + ")[0]
            res_sub2 = res_sub.split(" + ")[1]
            two_res_sub = True   
        else:
            res_sub1 = res_sub
            two_res_sub = False
        rate = l.split("->")[1].split(" , ")[1].strip()
        if two_sub:
            if sub1 in translations and sub2 in translations:
                if rate == "2.0 * pkp1":
                    if type(translations[sub1]) != list and type(translations[sub2]) != list:
                        res_meaning = [translations[sub1], translations[sub2]]
                    else:
                        if type(translations[sub1]) != list and type(translations[sub2]) == list:
                            res_meaning = copy.deepcopy(translations[sub2])
                            res_meaning.append(translations[sub1])
                        if type(translations[sub2]) != list and type(translations[sub1]) == list:
                            res_meaning = copy.deepcopy(translations[sub1])
                            res_meaning.append(translations[sub2])
                        if type(translations[sub1]) == list and type(translations[sub2]) == list:
                            res_meaning = copy.deepcopy(translations[sub1])
                            for item in translations[sub2]:
                                res_meaning.append(item)
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "pkp2":
                    if type(translations[sub1]) != list and type(translations[sub2]) != list:
                        res_meaning = [translations[sub1], translations[sub2]]
                    else:
                        if type(translations[sub1]) != list and type(translations[sub2]) == list:
                            res_meaning = copy.deepcopy(translations[sub2])
                            res_meaning.append(translations[sub1])
                        if type(translations[sub2]) != list and type(translations[sub1]) == list:
                            res_meaning = copy.deepcopy(translations[sub1])
                            res_meaning.append(translations[sub2])
                        if type(translations[sub1]) == list and type(translations[sub2]) == list:
                            res_meaning = copy.deepcopy(translations[sub1])
                            for item in translations[sub2]:
                                res_meaning.append(item)
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "pkpL":
                    if translations[sub1] == "Lyn":
                        if type(translations[sub2]) != list:
                            res_meaning = [translations[sub2], "cLyn"]
                        else:
                            res_meaning = copy.deepcopy(translations[sub2])
                            res_meaning.append("cLyn")
                    else:
                        if type(translations[sub1]) != list:
                            res_meaning = [translations[sub1], "cLyn"]
                        else:
                            res_meaning = copy.deepcopy(translations[sub1])
                            res_meaning.append("cLyn")
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "pkpLs":
                    if translations[sub1] == "Lyn":
                        if type(translations[sub2]) != list:
                            res_meaning = [translations[sub2], "SH2Lyn"]
                        else:
                            res_meaning = copy.deepcopy(translations[sub2])
                            res_meaning.append("SH2Lyn")
                    else:
                        if type(translations[sub1]) != list:
                            res_meaning = [translations[sub1], "SH2Lyn"]
                        else:
                            res_meaning = copy.deepcopy(translations[sub1])
                            res_meaning.append("SH2Lyn")
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "2.0 * pkpS":
                    if translations[sub1] == "Syk":
                        if type(translations[sub2]) != list:
                            res_meaning = [translations[sub2], "Syk"]
                        else:
                            res_meaning = copy.deepcopy(translations[sub2])
                            res_meaning.append("Syk")
                    else:
                        if type(translations[sub1]) != list:
                            res_meaning = [translations[sub1], "Syk"]
                        else:
                            res_meaning = copy.deepcopy(translations[sub1])
                            res_meaning.append("Syk")
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "pkpS":
                    if type(translations[sub1]) != list and type(translations[sub2]) != list:
                        res_meaning = [translations[sub1], translations[sub2]]
                    else:
                        if type(translations[sub1]) != list and type(translations[sub2]) == list:
                            res_meaning = copy.deepcopy(translations[sub2])
                            res_meaning.append(translations[sub1])
                        if type(translations[sub2]) != list and type(translations[sub1]) == list:
                            res_meaning = copy.deepcopy(translations[sub1])
                            res_meaning.append(translations[sub2])
                        if type(translations[sub1]) == list and type(translations[sub2]) == list:
                            res_meaning = copy.deepcopy(translations[sub1])
                            for item in translations[sub2]:
                                res_meaning.append(item)
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
        if not two_res_sub:
            if sub1 in translations and sub2 in translations:
                if rate == "ppLb":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("cbetaP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLg":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("cgammaP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLbs":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("SH2betaP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLgs":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("SH2gammaP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLS":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("cSykP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLSs":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("SH2SykP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppSS":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("NoLoopSykP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppSSs":
                    res_meaning = copy.deepcopy(translations[sub1])
                    res_meaning.append("LoopSykP")
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
        if two_res_sub:
            if sub1 in translations:
                if rate == "pkm1":
                    res_meaning_2 = copy.deepcopy(translations[sub1])
                    for item in translations[res_sub1]:
                        if item in res_meaning_2:
                            res_meaning_2.remove(item)
                        elif item in removal_dict:
                            res_meaning_2.remove(removal_dict[item])
                    if res_sub2 not in translations:
                        translations[res_sub2] = res_meaning_2
                if rate == "pkm2":
                    res_meaning_2 = copy.deepcopy(translations[sub1])
                    for item in translations[res_sub1]:
                        if item in res_meaning_2:
                            res_meaning_2.remove(item)
                        elif item in removal_dict:
                            res_meaning_2.remove(removal_dict[item])
                    if res_sub2 not in translations:
                        translations[res_sub2] = res_meaning_2
                if rate == "pkmS":
                    res_meaning_2 = copy.deepcopy(translations[sub1])
                    print(res_meaning_2)
                    print(n)
                    print(translations[res_sub1])
                    if type(translations[res_sub1]) != list:
                        res_meaning_2.remove(translations[res_sub1])
                    else:
                        for item in translations[res_sub1]:
                            if item in res_meaning_2:
                                res_meaning_2.remove(item)
                            elif item in removal_dict:
                                res_meaning_2.remove(removal_dict[item])
                    if res_sub2 not in translations:
                        translations[res_sub2] = res_meaning_2
    return translations

read_variables()
with open('fceri_meanings','w') as afile:
    translations = read_variables()
    #print(translations["S4"])
    list_of_strings = [ f'{key} : {translations[key]}' for key in translations]
    [afile.write(f'{st}\n') for st in list_of_strings ]
