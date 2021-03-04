def read_variables():
    #to 1479
    n = 0
    f = open("fceri_ode.txt", "r")
    lines = f.readlines()
    translations = {"S0": "IgE", "S1": "Lyn", "S2": "Syk", "S3": "r"}
    for l in lines:
        n += 1
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
                    res_meaning = translations[sub1] + ", " + translations[sub2]
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "pkp2":
                    res_meaning = translations[sub1] + ", " + translations[sub2]
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
                if rate == "pkpL":
                    if translations[sub1] == "Lyn":
                        res_meaning = translations[sub2] + ", cLyn"
                        if res_sub not in translations:
                            translations[res_sub] = res_meaning
                    else:
                        res_meaning = translations[sub1] + ", cLyn"
                        if res_sub not in translations:
                            translations[res_sub] = res_meaning
                if rate == "pkpLs":
                    if translations[sub1] == "Lyn":
                        res_meaning = translations[sub2] + ", SH2Lyn"
                        if res_sub not in translations:
                            translations[res_sub] = res_meaning
                    else:
                        res_meaning = translations[sub1] + ", SH2Lyn"
                        if res_sub not in translations:
                            translations[res_sub] = res_meaning
                if rate == "pkpS":
                    res_meaning = translations[sub1] + ", " + translations[sub2]
                    if res_sub not in translations:
                        translations[res_sub] = res_meaning
        if not two_res_sub:
            if sub1 in translations and sub2 in translations:
                if rate == "ppLb":
                    res_meaning = translations[sub1] + ", cbetaP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLg":
                    res_meaning = translations[sub1] + ", cgammaP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLbs":
                    res_meaning = translations[sub1] + ", SH2betaP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLgs":
                    res_meaning = translations[sub1] + ", SH2gammaP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLS":
                    res_meaning = translations[sub1] + ", cSykP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppLSs":
                    res_meaning = translations[sub1] + ", SH2SykP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppSS":
                    res_meaning = translations[sub1] + ", NotLoopSykP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
                if rate == "ppSSs":
                    res_meaning = translations[sub1] + ", LoopSykP"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning
        if two_res_sub:
            if sub1 in translations:
                if rate == "ppkm1":
                    res_meaning_1 = translations[sub1].replace(", r", "")
                    res_meaning_2 = "r"
                    if res_sub1 not in translations:
                        translations[res_sub1] = res_meaning_1
                    if res_sub2 not in translations:
                        translations[res_sub2] = res_meaning_2
                if rate == "pkm2":
                    res_meaning_2 = translations[sub1].replace(translations[res_sub1] + ", ", "")
                    if res_sub2 not in translations:
                        translations[res_sub2] = res_meaning_2
    print(translations)
read_variables()
