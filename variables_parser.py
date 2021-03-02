def read_variables():
    #to 1479
    f = open("fceri_ode.txt", "r")
    lines = f.readlines()
    translations = {"S0": "IgE", "S1": "Lyn", "S2": "Syk", "S3": "r"}
    for l in lines:
        substrates = l.split("->")[0].strip()
        if " + " in substrates:
            sub1 = substrates.split(" + ")[0]
            sub2 = substrates.split(" + ")[1]
            two_sub = True   
        else:
            sub1 = substrates
        res_sub = l.split("->")[1].split(" ,")[0].strip()
        rate = l.split("->")[1].split(" , ")[1].strip()
        if two_sub:
            res_meaning = translations[sub1] + " " + translations[sub2]
            pass
        else:
            pass
        print(rate)
    
read_variables()
