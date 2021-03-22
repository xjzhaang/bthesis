def checkdiff():
    f = open("2m_macrovariables1.txt", "r")
    lines = f.readlines()
    lin = []
    for l in lines:
        name = l.split(" = ")[1].strip()
        lin.append(name)
    f1 = open("diff.txt", "r")
    lines1 = f1.readlines()
    lin1 = []
    for l in lines1:
        name = l.split(" = ")[1].strip()
        lin1.append(name)
    ll = list(set(lin1) - set(lin))
    print(ll)

checkdiff()
