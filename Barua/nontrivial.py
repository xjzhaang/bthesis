def nontrivial():
    f = open("2m_macrovariables.txt", "r")
    lines = f.readlines()
    new_lines = []
    for l in lines:
        if " + " in l:
            new_lines.append(l)
    return new_lines

with open('2m_nontrivial.txt','w') as afile:
    new_lines = nontrivial()
    [afile.write(f'{st}') for st in new_lines]