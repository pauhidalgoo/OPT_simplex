import numpy as np
import re

def read_dades(num: int, prob: int):
    with open("OPT23-24_Datos práctica 1.txt", "r") as data:
        line = data.readline()
        while line != None:
            while not("datos"+ str(num) in line.replace(" ", "") and "problemaPL" + str(prob) in line.replace(" ", "")):
                line = data.readline()
            _ = data.readline()
            cost = []
            while "A=" not in line:
                pattern = r'^[0-9\- ]*$'
                if re.match(pattern, line):
                    cost += list(map(int, re.findall(r'-?\d+', line)))
                    line = data.readline()
                else:
                    line = data.readline()
            A = []
            firstcolumn = True
            while "b=" not in line:
                line = data.readline()
                pattern = r'^[0-9\- ]*$'
                if re.match(pattern, line) and line != "\n":
                    A.append( list(map(int, re.findall(r'-?\d+', line))))
                
                if "Column" in line:
                    if firstcolumn:
                        firstcolumn = False
                    else:
                        _ = data.readline()
                        for i in range(len(A)):
                            line = data.readline()
                            A[i] = A[i]+list(map(int, re.findall(r'-?\d+', line)))
            
            b = []
            line = data.readline()
            b += list(map(int, re.findall(r'-?\d+', line)))
            while "z*=" not in line or "---" in line:
                line = data.readline()
                if "---" in line:
                    z = None
                    v  = None
                    break
            else:
                line = data.readline()
                z = float(line)
                line = data.readline()
                line = data.readline()
                line = data.readline()
                v = list(map(int, re.findall(r'-?\d+', line)))
            return np.array(cost), np.array(A), np.array(b), z, v
        
c, A, b, z ,v = read_dades(1,1)

import time
def simplex(cost: np.array, A: np.array, b: np.array, z = None, v = None, inversa = None, fase1 = None):
    z_original = z
    m = len(b)
    n = len(A[0])
    no_basiques = [a for a in range(n-m)]
    basiques = [a for a in range(n) if a not in no_basiques]
    basiques_noves = None
    if fase1 == None:
        nova_A = np.hstack((A, np.eye(m))) # horizontal stack
        nou_cost = np.array([0 for _ in range(n)] + [1 for _ in range(m)])
        _, z_f1, basiques_noves, _, inv = simplex(nou_cost, nova_A, b,None,None, np.eye(m), fase1 = True)
        inversa = inv
        if z_f1 == None or np.round(z_f1,10) > 0:
            print("Infactible")
            return None, None, [], None, None
        
        basiques = basiques_noves
        no_basiques = [a for a in range(n) if a not in basiques]

    cost_b = cost[basiques]
    B_inv = inversa
    x = np.dot(B_inv, b)
    z = np.dot(cost_b, x)
    while True:
        degenerat = False
        B = A[:,basiques]
        A_n = A[:,no_basiques]
        B_inv = inversa
        cost_b = cost[basiques]
        cost_n = cost[no_basiques]
        if min(x) == 0:
            print("Degenerat")
            degenerat = True
            pass
        elif min(x) < 0:
            print("Error")
            return None, None, None, None, None
        # és optim?
        # és el MÉS ÒPTIM?
        r = np.subtract(cost_n, np.dot(np.dot(cost_b, B_inv),A_n))
        if min(r) < 0:
            #print("No és òptim burru")
            pass
        else:
            print("optim")
            return x, z, basiques, z_original, inversa
        for e in range(len(r)):
            if r[e] < 0:
                break
        entra = no_basiques[e]
        # direcció bàsica factible
        d_B = -np.dot(B_inv, A_n[:,e])

        if min(d_B) >= 0:
            print("Raig extrem no acotat😈")
            return x, float("-inf"), basiques, z_original, inversa
        # longitud de pas
        theta, p = min([(np.divide((-x[i]),d_i), i) for i, d_i in enumerate(d_B) if d_i < 0], key=lambda x: x[0])

        if degenerat and theta == 0:
            print("Insatisfier")
            return None, None, None, None, None

        marxa = basiques[p]
        basiques[p] = entra
        no_basiques[no_basiques.index(entra)] = marxa
        # actualitzacions
        transformacio = np.eye(m)
        transformacio[:,p] = [np.divide((-d_B[i]), d_B[p]) if i!=p else np.divide((-1),d_B[p]) for i in range(m)]
        inversa = np.dot(transformacio, B_inv)

        x += np.dot(theta,d_B)
        x[p] = theta

        z += np.dot(r[e],theta)

        print(f"q = {entra}, p={marxa}, ")


for alumne in range(1, 45):
    for problema in range(1, 5):
        c, A, b, z ,v = read_dades(alumne,problema)
        a = simplex(c, A, b, z, v, None)
        if z != None:
            result = f"{a[1]:.4f}" == f"{z:.4f}"
            print(result, f"{a[1]:.4f}")

def print_simplex(simplex):
    print("\n Solució: \n")
    x, z, basiques, z_original, _= simplex
    str_b = ""
    for basica in basiques:
        str_b +=str(basica) + " "
    print("z* = ", z)
    print("vb* =",str_b)
    if z_original != None:
        if float("%.4f" % z) == z_original:
            print("La solució és òptima")
            print("""        .
       ,O,
      ,OOO,
'oooooOOOOOooooo'
  `OOOOOOOOOOO`
    `OOOOOOO`
    OOOO'OOOO
   OOO'   'OOO
  O'         'O""")
        else:
            print("La solució no és òptima")

print()
print_simplex(simplex(c, A, b, z, v))

            
            
                    
        
            

