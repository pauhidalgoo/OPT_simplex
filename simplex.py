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
                print(line)
                z = float(line)
                line = data.readline()
                line = data.readline()
                line = data.readline()
                v = list(map(int, re.findall(r'-?\d+', line)))
            return np.array(cost), np.array(A), np.array(b), z, v
c, A, b, z ,v = read_dades(1,2)
import time
def simplex(cost: np.array, A: np.array, b: np.array, z = None, v = None, inversa = None):
    m = len(b)
    n = len(A[0])
    print(m,n)
    no_basiques = [a for a in range(n-m)]
    basiques = [a for a in range(n) if a not in no_basiques]
    basiques_noves = None
    while True:
        time.sleep(1)
        if basiques_noves != None:
            basiques = basiques_noves
            no_basiques = [a for a in range(n) if a not in basiques]
            basiques_noves = None
        B = A[:,basiques]
        A_n = A[:,no_basiques]
        B_inv = np.linalg.inv(B)
        x = np.dot(B_inv, b)
        cost_b = cost[basiques]
        cost_n = cost[no_basiques]
        z = np.dot(cost_b, x)
        if min(x) < 0:
            print(B)
            print(x)
            print("No és SBF burru, fes Fase I")
            nova_A = np.hstack((A, np.eye(m))) # horizontal stack
            nou_cost = np.array([0 for _ in range(n)] + [1 for _ in range(m)])
            basiques_noves, _ = simplex(nou_cost, nova_A, b)
            basiques_noves = np.nonzero(basiques_noves)
            continue
        elif min(x) == 0:
            #print("Degenerat")
            pass
        # és optim?
        # és el MÉS ÒPTIM?
        r = cost_n - np.dot(np.dot(cost_b, B_inv),A_n)
        if min(r) < 0:
            #print("No és òptim burru")
            pass
        else:
            return x, z
        for e in range(len(r)):
            if r[e] < 0:
                break
        entra = no_basiques[e]
        # direcció bàsica factible
        d_B = -np.dot(B_inv, A_n[:,e])
        # longitud de pas
        theta, marxa = min([((-x[i])/ d_i, i) for i, d_i in enumerate(d_B) ], key=lambda x: x[0])
        marxa = basiques[marxa]
        basiques[basiques.index(marxa)] = entra
        no_basiques.remove(entra)
        no_basiques.append(marxa)
        no_basiques = sorted(no_basiques)
        print(z, end="\r")
    
print(simplex(c, A, b, z, v))




            
            
                    
        
            

