import numpy as np
import re

def read_dades(num: int, prob: int, fitxer: str = "OPT23-24_Datos práctica 1.txt"):
    with open(fitxer, "r") as data:
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

import time
def simplex(cost: np.array, A: np.array, b: np.array, inversa = None, fase1 = None):
    with open ("output.txt", "a") as doc:
        m = len(b)
        n = len(A[0])
        if(n < m):
            return None, "Sistema incorrecte", [], None, None
        no_basiques = [a for a in range(n-m)]
        basiques = [a for a in range(n) if a not in no_basiques]
        basiques_noves = None
        if fase1 == None:
            doc.write("Fase 2\n")
            nova_A = np.hstack((A, np.eye(m))) # horizontal stack
            nou_cost = np.array([0 for _ in range(n)] + [1 for _ in range(m)])
            x_f1, z_f1, basiques_noves, inv = simplex(nou_cost, nova_A, b, np.eye(m), fase1 = True)
            inversa = inv
            if z_f1 in [None, "Infactible", "No acotat"] or np.round(z_f1,10) > 0:
                    doc.write("No hi ha solucio factible\n\n")
                    return None, "Infactible", [], None, None
            if min(x_f1) == 0 and basiques_noves[np.argmin(x_f1)] >= n:
                # Cas en el que fase I acaba amb degeneració en una de les variables artificials.
                marxa = basiques_noves[np.argmin(x_f1)]
                basiques_per_iterar = [a for a in basiques if a not in basiques_noves]
                for idx in basiques_per_iterar:
                    try:
                        basiques_temp = basiques_noves.copy()
                        basiques_temp[np.argmin(x_f1)] = idx
                        B = A[:,basiques_temp]
                        inversa = np.linalg.inv(B)
                        basiques_noves[np.argmin(x_f1)] = idx
                        break
                    except:
                        continue
                else:
                    cost = np.append(cost, (np.zeros(m)))
                    z = np.dot(cost[basiques_noves], x_f1)
                    return x_f1, z, basiques_noves, None
            basiques = basiques_noves
            no_basiques = [a for a in range(n) if a not in basiques]
        else:
            doc.write("Fase 1\n")
        cost_b = cost[basiques]
        B_inv = inversa
        x = np.dot(B_inv, b)
        z = np.dot(cost_b, x)

        if m == n: # Si hi ha el mateix nombre de variables que restriccions, la SBF és la única
            x = np.dot(np.linalg.inv(A), b)
            z = np.dot(cost, x)
            return x, z, basiques, None
        
        while True:
            degenerat = False
            B = A[:,basiques]
            A_n = A[:,no_basiques]
            B_inv = inversa
            cost_b = cost[basiques]
            cost_n = cost[no_basiques]
            if min(x) == 0:
                doc.write("Solucio degenerada\n\n")
                degenerat = True
                pass
            elif min(x) < 0:
                return x, None, None, None
            
            # és optim?
            r = np.subtract(cost_n, np.dot(np.dot(cost_b, B_inv),A_n))
            if min(r) < 0:
                pass
            elif min(r) > 0:
                doc.write("Solucio optima trobada\n\n") if fase1 != True else doc.write("SBF inicial trobada\n\n")
                doc.write(f"x = {x}\nz* = {z}\nbasiques = {basiques}\n\n") if fase1 != True else doc.write(f"basiques = {basiques}\n\n")
                return x, z, basiques, inversa
            else:
                doc.write("Una de les solucions òptimes trobada\n\n") if fase1 != True else doc.write("SBF inicial trobada\n\n")
                doc.write(f"x = {x}\nz* = {z}\nbasiques = {basiques}\n\n") if fase1 != True else doc.write(f"basiques = {basiques}\n\n")
                return x, z, basiques, inversa
            for e in range(len(r)):
                if r[e] < 0:
                    break
            entra = no_basiques[e]
            # direcció bàsica factible
            d_B = -np.dot(B_inv, A_n[:,e])

            if min(d_B) >= 0:
                doc.write("Optim no acotat (raig)\n\n")
                return x, "No acotat", basiques, inversa
            # longitud de pas
            theta, p = min([(np.divide((-x[i]),d_i), i) for i, d_i in enumerate(d_B) if d_i < 0])

            if degenerat and theta == 0:
                if max([np.divide((-x[i]),d_i) for i, d_i in enumerate(d_B) if d_i < 0]) == 0:
                    # No hi ha cap theta > 0 (lo de la presentació?)
                    # return x, "Infactible", basiques, None
                    pass
                else:
                    pass
                    

            marxa = basiques[p]
            basiques[p] = entra
            no_basiques = [a for a in range(n) if a not in basiques] # Ordenades
            doc.write(f"iout: {entra}, q = {p}, theta = {theta}, z = {z}\n")
            
            # actualitzacions
            transformacio = np.eye(m)
            transformacio[:,p] = [np.divide((-d_B[i]), d_B[p]) if i!=p else np.divide((-1),d_B[p]) for i in range(m)]
            inversa = np.dot(transformacio, B_inv)

            x += np.dot(theta,d_B)
            x[p] = theta

            z += np.dot(r[e],theta)
for alumne in range(1, 67):
    for problema in range(1, 5):
        print(f"Alumne {alumne}, problema {problema}")
        c, A, b, z ,v = read_dades(alumne,problema)
        with open ("output.txt", "a") as doc:
            doc.write(f"Alumne {alumne}, problema {problema}\n")
        a = simplex(c, A, b, None)
        if z != None:
            result = f"{a[1]:.4f}" == f"{z:.4f}"
            print(result, f"{a[1]:.4f}")
        else:
            print(a[1])

print("------------- \n Test \n --------------")
for alumne in range(67, 71):
    for problema in range(1, 5):
        print(f"Alumne {alumne}, problema {problema}")
        c, A, b, z ,v = read_dades(alumne,problema, fitxer="Datos_práctica_1_test.txt")
        with open ("output.txt", "a") as doc:
            doc.write(f"Alumne {alumne}, problema {problema}\n")
        a = simplex(c, A, b, None)
        print(a[1])