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
        iteracio = 0
        m = len(b)
        n = len(A[0])
        if(n < m):
            return None, "Sistema incorrecte", [], None, None
        no_basiques = [a for a in range(n-m)]
        basiques = [a for a in range(n) if a not in no_basiques]
        basiques_noves = None
        if fase1 == None:
            nova_A = np.hstack((A, np.eye(m))) # horizontal stack
            nou_cost = np.array([0 for _ in range(n)] + [1 for _ in range(m)])
            x_f1, z_f1, basiques_noves, inv, iteracions = simplex(nou_cost, nova_A, b, np.eye(m), fase1 = True)
            inversa = inv
            iteracio = iteracions + 1
            if z_f1 in [None, "Infactible", "No acotat"] or np.round(z_f1,10) > 0:
                    doc.write("No hi ha solucio factible\n\n")
                    return None, "Infactible", [], None, iteracio
            if min(x_f1) == 0 and basiques_noves[np.argmin(x_f1)] >= n:
                # Cas en el que fase I acaba amb degeneració en una de les variables artificials.
                """
                Una alternativa és buscar qualsevol no bàsica amb cost reduït 0 i l'element de la
                qual en la fila de la variable artificial sigui diferent a 0.
                """
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
                    doc.write("SBF inicial contenia una variable artificial que no s'ha pogut treure.\n")
                    cost = np.append(cost, (np.zeros(m)))
                    z = np.dot(cost[basiques_noves], x_f1)
                    doc.write(f"x = {x_f1}\nz* = {z}\nbasiques = {basiques_noves}\n\n")
                    return x_f1, z, basiques_noves, None, iteracio
            basiques = basiques_noves
            no_basiques = [a for a in range(n) if a not in basiques]
            doc.write("SBF inicial trobada.\n\n")
            doc.write("Fase 2\n")
        else:
            doc.write("Fase 1\n")
        cost_b = cost[basiques]
        B_inv = inversa
        x = np.dot(B_inv, b)
        z = np.dot(cost_b, x)

        if m == n: # Si hi ha el mateix nombre de variables que restriccions, la SBF és la única
            x = np.dot(np.linalg.inv(A), b)
            z = np.dot(cost, x)
            return x, z, basiques, None, iteracio
        
        while True:
            degenerat = False
            B = A[:,basiques]
            A_n = A[:,no_basiques]
            B_inv = inversa
            cost_b = cost[basiques]
            cost_n = cost[no_basiques]
            if np.round(min(x),10) == 0:
                doc.write("Solucio degenerada\n")
                degenerat = True
                pass
            elif min(x) < 0:
                return x, None, None, None, iteracio
            
            # és optim?
            r = np.subtract(cost_n, np.dot(np.dot(cost_b, B_inv),A_n))
            if min(r) < 0:
                pass
            elif min(r) > 0:
                doc.write("Solucio optima trobada\n\n") if fase1 != True else doc.write("Fi Fase I\n\n")
                doc.write(f"x = {x}\nz* = {z}\nbasiques = {basiques}\nr = {r}\n\n") if fase1 != True else doc.write(f"basiques = {basiques}\n\n")
                doc.write(f"Nombre d'iteracions: {iteracio + 1}\n\n")
                return x, z, basiques, inversa, iteracio
            else:
                doc.write("Una de les solucions optimes trobada\n\n") if fase1 != True else doc.write("Fi Fase I\n\n")
                doc.write(f"x = {x}\nz* = {z}\nbasiques = {basiques}\nr = {r}\n\n") if fase1 != True else doc.write(f"basiques = {basiques}\n\n")
                doc.write(f"Nombre d'iteracions: {iteracio + 1}\n\n")
                return x, z, basiques, inversa, iteracio
            for e in range(len(r)):
                if r[e] < 0:
                    break
            entra = no_basiques[e]
            # direcció bàsica factible
            d_B = -np.dot(B_inv, A_n[:,e])

            if min(d_B) >= 0:
                doc.write("Optim no acotat (raig)\n\n")
                doc.write(f"basiques = {basiques}\ndB = {d_B}\n\n")
                doc.write(f"Nombre d'iteracions: {iteracio + 1}\n\n")
                return x, "No acotat", basiques, inversa, iteracio
            # longitud de pas
            theta, marxa = min([(np.divide((-x[i]),d_i), basiques[i]) for i, d_i in enumerate(d_B) if d_i < 0])

            if degenerat and theta == 0:
                """
                Aquest condicional està per la propietat ii de la diapositiva 30.
                No vam acabar d'entendre exactament en quins casos es donava aquest
                problema
                """
                if max([np.divide((-x[i]),d_i) for i, d_i in enumerate(d_B) if d_i < 0]) == 0:
                    # No hi ha cap theta > 0 (lo de la presentació?)
                    # return x, "Infactible", basiques, None
                    pass
                else:
                    pass
                    

            p = basiques.index(marxa)
            basiques[p] = entra
            no_basiques = [a for a in range(n) if a not in basiques] # Ordenades
            doc.write(f"q = {entra}, out = {marxa}, p = {p}, theta = {theta}, z = {z}\n")
            
            # actualitzacions
            transformacio = np.eye(m)
            transformacio[:,p] = [np.divide((-d_B[i]), d_B[p]) if i!=p else np.divide((-1),d_B[p]) for i in range(m)]
            inversa = np.dot(transformacio, B_inv)

            x += np.dot(theta,d_B)
            x[p] = theta

            z += np.dot(r[e],theta)
            iteracio += 1

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
        print(f"Iteracions: {a[4]}")


print("------------- \n Test \n --------------")
with open ("output.txt", "a") as doc:
    doc.write("------------------------------------------------------------------------------\n")
    doc.write("                                 Test\n")
    doc.write("------------------------------------------------------------------------------\n\n")
for alumne in range(67, 71):
    for problema in range(1, 5):
        print(f"Alumne {alumne}, problema {problema}")
        c, A, b, z ,v = read_dades(alumne,problema, fitxer="Datos_práctica_1_test.txt")
        with open ("output.txt", "a") as doc:
            doc.write(f"Alumne {alumne}, problema {problema}\n")
        a = simplex(c, A, b, None)
        print(a[1])
"""
Per comprovar:

from scipy.optimize import linprog
r = linprog(c, A_eq = A, b_eq = b, method='highs')
if r['status']==0:
    if a[1] != None and a[1] != "?":
        print(f"{r['fun']:.4f}", f"{a[1]:.4f}")
    else:
        print("a was none and r was ", r['fun'])
elif r['status'] == 2:
    print("infactible")
    if a[1] == "Infactible":
        print(True)
    else:
        print(False)
elif r['status'] == 3:
    print("raig")
    if a[1] == "No acotat":
        print(True)
    else:
        print(False)
"""