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

def simplex(c, A, b, z, v):
    m, n = A.shape
    B = np.eye(m)
    N = np.eye(n)
    A = np.hstack((A, B))
    c = np.hstack((c, np.zeros(m)))
    c = -c
    v = np.array(v)
    while True:
        B = np.eye(m)
        N = np.eye(n)
        A = np.hstack((A, B))
        c = np.hstack((c, np.zeros(m)))
        c = -c
        v = np.array(v)
        cN = c
        cB = np.zeros(m)
        A = np.hstack((A, B))
        c = np.hstack((c, np.zeros(m)))
        c = -c
        v = np.array(v)
        cN = c
        cB = np.zeros(m)
        while True:
            cN = c[:n]
            cB = c[n:]
            B = A[:, n:]
            N = A[:, :n]
            Binv = np.linalg.inv(B)
            y = cB @ Binv
            cN = cN - y @ N
            if np.all(cN <= 0):
                break
            k = np.argmax(cN)
            d = Binv @ A[:, k]
            if np.all(d <= 0):
                return "El problema no tiene solución óptima"
            theta = np.array([b[i] / d[i] if d[i] > 0 else np.inf for i in range(m)])
            l = np.argmin(theta)
            x = np.zeros(n)
            x[k] = theta[l]
            x = x @ Binv
            x = np.hstack((x, np.zeros(m)))
            A = A - np.outer(A[:, k], x)
            c = c - c[k] * x
            b = b - A[:, k] * x
            z = z - c[k] * x
            v = v - c[k] * x
            if np.all(v >= 0):
                return z, v

print(simplex(c, A, b, z, v))




            
            
                    
        
            

