
import math
import numpy as np


class MiembroArmadura():
    def __init__(self, elemento, area, modElasticidad, coordenada_xi, coordenada_yi, coordenada_xf, coordenada_yf, vec_coord):
        self.elem = elemento
        self.A = area
        self.E = modElasticidad
        self.xi = coordenada_xi
        self.yi = coordenada_yi
        self.xf = coordenada_xf
        self.yf = coordenada_yf
        self.vc = vec_coord
        self.L = self.Longitud()
        self.k_loc = self.Rig_Loc()
        self.c = self.cos()
        self.s = self.sen()
        self.T = self.Transf()
        self.k_glob = self.Rig_Glob()

    def __str__(self):
        print("Área : ", self.A)
        print("Módulo de elasticidad : ", self.E)
        print("Coordenadas iniciales : ({}, {})".format(self.xi, self.yi))
        print("Coordenadas finales : ({}, {})".format(self.xf, self.yf))
        print("Vector de coordenadas : ", self.vc)
        print("Longitud de la barra : ", self.L)
        print("Matriz de rígidez de la barra : ", self.k_loc)
        print("Coseno : ", self.c)
        print("Seno : ", self.s)
        print("Matriz de transformación : ", self.T)
        print("")
        print("Matriz global de la barra : ", self.k_glob)
        return ""

    def Longitud(self):
        return math.sqrt((self.xf-self.xi)**2+(self.yf-self.yi)**2)

    def Rig_Loc(self):
        return (self.A*self.E/self.L)*np.array([[1, -1], [-1, 1]])

    def cos(self):
        return (self.xf-self.xi)/self.L

    def sen(self):
        return (self.yf-self.yi)/self.L

    def Transf(self):
        return np.array([[self.c, self.s, 0, 0], [0, 0, self.c, self.s]])

    def Rig_Glob(self):
        Tt = np.transpose(self.T)
        Tt_k_loc = np.matmul(Tt, self.k_loc)
        return np.matmul(Tt_k_loc, self.T)


class AnalisisMatricial():
    def __init__(self, tbl_Elem, tbl_Nods, tbl_Frza, tbl_Desp):
        self.tE = tbl_Elem
        self.tN = tbl_Nods
        self.tF = tbl_Frza
        self.tD = tbl_Desp
        self.nE = len(self.tE)
        self.nN = len(self.tN)
        self.nGl = self.nN * 2
        [self.nR, self.nFc, self.N] = self.VectorCoordenadasGlobales()
        self.PI = np.array(self.Matriz_PI())
        self.Armad = self.Armadura()
        self.Kg = self.MatrizRigidezGlobal()
        [self.k11, self.k12, self.k21, self.k22] = self.MatrizRigidezGlobalParciales()
        self.Fc = self.VecFrzasGlobalesConocidas()
        self.Dc = self.VecDespGlobalesConocidas()
        [self.Dd, self.Fd, self.Fg, self.Dg] = self.ReaccionesDesplazamientos()
        self.TC = self.TensionCompresion()

    def __str__(self):
        print("Cantidad de barras : ", self.nE)
        print("Cantidad de nodos : ", self.nN)
        print("Grados de libertad : ", self.nGl)
        print("Número de reacciones : ", self.nR)
        print("Número de gdl libres : ", self.nFc)
        print("Diccionario de nodos : ", self.N)
        print("Matriz PI : ", self.PI)
        print("Elemento 1: ", self.Armad[0])
        print("Matriz de rígidez global de la armadura : \n")
        print(self.Kg)
        print("")
        print("Matriz de rígidez k11: \n", self.k11)
        print("")
        print("Matriz de rígidez k12: \n", self.k12)
        print("")
        print("Matriz de rígidez k21: \n", self.k21)
        print("")
        print("Matriz de rígidez k22: \n", self.k22)
        print("")
        print("Vector de fuerzas conocidas: \n", self.Fc)
        print("")
        print("Vector de desplazamientos conocidas: \n", self.Dc)
        print("")
        print("Vector de desplazamientos desconocidas: \n", self.Dd)
        print("")
        print("Vector de fuerzas desconocidas: \n", self.Fd)
        print("")
        print("Vector de fuerzas globales: \n", self.Fg)
        print("")
        print("Vector de desplazamientos globales: \n", self.Dg)
        print("")
        print("Vector de fuerzas internas: \n", self.TC)
        return ""

    def VectorCoordenadasGlobales(self):
        numReacciones = 0
        DiccionarioNodos = {}
        c = 1
        gdl = self.nGl
        for i in range(self.nN):
            tipo = self.tN[i][3]
            NODO_key = self.tN[i][0]
            cx = self.tN[i][1]
            cy = self.tN[i][2]
            if tipo == "Libre":
                numReacciones += 0
                DiccionarioNodos.setdefault(NODO_key, [cx, cy, c, c + 1])
                c += 2
            elif tipo == "Fijo":
                numReacciones += 2
                DiccionarioNodos.setdefault(NODO_key, [cx, cy, gdl - 1, gdl])
                gdl -= 2
            elif tipo == "DX":
                numReacciones += 1
                DiccionarioNodos.setdefault(NODO_key, [cx, cy, c, gdl])
                c += 1
                gdl -= 1
            elif tipo == "DY":
                numReacciones += 1
                DiccionarioNodos.setdefault(NODO_key, [cx, cy, gdl, c])
                c += 1
                gdl -= 1
        numFuerzasConocidas = self.nGl - numReacciones
        return numReacciones, numFuerzasConocidas, DiccionarioNodos

    def Matriz_PI(self):
        PI = []
        for i in range(self.nE):
            ni = self.tE[i][3]
            nf = self.tE[i][4]
            PI.append([self.N[ni][2], self.N[ni][3], self.N[nf][2], self.N[nf][3]])
        return PI

    def Armadura(self):
        Elem = []
        for i in range(self.nE):
            el = self.tE[i][0]
            a = self.tE[i][1]
            me = self.tE[i][2]
            NI = self.tE[i][3]
            NF = self.tE[i][4]
            xi = self.N[NI][0]
            yi = self.N[NI][1]
            xf = self.N[NF][0]
            yf = self.N[NF][1]
            Clix = self.N[NI][2]
            Cliy = self.N[NI][3]
            Clfx = self.N[NF][2]
            Clfy = self.N[NF][3]
            Elem.append(MiembroArmadura(el, a, me, xi, yi, xf, yf, [Clix, Cliy, Clfx, Clfy]))
        return Elem

    def MatrizRigidezGlobal(self):
        K = np.zeros((self.nGl, self.nGl))
        for e in range(self.nE):
            ke_global = self.Armad[e].k_glob
            for i in range(4):
                for j in range(4):
                    a = self.PI[e, i] - 1
                    b = self.PI[e, j] - 1
                    K[a, b] = ke_global[i, j] + K[a, b]
        return K

    def MatrizRigidezGlobalParciales(self):
        k11 = self.Kg[0:self.nFc, 0:self.nFc]
        k12 = self.Kg[0:self.nFc, self.nFc:self.nGl]
        k21 = self.Kg[self.nFc:self.nGl, 0:self.nFc]
        k22 = self.Kg[self.nFc:self.nGl, self.nFc:self.nGl]
        return k11, k12, k21, k22

    def VecFrzasGlobalesConocidas(self):
        Fc = np.zeros((self.nFc, 1))
        for i in range(len(self.tF)):
            direccion = self.tF[i][2]
            nodo = self.tF[i][1]
            if direccion == "DX":
                j = self.N[nodo][2] - 1
            elif direccion == "DY":
                j = self.N[nodo][3] - 1
            Fc[j] = self.tF[i][0]
        return Fc

    def VecDespGlobalesConocidas(self):
        Dc = np.zeros((self.nGl - self.nFc, 1))
        for i in range(len(self.tD)):
            direccion = self.tD[i][2]
            nodo = self.tD[i][1]
            if direccion == "DX":
                j = self.N[nodo][2] - 1 - self.nFc
            elif direccion == "DY":
                j = self.N[nodo][3] - 1 - self.nFc
            Dc[j] = self.tD[i][0]
        return Dc

    def ReaccionesDesplazamientos(self):
        k11_inv = np.linalg.inv(self.k11)
        Dd = np.matmul(k11_inv, (self.Fc - np.matmul(self.k12, self.Dc)))
        Fd = np.matmul(self.k21, Dd) + np.matmul(self.k22, self.Dc)
        Fg = np.concatenate((self.Fc, Fd), axis=0)
        Dg = np.concatenate((Dd, self.Dc), axis=0)
        return Dd, Fd, Fg, Dg

    def TensionCompresion(self):
        TC = []
        for e in range(self.nE):
            Te = self.Armad[e].T
            ke_loc = self.Armad[e].k_loc
            dix = self.Dg[self.Armad[e].vc[0] - 1]
            diy = self.Dg[self.Armad[e].vc[1] - 1]
            dfx = self.Dg[self.Armad[e].vc[2] - 1]
            dfy = self.Dg[self.Armad[e].vc[3] - 1]
            D = np.array([dix, diy, dfx, dfy])
            ke_loc_Te = np.matmul(ke_loc, Te)
            ke_loc_Te_D = np.matmul(ke_loc_Te, D)
            TC.append(round(float(ke_loc_Te_D[1]), 2))
        return TC