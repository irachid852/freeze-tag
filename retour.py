# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 11:34:12 2025

@author: racni
"""

import freeze_tag_mini as ft
import numpy as np
import matplotlib.pyplot as plt
from random import randrange,uniform
import gc
import copy as cp
ft.SEED = -1
ft.PROCESS_TIME = 0
def pol2car(ls):
    
    return [np.cos(ls[0])*ls[1],np.sin(ls[0])*ls[1]]

def draw_points(P, c='gray', s=5, m='o'):
    """
    Draw a set of points P, each with colour c, size s, and shape m, in the current plt drawing. For a simple display of the points, a following extra plt.show() or draw_all() is required.
    """
    for (x,y) in P:
        plt.plot([x], [y], color = c, marker = m, ms = s, ls = 'none')
class Noeud:
    def __init__(self,coord,suivant,reveil):
        self.coord = np.array(coord,dtype=float)
        self.suivant = suivant
        self.reveil = reveil
    def aller(self,v):
        if self.reveil == True and len(self.suivant)>0:
            if np.linalg.norm(self.coord-self.suivant[0].coord)**2<=2*v**2:
                self.coord = self.suivant[0].coord
                self.suivant[0].reveil = True
                self.suivant.pop(0)
                
            else:
                self.coord[1]+=v*np.sign(self.suivant[0].coord[1]-self.coord[1])
                self.coord[0]+=v*np.sign(self.suivant[0].coord[0]-self.coord[0])
        elif len(self.suivant) == 0:
            return True
        return False
    def aller(self, v):
        if self.reveil and len(self.suivant) > 0:
            cible = self.suivant[0]
            direction = cible.coord - self.coord
            distance = np.linalg.norm(direction)
            
            if distance <= v:
                self.coord = cible.coord.copy()
                cible.reveil = True
                self.suivant.pop(0)
            else:
                self.coord += v * direction / distance
        elif len(self.suivant) == 0:
            return True
        return False
                

from collections import deque

def obtenir_chemins_robots(arbre):
    from collections import deque
    file = deque()
    chemins_robots = {}

    if not arbre:
        return chemins_robots

    racine_valeur = arbre[0]
    file.append((0, [], arbre))  # (id_robot, chemin_courant, sous_arbre)

    while file:
        id_robot, chemin_courant, sous_arbre = file.popleft()
        valeur_courante = sous_arbre[0]
        nouveau_chemin = chemin_courant + [valeur_courante]

        enfants = sous_arbre[1:]

        if enfants:
            premier_enfant = enfants[0]
            file.append((id_robot, nouveau_chemin, premier_enfant))

            for enfant in enfants[1:]:
                nouveau_id_robot = valeur_courante
                nouveau_chemin_robot = [valeur_courante]
                file.append((nouveau_id_robot, nouveau_chemin_robot, enfant))
        else:
            if id_robot in chemins_robots:
                chemin_existant = chemins_robots[id_robot]
                chemin_existant.extend(nouveau_chemin[len(chemin_existant):])
            else:
                chemins_robots[id_robot] = nouveau_chemin

    return chemins_robots

class File:
    def __init__(self):
        self.elements = []

    def enfiler(self, elem):
        self.elements.append(elem)  # ajoute à la fin

    def defiler(self):
        if not self.est_vide():
            return self.elements.pop(0)  # retire le premier
        else:
            raise IndexError("File vide")
    
    def premier(self):
        if not self.est_vide():
            return self.elements[0]
        else:
            raise IndexError("File vide")
    
    def est_vide(self):
        return len(self.elements) == 0

    def taille(self):
        return len(self.elements)

def temps_rev(lsr):
    maxi = 0 
    maxval = None
    for r in lsr:
        dist = 0
        elt = np.array(r.coord)
        for e in r.suivant:
            val = np.array(pol2car(e.coord))
            dist+=np.linalg.norm(val-elt)
            elt = val[:]
        if dist>maxi:
            maxi = dist
            maxval = r
    return maxi, maxval

def temps_rev(lsr):
    maxi = 0 
    maxval = None
    for r in lsr:
        dist = 0
        elt = np.array(r.coord)  # cartésien
        for e in r.suivant:
            val = np.array(e.coord)  # aussi cartésien
            dist += np.linalg.norm(val - elt)
            elt = val
        if dist > maxi:
            maxi = dist
            maxval = r
    return maxi, maxval

from itertools import permutations

def retour(dico):
    liste = [ Noeud(i.coord[:],i.suivant[:],i.reveil) for i in list(dico.values())]
    
    print(len(liste))
    mini = np.inf
    valinf = None
    print(len(list(permutations(liste))))
    for i, el in enumerate(permutations(liste)):
        #print(i)
        for i, k in enumerate(dico.keys()):
            dico[k].suivant.append(el[i])
        v,_ = temps_rev(liste)
        if v<mini :
            mini = v
            valinf = el[:]
        for i, k in enumerate(dico.keys()):
            dico[k].suivant.pop(-1)
    for i, k in enumerate(dico.keys()):
        dico[k].suivant.append(valinf[i])
    
    print("fini")
            

from PIL import Image
import matplotlib.pyplot as plt
import io
import matplotlib.patches as patches

def etape(P, v, gifname='freeze_tag.gif'):
    x, T = ft.optimal_tree(0, P, ft.dist_L2)
    #ft.draw_all('optimal', x, T)
    arbre_chemin = obtenir_chemins_robots(T)
   
    dic = {}
    
    ls = [False] * len(P)
    for i, el in enumerate(P):
        dic[f'Robot {i}'] = Noeud(el, None, False)
    for i in range(len(P)):
        if not i in arbre_chemin.keys(): arbre_chemin[i] = []
        dic[f'Robot {i}'].suivant = [dic[f'Robot {j}'] for j in arbre_chemin[i] ]
    retour(dic)
    ls = list(dic.values())
    dic['Robot 0'].reveil = True
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    #plt.axis('equal')
    ft.draw_disc()
    mk, f = temps_rev(ls)
    print(mk)
    for r in ls:
        x, y = r.coord
        ax.plot(x, y, 'o')  # point du robot
        ax.text(x + 0.01, y + 0.01, f"{ls.index(r)}", fontsize=8)
    
        for i in range(len(r.suivant)-1):
            x1, y1 = r.suivant[i].coord
            x2, y2 = r.suivant[i+1].coord
            if r != f:
                ax.annotate("",
                    xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="blue",linewidth=3)
                )
            else:
                ax.annotate("",
                    xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="green",linewidth=3)
                )
    
    plt.show()
    frames = []

    while ls != [True] * len(P) and [i.suivant for i in dic.values()] !=[[] for i in range(len(list(dic.values())))]:
        if ls != [True] * len(P):
            for i in range(len(P)):
                ls[i] = dic[f'Robot {i}'].aller(v)
            lst = [dic[f'Robot {i}'].coord for i in range(len(P))]
    
            
        else:
            for i in range(len(P)):
                ls[i] = dic[f'Robot {i}'].aller(v)
            lst = [dic[f'Robot {i}'].coord for i in range(len(P))]
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ft.draw_disc()
        dicc= {True:'green', False:'red'}
        for r in dic.keys():
            plt.plot([dic[r].coord[0]], [dic[r].coord[1]],color = dicc[dic[r].reveil] , marker = 'o', ms = 7, ls = 'none')
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        frames.append(Image.open(buf).convert('RGB'))
        plt.close(fig)
    frames[0].save(gifname, save_all=True, append_images=frames[1:], duration=100, loop=0)
    print(f"GIF enregistré : {gifname}")

def pol2carl(ls):
    lsc = []
    for (a,r) in ls:
        lsc.append([r*np.cos(a),r*np.sin(a)])
    return lsc
octo = [[0,0],[1.0, 0.0],
 [0.7071, 0.7071],
 [0.0, 1.0],
 [-0.7071, 0.7071],
 [-1.0, 0.0],
 [-0.7071, -0.7071],
 [0.0, -1.0],
 [0.7071, -0.7071]]
la = [[0,0]]+pol2carl([[ uniform(0,2*np.pi),uniform(0.5,1)] for i in range(12)])
pts = [[0,0]]+[(np.cos(2*np.pi*k/10), np.sin(2*np.pi*k/10)) for k in range(10)]
etape(pts,0.05,'freeze_tag.gif')
