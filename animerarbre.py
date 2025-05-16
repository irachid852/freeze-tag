# -*- coding: utf-8 -*-
"""
Created on Fri May 16 10:16:21 2025

@author: ismail rachid
"""

import freeze_tag_mini as ft
import numpy as np
import matplotlib.pyplot as plt
from random import randrange,uniform
import gc
ft.SEED = -1
ft.PROCESS_TIME = 0


def draw_points(P, c='gray', s=5, m='o'):
    """
    Draw a set of points P, each with colour c, size s, and shape m, in the current plt drawing. For a simple display of the points, a following extra plt.show() or draw_all() is required.
    """
    for (x,y) in P:
        plt.plot([x], [y], color = c, marker = m, ms = s, ls = 'none')
        
        
lio= [[0.0, 1], [0.785, 1], [1.571, 1], [2.356, 1],
 [3.142, 1], [3.927, 1], [4.712, 1], [5.498, 1]]
li= [[0,0],[1.0, 0.0],[0.6234898018587336, 0.7818314824680298],[-0.22252093395631434, 0.9749279121818236],[-0.900968867902419, 0.43388373911755823],[-0.9009688679024191, -0.433883739117558],[-0.2225209339563146, -0.9749279121818236],[0.6234898018587334, -0.7818314824680299],[0.5,0.25]]
octo = [[0,0],[1.0, 0.0],
 [0.7071, 0.7071],
 [0.0, 1.0],
 [-0.7071, 0.7071],
 [-1.0, 0.0],
 [-0.7071, -0.7071],
 [0.0, -1.0],
 [0.7071, -0.7071]]
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


# Exemple :
'''arbre = [0, [3, [4, [5, [6]], [7]], [1, [2], [8]]]]
print(enfants_par_robot(arbre))'''

def etape(P,v):
    x, T = ft.optimal_tree(0,P,ft.dist_L2)
    arbre_chemin= obtenir_chemins_robots(T)
    dic = {}
    print(arbre_chemin)
    ls = [False for i in range(9)]
    for i, el in enumerate(P):
        dic[f'Robot {i}'] = Noeud(el,None,False)
    for i in range(9):
        dic[f'Robot {i}'].suivant = [dic[f'Robot {j}'] for j in arbre_chemin[i]]
    dic['Robot 0'].reveil = True
    while ls != [True for i in range(9)]:
        for i in range(9):
            ls[i]= dic[f'Robot {i}'].aller(v)
        lst = [dic[f'Robot {i}'].coord for i in range(9)]
        #x = [dic[f'Robot {i}'].coord[0] for i in range(9)]
        #y = [dic[f'Robot {i}'].coord[1] for i in range(9)]
        plt.clf()
        plt.axis('equal')
        ft.draw_disc()
        ft.draw_points(lst,'red',10)
        #plt.savefig(f'{idx}.png')
        
        plt.show()
        plt.close('all')
        gc.collect()
    print(ls)

from PIL import Image
import matplotlib.pyplot as plt
import io

def etape(P, v, gifname='freeze_tag.gif'):
    x, T = ft.optimal_tree(0, P, ft.dist_L2)
    ft.draw_all('optimal', x, T)
    arbre_chemin = obtenir_chemins_robots(T)
   
    dic = {}
    
    ls = [False] * len(P)
    for i, el in enumerate(P):
        dic[f'Robot {i}'] = Noeud(el, None, False)
    for i in range(len(P)):
        if not i in arbre_chemin.keys(): arbre_chemin[i] = []
        dic[f'Robot {i}'].suivant = [dic[f'Robot {j}'] for j in arbre_chemin[i] ]
    dic['Robot 0'].reveil = True

    frames = []

    while ls != [True] * len(P):
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

def pol2car(ls):
    lsc = []
    for (a,r) in ls:
        lsc.append([r*np.cos(a),r*np.sin(a)])
    return lsc
la = [[0,0]]+pol2car([[ uniform(0,2*np.pi),uniform(0.5,1)] for i in range(12)])
etape(octo,0.05,'freeze_tag.gif')
            




# Lecture et analyse de l'entrée



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
