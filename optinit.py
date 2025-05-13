# -*- coding: utf-8 -*-
"""
Created on Tue May 13 10:51:30 2025

@author: ismail rachid
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 13 10:21:02 2025

@author: ismail rachid
"""

import matplotlib.pyplot as plt
from random import randrange
import freeze_tag_mini as ft
from itertools import combinations
import numpy as np
import gc

def ptposs(ep):
    ls =[]
    max_i = int(1 / ep)
    lx = [i*ep for i in range(-max_i, max_i + 1)]
    for x in lx:
        for y in lx:
            if x**2 + y**2 <= 1:
                ls.append([x, y])
    return ls

def generer(ep):
    ls = []
    max_i = int(1 / ep)
    for _ in range(8):
        while True:
            xi = randrange(-max_i, max_i + 1)
            yi = randrange(-max_i, max_i + 1)
            x = xi * ep
            y = yi * ep
            if x**2 + y**2 <= 1:
                ls.append([x, y])
                break
    return ls

ft.SEED = -1
ft.PROCESS_TIME = 0

def homothetie(pts):
    l =np.sqrt( max((x**2 + y**2) for (x, y) in pts))
    return [(x/l, y/l) for (x, y) in pts]

somd = []
m = 0
t = [None]
pts=[]
lsv = ptposs(0.01)
lf = [[0,0]]
lsv.remove([0,0])
lfc =[]
lsd = set()
li= [[1.0, 0.0],[0.6234898018587336, 0.7818314824680298],[-0.22252093395631434, 0.9749279121818236],[-0.900968867902419, 0.43388373911755823],[-0.9009688679024191, -0.433883739117558],[-0.2225209339563146, -0.9749279121818236],[0.6234898018587334, -0.7818314824680299],[0.5,0.25]]
dic ={}

class Noeud:
    def __init__(self,combinaison=homothetie(li),fils=None,parent=None):
        self.combinaison = combinaison
        self.fils = [] if fils is None else fils
        self.maxi , _ =  ft.optimal_tree(0,[[0,0]]+self.combinaison,ft.dist_L2)
        self.parent = parent

    def suivant(self,ind):
        for i in lsv:
            if 1==1 and \
               (i not in [[0,0]]+self.combinaison) and \
               0<((i[1]-self.combinaison[ind][1])**2 +(i[0]-self.combinaison[ind][0])**2) <4*0.01**2:

                config = frozenset(map(tuple, self.combinaison[:ind] + [i] + self.combinaison[(ind+1):]))
                if config not in lsd:
                    ls = [[0,0]] + homothetie(self.combinaison[:ind]+[i]+self.combinaison[(ind+1):])
                    x, _ = ft.optimal_tree(0, ls, ft.dist_L2)
                    lsd.add(config)
                    if 1==1:#self.maxi < x:
                        self.fils.append(Noeud(self.combinaison[:ind]+[i]+self.combinaison[(ind+1):], parent=self))

    def cons(self):
        for ind in range(len(self.combinaison)):
            self.suivant(ind)

from time import time 
import heapq



class FilePrioriteMax:
    def __init__(self):
        self.file = []
        self.compteur = 0  # pour gérer les égalités

    def ajouter(self, noeud):
        heapq.heappush(self.file, (-noeud.maxi, self.compteur, noeud))
        self.compteur += 1

    def retirer(self):
        return heapq.heappop(self.file)[2]

    def est_vide(self):
        return len(self.file) == 0

def trouver(racine,seuil,option='seuil'):
    racine.cons()
    pile = racine.fils
    meilleur = racine
    t0 = time()
    while pile and time()-t0<1000:
        courant = pile.pop()
        courant.cons()
        print(courant.maxi)

        if courant.maxi > meilleur.maxi:
            print(f'nouveau max : {courant.maxi}')
            meilleur = courant

        if option == 'seuil' and meilleur.maxi > seuil:
            return meilleur
        elif option == 'local' and abs(meilleur.maxi - courant.maxi) < seuil:
            return meilleur

        pile.extend(courant.fils)
        courant.fils = []
        gc.collect()

    return meilleur

def trouver(racine, seuil, option='seuil'):
    racine.cons()
    file = FilePrioriteMax()
    for f in racine.fils:
        file.ajouter(f)

    meilleur = racine
    t0 = time()

    while not file.est_vide() and time() - t0 < 500:
        courant = file.retirer()
        courant.cons()
        print(courant.maxi)

        if courant.maxi > meilleur.maxi:
            print(f'nouveau max : {courant.maxi}')
            meilleur = courant

        if option == 'seuil' and meilleur.maxi > seuil:
            return meilleur
        elif option == 'local' and abs(meilleur.maxi - courant.maxi) < seuil:
            return meilleur

        for f in courant.fils:
            file.ajouter(f)

    return meilleur

rac = Noeud()
x, T = ft.optimal_tree(0, [[0,0]]+rac.combinaison, ft.dist_L2)
ft.draw_all('optimal', x, T)

rep = trouver(rac, 3.48)
x, T = ft.optimal_tree(0, [[0,0]]+rep.combinaison, ft.dist_L2)
ft.draw_all('optimal', x, T)

noeudac = rep
evolution = []
while noeudac:
    evolution.append([[0,0]] + homothetie(noeudac.combinaison))
    noeudac = noeudac.parent
evolution.reverse()

image_files = []
for idx, evo in enumerate(evolution):
    x, T = ft.optimal_tree(0, evo, ft.dist_L2)
    plt.clf()
    ft.draw_all('optimal', x, T, save=str(idx)+'.png')
    plt.close('all')
    gc.collect()
