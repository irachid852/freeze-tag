# -*- coding: utf-8 -*-
"""
Created on Tue May 13 15:43:40 2025

@author: ismail rachid
"""
import freeze_tag_mini as ft
import numpy as np
import matplotlib.pyplot as plt
from random import randrange,uniform
ft.SEED = -1
ft.PROCESS_TIME = 0

def ptposs(c):
    ls = []
    for a in range(0, int(2 / c)+1):
        for r in range(0, int(1 / c) + 1):
            if r!=0:
                ls.append([a * c * np.pi, r * c])
    return ls

        
lsv = ptposs(0.1)
'''
absi = [r*np.cos(a) for (a,r) in lsv]
ordo = [r*np.sin(a) for (a,r) in lsv]
plt.plot(absi,ordo,'.r')'''
plt.show()
lsd = set()



def pol2car(ls):
    lsc = []
    for (a,r) in ls:
        lsc.append([r*np.cos(a),r*np.sin(a)])
    return lsc

#ra = random()
li = [[0.0, 1.0],
 [0.8961, 1.0],
 [1.7952, 1.0],
 [2.6779, 1.0],
 [3.6053, 1.0],
 [4.4870, 1.0],
 [5.3871, 1.0],
 [0.4636, 0.5590]]
lir=[[randrange(0,int(2 / 0.1)+1)*0.1*np.pi,uniform(0.5,1)] for i in range(8)]
class Noeud:
    def __init__(self,combinaison=lir,fils=None,parent=None):
        self.combinaison = combinaison
        self.fils = [] if fils is None else fils
        self.maxi , _ =  ft.optimal_tree(0,[[0,0]]+pol2car(self.combinaison),ft.dist_L2)
        self.parent = parent

    def suivant(self,ind):
        for i in lsv:
            if 1==1 and \
               (i not in [[0,0]]+self.combinaison) and \
               ((abs(i[1]-self.combinaison[ind][1])<0.1*np.pi+1e-5 and (i[0]-self.combinaison[ind][0])==0)) or (abs(i[1]-self.combinaison[ind][1])==0 and abs(i[0]-self.combinaison[ind][0]) <0.1*np.pi+1e-5):

                config = frozenset(map(tuple, self.combinaison[:ind] + [i] + self.combinaison[(ind+1):]))
                if config not in lsd:
                    ls = [[0,0]] + pol2car(self.combinaison[:ind]+[i]+self.combinaison[(ind+1):])
                    x, _ = ft.optimal_tree(0, ls, ft.dist_L2)
                    lsd.add(config)
                    if self.maxi < x:
                        self.fils.append(Noeud(self.combinaison[:ind]+[i]+self.combinaison[(ind+1):], parent=self))

    def cons(self):
        for ind in range(len(self.combinaison)):
            self.suivant(ind)
            
            
from time import time 
import heapq
import gc
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
    while pile and time()-t0<5000:
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
    try:
        while not file.est_vide() and time() - t0 < 10000:
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
            
    except KeyboardInterrupt:
        # Si on interrompt l'exécution avec Ctrl+C, cette partie sera exécutée
        print("Exécution interrompue manuellement.")
        return meilleur
    return meilleur

rac = Noeud()
x, T = ft.optimal_tree(0, [[0,0]]+pol2car(rac.combinaison), ft.dist_L2)
ft.draw_all('optimal', x, T)

rep = trouver(rac, 3.48)
x, T = ft.optimal_tree(0, [[0,0]]+pol2car(rep.combinaison), ft.dist_L2)
ft.draw_all('optimal', x, T)

noeudac = rep
evolution = []
while noeudac:
    evolution.append([[0,0]] + pol2car(noeudac.combinaison))
    noeudac = noeudac.parent
evolution.reverse()

image_files = []
for idx, evo in enumerate(evolution):
    x, T = ft.optimal_tree(0, evo, ft.dist_L2)
    plt.clf()
    ft.draw_all('optimal', x, T, save=str(idx)+'.png')
    plt.close('all')
    gc.collect()
