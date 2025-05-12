import matplotlib.pyplot as plt
from random import randrange
import freeze_tag_mini as ft
from itertools import combinations
import numpy as np
def ptposs(ep):
    ls =[]
    max_i = int(1 / ep)
    lx = [i*ep for i in range(-max_i, max_i + 1)]
    #print(lx)
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
ft.SEED = -1 # do not display SEED
ft.PROCESS_TIME = 0 # elapsed time temps to display in draw_all()


def homothetie(pts):
    l =np.sqrt( max((x**2 + y**2) for (x, y) in pts))
    return [(x/l, y/l) for (x, y) in pts]
import math

somd = []
m = 0
t = [None]
pts=[]
lsv = ptposs(0.1)
lf = [[0,0]]
lsv.remove([0,0])

#lfc = list(combinations(lsv,8))
lfc =[]
lsd =[]
li= [[0,0],(1.0, 0.0),(0.6234898018587336, 0.7818314824680298),(-0.22252093395631434, 0.9749279121818236),(-0.900968867902419, 0.43388373911755823),(-0.9009688679024191, -0.433883739117558),(-0.2225209339563146, -0.9749279121818236),(0.6234898018587334, -0.7818314824680299),(0.5,0.25)]


class Noeud:
    def __init__(self,combinaison=homothetie(li),fils=[],parent=None):
        self.combinaison = combinaison
        self.fils= fils
        self.maxi , _ =  ft.optimal_tree(0,[[0,0]]+self.combinaison,ft.dist_L2)
        self.parent = parent
    def suivant(self,ind):
        deplacement= [i for i in lsv if (self.combinaison[ind][0] ==i[0]  or self.combinaison[ind][1] ==i[1] ) and (not i in [[0,0]]+self.combinaison) and not self.combinaison[:ind]+[i]+self.combinaison[(ind+1):] in lsd ]
        #print( [[0,0]]+self.combinaison)
        i=0
        #print(len(deplacement))
        #print(deplacement)
        for dep in deplacement:
            #print(dep)
            #print(f'{i}/{len(deplacement)}')
            i+=1
            
            ls = [[0,0]]+homothetie(self.combinaison[:ind]+[dep]+self.combinaison[(ind+1):])
            #print(f'long : {len(ls)}')
            #print(ls)
            x, _ = ft.optimal_tree(0,ls,ft.dist_L2)
            if self.maxi<x:
                
                self.fils.append(Noeud(self.combinaison[:ind]+[dep]+self.combinaison[(ind+1):],parent=self))
            else:
                lsd.append(self.combinaison[:ind]+[dep]+self.combinaison[(ind+1):])
                
        #print('boom')
    def cons(self):
        for ind in range(len(self.combinaison)):
            self.suivant(ind)
            
    
            
from time import time 
def trouver(racine,seuil,option = 'seuil'):
    racine.cons()
    pile = racine.fils
    meilleur = racine
    t0 = time()
    while pile and time()-t0<300:
        courant = pile.pop()
        courant.cons()
        print(courant.maxi)
        if courant.maxi > meilleur.maxi:
            print(f'nouveau max : {courant.maxi}')
            meilleur = courant
        if option == 'seuil':
            if  meilleur.maxi > seuil:
                return meilleur  # cas d'arrêt
        elif option == 'local':
            if abs(meilleur.maxi-courant.maxi)<seuil:
                return meilleur
        pile.extend(courant.fils)

    return meilleur
            
     
rac = Noeud()
x,T = ft.optimal_tree(0,[[0,0]]+rac.combinaison,ft.dist_L2)
ft.draw_all('optimal',x,T)
rep = trouver(rac,3.4)
x,T = ft.optimal_tree(0,[[0,0]]+rep.combinaison,ft.dist_L2)
ft.draw_all('optimal',x,T)
noeudac = rep
evolution = []
while noeudac !=None:
    evolution.append([[0,0]]+homothetie(noeudac.combinaison))
    noeudac = noeudac.parent
evolution.reverse()
'''
for evo in evolution:
    x,T = ft.optimal_tree(0,evo,ft.dist_L2)
    ft.draw_all('optimal',x,T)'''
    
import matplotlib.pyplot as plt
image_files = []

for idx, evo in enumerate(evolution):
    # Crée un plot pour chaque étape
    x, T = ft.optimal_tree(0, evo, ft.dist_L2)
    
    # Sauvegarder l'image
    plt.clf()
    ft.draw_all('optimal', x, T,save=str(idx)+'.png')
   



























