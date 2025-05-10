import matplotlib.pyplot as plt
from random import randrange
import freeze_tag_mini as ft
from itertools import combinations

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



import math

somd = []
m = 0
t = [None]
pts=[]
lsv = ptposs(0.3)
lf = [[0,0]]
lsv.remove([0,0])

lfc = [list(i) for i in list(combinations(lsv,8))]
lsd =[]
class Noeud:
    def __init__(self,combinaison=lfc[0],fils=[]):
        self.combinaison = combinaison
        self.fils= fils
        self.maxi , _ =  ft.optimal_tree(0,[[0,0]]+self.combinaison,ft.dist_L2)
    def suivant(self,ind):
        deplacement= [i for i in lsv if (self.combinaison[ind][0] ==i[0]  or self.combinaison[ind][1] ==i[1] ) and not i in [[0,0]]+self.combinaison and not self.combinaison[:ind]+[i]+self.combinaison[(ind+1):] in lsd ]
        #print( [[0,0]]+self.combinaison)
        i=0
        #print(len(deplacement))
        #print(deplacement)
        for dep in deplacement:
            #print(dep)
            #print(f'{i}/{len(deplacement)}')
            i+=1
            
            ls = [[0,0]]+self.combinaison[:ind]+[dep]+self.combinaison[(ind+1):]
            #print(f'long : {len(ls)}')
            #print(ls)
            x, _ = ft.optimal_tree(0,ls,ft.dist_L2)
            if self.maxi<x:
                
                self.fils.append(Noeud(self.combinaison[:ind]+[dep]+self.combinaison[(ind+1):]))
            else:
                lsd.append(dep)
                
        #print('boom')
    def cons(self):
        for ind in range(len(self.combinaison)):
            self.suivant(ind)
            
    
            

def trouver(racine,seuil):
    racine.cons()
    pile = racine.fils
    meilleur = racine

    while pile:
        courant = pile.pop()
        courant.cons()
        print(courant.maxi)
        if courant.maxi > meilleur.maxi:
            print(f'nouveau max : {courant.maxi}')
            meilleur = courant
        elif  meilleur.maxi > seuil:
            return courant  # cas d'arrêt

        pile.extend(courant.fils)

    return meilleur
            
            
rac = Noeud()
x,T = ft.optimal_tree(0,[[0,0]]+rac.combinaison,ft.dist_L2)
ft.draw_all('optimal',x,T)
rep = trouver(rac,2.9)
x,T = ft.optimal_tree(0,[[0,0]]+rep.combinaison,ft.dist_L2)
ft.draw_all('optimal',x,T)

































