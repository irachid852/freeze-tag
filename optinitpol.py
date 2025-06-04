# -*- coding: utf-8 -*-
"""
Created on Tue May 13 15:43:40 2025

@author: ismail rachid
"""
import freeze_tag_mini as ft
import numpy as np
import matplotlib.pyplot as plt
from random import randrange,uniform
import algo2f as af
ft.SEED = -1
ft.PROCESS_TIME = 0

def ptposs(c):
    ls = []
    for a in range(0, int(2 / c)+1):
        for r in range(0, int(1 / c) + 1):
            if r!=0:
                ls.append([a * c * np.pi, r * c])
    return ls

        
lsv = ptposs(0.05)
lsv.remove([0,1])
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
lir=[[randrange(0,int(2 / 0.1)+1)*0.1*np.pi,uniform(0.5, 0.8)] for i in range(7)]
#lir=[[randrange(0,int(2 / 0.1)+1)*0.1*np.pi,0.45]]+[[randrange(0,int(2 / 0.1)+1)*0.1*np.pi,uniform(0.5, 0.8)] for i in range(6)]
lio= [[0.0, 1], [0.785, 1], [1.571, 1], [2.356, 1],
 [3.142, 1], [3.927, 1], [4.712, 1], [5.498, 1]]
class Noeud:
    def __init__(self,combinaison=lir,fils=None,parent=None):
        self.combinaison = combinaison
        self.fils = [] if fils is None else fils
        self.maxi , self.arbre =  ft.optimal_tree(0,[[0,0],[1,0]]+pol2car(self.combinaison),ft.dist_L2)
        #self.maxi , _ =  af.temps_rev(af.recherche(8,[[0,0],[0,1]]+self.combinaison))
        self.parent = parent

    def suivant(self,ind):
        for i in lsv:
            if 1==1 and \
               ([round(i[0],3),round(i[1],3)] not in [[round(j[0],3),round(j[1],3)] for j in [[0,0],[0,1]]+self.combinaison]) and \
               ((abs(i[1]-self.combinaison[ind][1])<0.1 * np.pi + 1e-5 and (i[0]-self.combinaison[ind][0])==0)) or (abs(i[1]-self.combinaison[ind][1])==0 and abs(i[0]-self.combinaison[ind][0]) <0.1 * np.pi + 1e-5):

                config = frozenset(map(tuple, self.combinaison[:ind] + [i] + self.combinaison[(ind+1):]))
                if config not in lsd:
                    #ls = [[0,0]] + pol2car(self.combinaison[:ind]+[i]+self.combinaison[(ind+1):])
                    ls = [[0,0],[0,1]] + self.combinaison[:ind]+[i]+self.combinaison[(ind+1):]
                    x, _ = ft.optimal_tree(0, ls, ft.dist_L2)
                    #x, _ = af.temps_rev(af.recherche(8,ls))
                    lsd.add(config)
                    if 1==1 and self.maxi < x:
                        self.fils.append(Noeud(self.combinaison[:ind]+[i]+self.combinaison[(ind+1):], parent=self))
    def cons(self):
        for ind in range(len(self.combinaison)):
            self.suivant(ind)
""" def suivant(self, ind):
        for i in lsv:
            if [round(i[0], 3), round(i[1], 3)] not in [[round(j[0], 3), round(j[1], 3)] for j in [[0, 0], [1, 0]] + self.combinaison]:
    
                a0, r0 = self.combinaison[ind]
                a1, r1 = i
                delta_a = abs(a1 - a0)
                delta_r = abs(r1 - r0)
    
                if (ind == 0 and delta_r < 1e-10 and delta_a < 0.1*np.pi + 1e-5) or \
                   (ind != 0 and ((delta_r < 0.1 + 1e-5 and delta_a < 1e-5) or (delta_a < 0.1*np.pi + 1e-5 and delta_r < 1e-5))):
    
                    config = frozenset(map(tuple, self.combinaison[:ind] + [i] + self.combinaison[ind+1:]))
                    if config not in lsd:
                        ls = [[0, 0]] + pol2car(self.combinaison[:ind] + [i] + self.combinaison[ind+1:])
                        x, _ = ft.optimal_tree(0, ls, ft.dist_L2)
                        lsd.add(config)
                        if self.maxi < x:
                            self.fils.append(Noeud(self.combinaison[:ind] + [i] + self.combinaison[ind+1:], parent=self))

    def suivant(self, ind):
        points_existants = set((round(p[0], 3), round(p[1], 3)) for p in [[0, 0], [1, 0]] + self.combinaison)
    
        for i in lsv:
            p_arrondi = (round(i[0], 3), round(i[1], 3))
            if p_arrondi not in points_existants:
                a0, r0 = self.combinaison[ind]
                a1, r1 = i
                delta_a = abs(a1 - a0)
                delta_r = abs(r1 - r0)
    
                if (ind == 0 and delta_r < 1e-10 and delta_a < 0.1 * np.pi + 1e-5) or \
                   (ind != 0 and ((delta_r < 0.1 + 1e-5 and delta_a < 1e-5) or 
                                 (delta_a < 0.1 * np.pi + 1e-5 and delta_r < 1e-5))):
    
                    nouvelle_combinaison = self.combinaison[:ind] + [i] + self.combinaison[ind + 1:]
                    config = frozenset(map(tuple, nouvelle_combinaison))
    
                    if config not in lsd:
                        ls = [[0, 0],[1,0]] + pol2car(nouvelle_combinaison)
                        x, _ = ft.optimal_tree(0, ls, ft.dist_L2)
                        lsd.add(config)
    
                        if 1==1:#self.maxi < x:
                            self.fils.append(Noeud(nouvelle_combinaison, parent=self))
                            points_existants.add(p_arrondi)  # On enregistre le point comme utilisé
"""
"""
    #on bouge tout selon l angle et le rayon sauf un ou on bouge pas le rayon
    def suivant(self, ind):
        points_existants = set((round(p[0], 3), round(p[1], 3)) for p in [[0, 0], [1, 0]] + self.combinaison)
    
        for i in lsv:
            p_arrondi = (round(i[0], 3), round(i[1], 3))
            if p_arrondi not in points_existants:
                a0, r0 = self.combinaison[ind]
                a1, r1 = i
                delta_a = abs(a1 - a0)
                delta_r = abs(r1 - r0)
    
                if (ind == 20 and delta_r < 1e-10 and delta_a < 0.1 * np.pi + 1e-5) or \
                   (ind != 20 and ((delta_r < 0.1 + 1e-5 and delta_a < 1e-5) or 
                                 (delta_a < 0.1 * np.pi + 1e-5 and delta_r < 1e-5))):
    
                    nouvelle_combinaison = self.combinaison[:ind] + [i] + self.combinaison[ind + 1:]
                    config = frozenset(map(tuple, nouvelle_combinaison))
    
                    if config not in lsd:
                        ls = [[0, 0]] + pol2car(nouvelle_combinaison)
                        x, _ = ft.optimal_tree(0, ls, ft.dist_L2)
                        lsd.add(config)
    
                        if 1==1:#self.maxi < x:
                            self.fils.append(Noeud(nouvelle_combinaison, parent=self))
                            points_existants.add(p_arrondi)  # On enregistre le point comme utilisé

"""
    
            
            
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
            print('vfin')
            return meilleur
        elif option == 'local' and abs(meilleur.maxi - courant.maxi) < seuil:
            return meilleur

        pile.extend(courant.fils)
        courant.fils = []
        gc.collect()
    print("fin")
    return meilleur

def trouver(racine, seuil, option='seuil'):
    racine.cons()
    file = FilePrioriteMax()
    for f in racine.fils:
        file.ajouter(f)
    
    meilleur = racine
    t0 = time()
    try:
        while not file.est_vide() and time() - t0 < 100000:
            courant = file.retirer()
            courant.cons()
            #print(courant.maxi)
    
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
    print('fin')
    return meilleur

rac = Noeud()
x, T = ft.optimal_tree(0, [[0,0],[1,0]]+pol2car(rac.combinaison), ft.dist_L2)
ft.draw_all('optimal', x, T,vis=True)

rep = trouver(rac, 50)
x, T = ft.optimal_tree(0, [[0,0],[1,0]]+pol2car(rep.combinaison), ft.dist_L2)
ft.draw_all('optimal', x, T,vis=True)
dics = {}
noeudac = rep
evolution = []
i = 1
ev= []
while noeudac:
    ev.append([noeudac.maxi,noeudac.arbre])
    #[f"config {i}"] = [noeudac.maxi,noeudac.arbre]
    i+=1
    evolution.append([[0,0],[1,0]] + pol2car(noeudac.combinaison))
    noeudac = noeudac.parent
ev.reverse()
#for h, elt in dics.items():
for j in range(len(ev)):
    dics[f"config {j+1}"] = ev[j]
evolution.reverse()
'''
image_files = []
for idx, evo in enumerate(evolution):
    x, T = ft.optimal_tree(0, evo, ft.dist_L2)
    plt.clf()
    ft.draw_all('optimal', x, T, save=str(idx)+'.png')
    plt.close('all')
    gc.collect()
'''

from PIL import Image
import matplotlib.pyplot as plt
import io
frames = []
option = int(input("0 pour representation avec les abres et 1 pour les representations sans les arbres : "))
for idx, evo in enumerate(evolution):
    
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    #plt.axis('equal')
    x,T = dics[f"config {idx+1}"]
    ft.draw_disc()
    if option == 0: ft.draw_all('optimal', x,T,P = evo )
    else : ft.draw_points(evo,'green')
    #plt.savefig(f'{idx}.png')
    buf = io.BytesIO()
    
    plt.savefig(buf, format='png')
    buf.seek(0)
    frames.append(Image.open(buf).convert('RGB'))
    if idx == 0:
        Image.open(buf).convert('RGB').save("etape_initiale.png")
    if idx == len(evolution) - 1:
        Image.open(buf).convert('RGB').save("etape_finale.png")
    plt.close(fig)
    #plt.show()
    #plt.close('all')
    gc.collect()
    
frames[0].save("pire.gif", save_all=True, append_images=frames[1:], duration=1000, loop=0)
while True:
    try:
        reponse = input(f"quelle image souhaitez vous visualiser et avec ou sans arbre sur {len(evolution)} images (ex : 12/0 donne image 12 avec arbres, on commence a 1) : ")
        option = reponse.split("/")
        x,T = dics[f"config {int(option[0])}"]
        if int(option[1]) == 0: ft.draw_all('optimal', x,T,P = evolution[int(option[0])-1],vis=True )
        else : 
            ft.draw_points(evolution[int(option[0])-1],'green')
            plt.show()
    except KeyboardInterrupt:
        # Si on interrompt l'exécution avec Ctrl+C, cette partie sera exécutée
        print("Exécution interrompue manuellement.")
        break
