# -*- coding: utf-8 -*-
"""
Created on Thu May 22 10:49:44 2025

@author: racni
"""

import random as rd
import numpy as np
import copy as cp
import matplotlib.pyplot as plt 
import animerarbre as aa
def pol2car(ls):
    
    return [np.cos(ls[0])*ls[1],np.sin(ls[0])*ls[1]]
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
    
class Robot:
    def __init__(self,coord,reveil):
        self.coord, self.reveil = coord, reveil
        self.groupe = None
        self.suivant = []
        self.borne_rayon_sup = None
        self.borne_rayon_inf = None
        
    def aller(self,coord):
        self.suivant.append(coord)
r1 = Robot([0,0],False)
lsr = [Robot([0,0],True)] + [Robot([rd.uniform(0,2*np.pi), np.sqrt(rd.uniform(0,1))], False) for i in range(500)] + [r1]

mini = np.Inf 
minval = None
mini2 = np.Inf
minval2 = None
for e in lsr:
    if e.coord !=[0,0]:
        dist = np.linalg.norm(np.array(pol2car(e.coord)))
        if dist<mini: 
            mini = dist; minval = e
        elif dist<mini2:
            mini2 = dist
            minval2 = e
    
 
lsr[0].aller(minval)
minval.reveil = True
lsr[0].aller(r1)
minval.aller(r1)
minval.aller(minval2)
minval2.reveil = True


dic = {'couronne 1': [], 'couronne 2': [], 'couronne 3': []}

for e in lsr:
    t, r = e.coord
    t = t % (2*np.pi)  # Normaliser l'angle entre 0 et 2π
    
    if 0 <= t < 2*np.pi/3:
        dic['couronne 1'].append(e)
        e.groupe = 'couronne 1'
    elif 2*np.pi/3 <= t < 4*np.pi/3:
        dic['couronne 2'].append(e)
        e.groupe = 'couronne 2'
    else:
        dic['couronne 3'].append(e)
        e.groupe = 'couronne 3'
aa.ft.draw_disc()
dicc = {'couronne 1': "ro", 'couronne 2': "bo", 'couronne 3': "go"}
for r in lsr:
    x,y = pol2car(r.coord)
    if r not in [lsr[0],minval,minval2]:
        
        plt.plot(x,y,dicc[r.groupe])
    else:
        plt.plot(x,y,"mo")
plt.show()
def recherche_couronne(d,couronne,nom):
    file = File()
    file.enfiler(d)
    d.borne_rayon_inf = minval2.coord[1] 
    d.borne_rayon_sup = 1
    if d in couronne : couronne.remove(d)
    
    while not file.est_vide() and not couronne == []:
        elt = file.premier()
        minc =np.Inf
        mincour = None
        print(elt.coord,elt.borne_rayon_inf,elt.borne_rayon_sup)
        for e in couronne:
            t, r = e.coord 
            if elt.suivant != []: dist = abs(t - elt.suivant[-1].coord[0])
            else : dist = abs(t - elt.coord[0])
            if dist < minc and elt.borne_rayon_inf<r<elt.borne_rayon_sup:
                minc = t
                mincour = e
            
        if mincour == None:
            file.defiler()
        else :
            elt.aller(mincour)
            mincour.reveil = True
            mincour.borne_rayon_sup = (elt.borne_rayon_inf + elt.borne_rayon_sup) / 2
            mincour.borne_rayon_inf = elt.borne_rayon_inf
            elt.borne_rayon_inf = mincour.borne_rayon_sup
            if mincour in couronne :couronne.remove(mincour)
            file.enfiler(mincour)
        #print(couronne)

def recherche_couronne(d,couronne,nom):
    file = File()
    file.enfiler(d)
    d.borne_rayon_inf = minval2.coord[1] 
    d.borne_rayon_sup = 1
    if d in couronne : couronne.remove(d)
    
    while not file.est_vide() and not couronne == []:
        elt = file.premier()
        minc = np.Inf
        mincour = None
        print(elt.coord,elt.borne_rayon_inf,elt.borne_rayon_sup)
        
        for e in couronne:
            t, r = e.coord 
            if elt.suivant != []: 
                dist = abs(t - elt.suivant[-1].coord[0])
            else : 
                dist = abs(t - elt.coord[0])
            
            if dist < minc and elt.borne_rayon_inf < r < elt.borne_rayon_sup:
                minc = dist  # ← Corriger ici : assigner 'dist' pas 't'
                mincour = e
            
        if mincour == None:
            file.defiler()  # ← Aussi corriger ici : defiler() pas premier()
        else :
            file.defiler()  # ← Retirer l'élément traité
            elt.aller(mincour)
            mincour.reveil = True
            mincour.borne_rayon_sup = (elt.borne_rayon_inf + elt.borne_rayon_sup) / 2
            mincour.borne_rayon_inf = elt.borne_rayon_inf
            elt.borne_rayon_inf = mincour.borne_rayon_sup
            if mincour in couronne : couronne.remove(mincour)
            file.enfiler(mincour)
            
def recherche_couronne(d, couronne, nom):
    file = File()
    file.enfiler(d)
    # Initialisation correcte des bornes (exemple avec 0 et 1)
    d.borne_rayon_inf = 0
    d.borne_rayon_sup = 1
    if d in couronne:
        couronne.remove(d)
    
    while not file.est_vide() and couronne:
        elt = file.defiler()  # Supposons que "defiler" retire et retourne l'élément
        # Trouver le nœud dans couronne dont le rayon est au milieu de l'intervalle
        milieu = (elt.borne_rayon_inf + elt.borne_rayon_sup) / 2
        mincour = None
        min_diff = np.Inf
        # Recherche du nœud le plus proche du milieu
        for e in couronne:
            r = e.coord[1]  # Supposons que coord[1] est le rayon
            diff = abs(r - milieu)
            if diff < min_diff and r!=0 and e.coord not in [minval.coord,minval2.coord] :
                min_diff = diff
                mincour = e
        if mincour:
            elt.aller(mincour)
            mincour.reveil = True
            # Mise à jour des bornes pour les enfants
            mincour.borne_rayon_inf = elt.borne_rayon_inf
            mincour.borne_rayon_sup = milieu
            elt.borne_rayon_inf = milieu  # Ajuster la borne du parent
            couronne.remove(mincour)
            file.enfiler(mincour)
            file.enfiler(elt)  # Conserver le parent pour l'autre moitié
'''  
def recherche_couronne(d, couronne, nom):
    file = File()
    file.enfiler(d)
    d.borne_rayon_inf = minval2.coord[1] 
    d.borne_rayon_sup = 1
    if d in couronne:
        couronne.remove(d)

    while not file.est_vide() and couronne:
        elt = file.defiler()
        rmin, rmax = elt.borne_rayon_inf, elt.borne_rayon_sup
        r_centre = (rmin + rmax) / 2

        best = None
        min_dist = np.inf
        for e in couronne:
            t, r = e.coord
            if rmin < r < rmax:
                dist = abs(r - r_centre)
                if dist < min_dist:
                    min_dist = dist
                    best = e

        if best:
            elt.aller(best)
            best.reveil = True
            best.borne_rayon_inf = rmin
            best.borne_rayon_sup = r_centre
            elt.borne_rayon_inf = r_centre  # maj pour la suite de l'arbre
            file.enfiler(best)
            couronne.remove(best)
def recherche_couronne(d, couronne, nom):
    file = File()
    file.enfiler(d)
    d.borne_rayon_inf = minval2.coord[1] 
    d.borne_rayon_sup = 1

    if d in couronne:
        couronne.remove(d)

    while not file.est_vide() and couronne:
        elt = file.defiler()
        r_inf = elt.borne_rayon_inf
        r_sup = elt.borne_rayon_sup
        r_milieu = (r_inf + r_sup) / 2

        # chercher le robot dont le rayon est le plus proche de r_milieu dans l’intervalle
        min_rayon_diff = np.inf
        candidat = None
        for e in couronne:
            _, r = e.coord
            if r_inf < r < r_sup:
                diff = abs(r - r_milieu)
                if diff < min_rayon_diff:
                    min_rayon_diff = diff
                    candidat = e

        if candidat is not None:
            elt.aller(candidat)
            candidat.reveil = True
            candidat.borne_rayon_inf = r_inf
            candidat.borne_rayon_sup = r_milieu
            elt.borne_rayon_inf = r_milieu
            file.enfiler(candidat)
            couronne.remove(candidat)
'''
'''

recherche_couronne_binaire(minval, dic['couronne 1'], 0, 2*np.pi/3, minval2.coord[1], 1)
recherche_couronne_binaire(minval2, dic['couronne 2'],  2*np.pi/3, 4*np.pi/3, minval2.coord[1], 1)
recherche_couronne_binaire(lsr[0], dic['couronne 3'],  4*np.pi/3, 6*np.pi/3, minval2.coord[1], 1)
'''

recherche_couronne(minval,dic['couronne 1'],'couronne 1')
recherche_couronne(minval2,dic['couronne 2'],'couronne 2')
recherche_couronne(lsr[0],dic['couronne 3'],'couronne 3')

lst = [i for i in lsr if i.reveil == True and i.groupe == 'couronne 1']
lsv = [i for i in lsr if i.groupe == 'couronne 1'] 
print(f"{len(lst)}/{len(lsv)}")



from PIL import Image
import matplotlib.pyplot as plt
import io

def pol2car(ls):
    
    return [np.cos(ls[0])*ls[1],np.sin(ls[0])*ls[1]]

def etape(lsr, v, gifname='freeze_tag.gif'):
    
   
    dic = {}
    
    ls = [False] * len(lsr)
    for i, el in enumerate(lsr):
        dic[f'Robot {i}'] = aa.Noeud(pol2car(el.coord), None, False)
    for i in range(len(lsr)):
        #if not i in arbre_chemin.keys(): arbre_chemin[i] = []
        dic[f'Robot {i}'].suivant = [dic[f'Robot {lsr.index(j)}'] for j in lsr[i].suivant ]
    dic['Robot 0'].reveil = True
    
    frames = []
    
    while ls != [True] * len(lsr):
        for i in range(len(lsr)):
            ls[i] = dic[f'Robot {i}'].aller(v)
        lst = [dic[f'Robot {i}'].coord for i in range(len(lsr))]

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        aa.ft.draw_disc()
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


#la = [[0,0]]+pol2car([[ uniform(0,2*np.pi),uniform(0.5,1)] for i in range(12)])
etape(lsr,0.05,'freeze_tag.gif')




























