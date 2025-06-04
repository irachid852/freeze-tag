# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 10:18:22 2025

@author: racni
"""

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

def initial(n,lsc=None):
    if lsc == None:
        
        lsr = [Robot([0,0],True)] + [Robot([rd.uniform(0,2*np.pi), np.sqrt(rd.uniform(0.2,1))], False) for i in range(n)] 
        print("bizarre")
    else : lsr = [Robot([0,0],True)]+ [Robot(lsc[i][:], False) for i in range(1,1+n)] 
    distances = [(np.linalg.norm(np.array(pol2car(e.coord))), e) for e in lsr if e.coord != [0,0]]
    distances.sort(key=lambda x: x[0])
    
    if len(distances) == 0:
        minval = None
        minval2 = None
    elif len(distances) == 1:
        minval = distances[0][1]
        minval2 = None
    else:
        minval = distances[0][1]
        minval2 = distances[1][1]
    
    """
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
            elif dist == mini:
                mini2 = dist
                minval2 = e"""
    dic = {'couronne 1': [], 'couronne 2': [], 'couronne 3': []}
    #aa.ft.draw_disc()
    dicc = {'couronne 1': "ro", 'couronne 2': "bo", 'couronne 3': "go"}
    #print(minval2)
    #print(minval2.coord)
    #print([r.coord for r in lsr])
    for e in lsr:
        t, r = e.coord
        t = t % (2*np.pi)  # Normaliser l'angle entre 0 et 2π
        
        if 0 <= (t-minval2.coord[0])%(2*np.pi) < 2*np.pi/3:
            dic['couronne 1'].append(e)
            e.groupe = 'couronne 1'
        elif 2*np.pi/3 <= (t-minval2.coord[0])%(2*np.pi) < 4*np.pi/3:
            dic['couronne 2'].append(e)
            e.groupe = 'couronne 2'
        else:
            dic['couronne 3'].append(e)
            e.groupe = 'couronne 3'

    """
    for r in lsr:
        x,y = pol2car(r.coord)
        if r not in [lsr[0],minval,minval2]:
            
            plt.plot(x,y,dicc[r.groupe])
        else:
            plt.plot(x,y,"mo")
    plt.show()"""
    l = minval2.coord[0]
    for r in lsr:
        r.coord[0] = (r.coord[0]-l)%(2*np.pi)
    lsr[0].aller(minval)
    minval.reveil = True
    #lsr[0].aller(r1)
    #minval.aller(r1)
    minval.aller(minval2)
    minval2.reveil = True
    


    """
    aa.ft.draw_disc()
    dicc = {'couronne 1': "ro", 'couronne 2': "bo", 'couronne 3': "go"}
    for r in lsr:
        x,y = pol2car(r.coord)
        if r not in [lsr[0],minval,minval2]:
            
            plt.plot(x,y,dicc[r.groupe])
        else:
            plt.plot(x,y,"mo")
    plt.show()"""
    return lsr, minval, minval2, dic

def recherche_couronne(lsd, couronne, nom,minval,lsr,minval2):
    angleinf = ((int(nom[-1])-1)*2*np.pi/3)%(2*np.pi)
    anglesup = ((int(nom[-1]))*2*np.pi/3)%(2*np.pi)
    mind = np.Inf
    valmind = None
    for eld in lsd:
        if eld== lsr[0] :
            t = minval.coord[0]
        elif eld == minval: t = minval2.coord[0]
        else: t = eld.coord[0]
        k = [abs(t-angleinf),abs(t-anglesup)]
        #print(k)
        m = min(k)
        vk = [k.index(m),eld]
        if m<mind:
            mind = m
            valmind = vk
    #print(valmind[0])
    d = valmind[1]
    file = File()
    file.enfiler(d)
    # Initialisation correcte des bornes (exemple avec 0 et 1)
    d.borne_rayon_inf = minval2.coord[1]
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
            t = e.coord[0]  # Supposons que coord[1] est le rayon
            diff = abs(t-(angleinf+valmind[0]*2*np.pi/3))
            if diff < min_diff and e.coord[1]!=0 and e.coord not in [minval.coord,minval2.coord] and e!=r1:
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
    return d

def recherche(n,lsc=None):
    
    lsr, minval, minval2, dic = initial(n,lsc)
    lsd = [minval,minval2,lsr[0]]
    #print("liste ", [pol2car(r.coord) for r in lsd])
    v = recherche_couronne(lsd,dic['couronne 1'],'couronne 1',minval,lsr,minval2)
    #print("couronne 1",pol2car(v.coord))
    
    lsd.remove(v)
    k = recherche_couronne(lsd,dic['couronne 2'],'couronne 2',minval,lsr,minval2)
    #print("couronne 2",pol2car(k.coord))
    lsd.remove(k)
    recherche_couronne(lsd,dic['couronne 3'],'couronne 3',minval,lsr,minval2)
    #print("couronne 3",pol2car(recherche_couronne(lsd,dic['couronne 3'],'couronne 3',minval,lsr,minval2).coord))
    lst = [i for i in lsr if i.reveil == True and i.groupe == 'couronne 1']
    lsv = [i for i in lsr if i.groupe == 'couronne 1'] 
    #print(f"{len(lst)}/{len(lsv)}")
    return lsr


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
#etape(lsr,0.045,'freeze_tag.gif')

#ls = [[rd.uniform(0,2*np.pi), np.sqrt(rd.uniform(0.2,1))] for i in range(49)]
#lsr = recherche(50,ls)
def temps_rev(lsr):
    maxi = 0 
    maxval = None
    for r in lsr:
        dist = 0
        elt = np.array(pol2car(r.coord))
        for e in r.suivant:
            val = np.array(pol2car(e.coord))
            dist+=np.linalg.norm(val-elt)
            elt = val
        if dist>maxi:
            maxi = dist
            maxval = r
    return maxi, maxval

#k,f = temps_rev(lsr)
#print(k)























