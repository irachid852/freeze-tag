# -*- coding: utf-8 -*-
"""
Created on Wed May 21 13:57:09 2025

@author: racni
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 19 09:46:50 2025

@author: ismail rachid
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
        self.borne_angle_sup = None
        self.borne_angle_inf = None
        self.borne_hauteur_sup = None
        self.borne_hauteur_inf = None
        self.bai = None
        self.bas = None
    def aller(self,coord):
        self.suivant.append(coord)
lsr = [Robot([0,0],True)] + [Robot([rd.uniform(0,2*np.pi),rd.uniform(0,1)],False) for i in range(1000)]
for ep in range(0,int(np.pi*100)+1):
    eps = ep/100
    lsf = []
    for e in lsr:
        r = e.coord[1]
        t = e.coord[0]
        if 0<=t<=eps:
            lsf.append(e)
    if len(lsf)>=10:
        break
file = File()
lst=lsf[:]
robot10 = Robot([0,1],False)
lsr.append(robot10)
file.enfiler(lsr[0])
lsr[0].bas = eps
lsr[0].bai = 0
if lsr[0] in lsf :lsf.remove(lsr[0])

while lsf !=[]:
    elt = file.premier()
    mini = np.Infinity
    minval = None
    for e in lsf:
        dist = e.coord[1]*np.cos(e.coord[0])
        if elt.bai <e.coord[0]< elt.bas and dist <mini:
            mini = dist 
            minval = e
    if minval == None:
        file.defiler()
        lst.append(elt)
    else:
        minval.reveil = True
        elt.aller(minval)
        minval.bas = (elt.bas+elt.bai)/2 
        minval.bai = elt.bai
        elt.bai = minval.bas
        file.enfiler(minval)
        lsf.remove(minval)
for r in lst:
    r.aller(robot10)
        
dic = {f"cone {i}": [] for i in range(3,len(lst)-3)}
dic["dome 1"] = []; dic["dome 2"] = []
pas = np.pi/len(lst)
angles = [pas*i for i in range(len(lst)+1)]
for e in lsr:
    if e not in lst:
        r = e.coord[1]
        t = e.coord[0]
        
        x, y = r * np.cos(t), r * np.sin(t)
        vec = np.array([x - 1, y])
        norm_vec = np.linalg.norm(vec)
        
        if norm_vec == 0:
            continue  # éviter division par zéro
        
        val = y / norm_vec
        val = np.clip(val, -1, 1)
        angle = np.arccos(val)
        for i in range(3,len(lst)-3):
            
    
            #angle = np.arccos(y/(r**2+1-2*x))  # angle entre [1,0] et le robot
            #(np.arctan(r*np.sin(t) / (r*np.cos(t)-1) ))
           
            #if angle<=0: angle+=2*np.pi
            #print(angle*360/(2*np.pi))
            if pas*(i+1)>=angle>=pas*i: 
                dic[f"cone {i}"].append(e)
                e.groupe = f"cone {i}"
                break
        if e.groupe == None :
            if angle<pas*3:
                dic["dome 1"].append(e)
                e.groupe = "dome 1"
            else:
                dic["dome 2"].append(e)
                e.groupe = "dome 2"


dicc = {f"cone {i}": [] for i in range(3,len(lst)-3)}
for c in dic.keys():
    dicc[c] = "#"+''.join([rd.choice('0123456789ABCDEF') for j in range(6)])

for r in lsr:
    if r.groupe !=None :
        plt.plot(r.coord[1]*np.cos(r.coord[0]),r.coord[1]*np.sin(r.coord[0]),'.',color =dicc[r.groupe] )
    else:
        plt.plot(r.coord[1]*np.cos(r.coord[0]),r.coord[1]*np.sin(r.coord[0]),'ob')
plt.show()

def recherche_cone(d: Robot, cone: list, nom: str):
    file = File()
    file.enfiler(d)
    d.borne_angle_inf = pas * int(nom.split(" ")[-1])
    d.borne_angle_sup = pas * (int(nom.split(" ")[-1]) + 1)

    if d in cone:
        cone.remove(d)

    while len(cone) > 0 and not file.est_vide():
        elt = file.premier()
        if elt.suivant != [] : rp, tp = elt.suivant[-1].coord[1], elt.suivant[-1].coord[0]
        else:rp, tp = elt.coord[1], elt.coord[0]
        xp, yp = rp * np.cos(tp), rp * np.sin(tp)

        mini = np.Infinity
        minval = None

        for e in cone:
            r, t = e.coord[1], e.coord[0]
            x, y = r * np.cos(t), r * np.sin(t)
            vec = np.array([x - 1, y])
            norm_vec = np.linalg.norm(vec)
            if norm_vec == 0:
                continue

            angle = np.arccos(np.clip(y / norm_vec, -1, 1))
            if elt.borne_angle_inf <= angle <= elt.borne_angle_sup:
                dist =np.linalg.norm([x - xp, y - yp]) #abs(rp**2+1-2*np.cos(tp) - r**2-1+2*np.cos(t)) 
                if dist < mini:
                    mini = dist
                    minval = e
        #print(mini,elt.borne_angle_inf,elt.borne_angle_sup)
        if minval is None:
            file.defiler()
        else:
            elt.aller(minval)
            minval.reveil = True
            minval.borne_angle_sup = (elt.borne_angle_inf + elt.borne_angle_sup) / 2
            minval.borne_angle_inf = elt.borne_angle_inf
            elt.borne_angle_inf = minval.borne_angle_sup

            file.enfiler(minval)
            cone.remove(minval)

#recherche_cone(lst[0],dic["cone 3"],"cone 3")
    
import numpy as np

def trouver_max_hauteur_dome(lsr, alpha):
    max_hauteur = -np.inf
    meilleur_robot = None

    cos_a, sin_a = np.cos(alpha), np.sin(alpha)

    for r in lsr:
        theta, rayon = r.coord
        x, y = rayon * np.cos(theta), rayon * np.sin(theta)
        hauteur = (x - 1) * cos_a + y * sin_a  # projection sur vecteur direction dôme

        if hauteur > max_hauteur:
            max_hauteur = hauteur
            meilleur_robot = r

    return meilleur_robot, max_hauteur




def trouver_max_hauteur_dome_inf(lsr, alpha):
    max_hauteur = -np.inf
    meilleur_robot = None

    cos_a, sin_a = np.cos(alpha), np.sin(alpha)

    for r in lsr:
        theta, rayon = r.coord
        x, y = rayon * np.cos(theta), rayon * np.sin(theta)
        hauteur = -( (x - 1) * cos_a + y * sin_a )  # projection opposée

        if hauteur > max_hauteur:
            max_hauteur = hauteur
            meilleur_robot = r

    return meilleur_robot, max_hauteur

alpha = 3*pas  # exemple : dôme incliné de 30°
robot, h = trouver_max_hauteur_dome(lsr, alpha)
print("Coordonnées :", robot.coord)
print("Hauteur max :", h)

alpha = (len(lst)-3)*pas  # même angle d’inclinaison
robot2, h2 = trouver_max_hauteur_dome_inf(lsr, alpha)
print("Robot plus bas du dôme inférieur :", robot2.coord)
print("Hauteur projetée maximale :", h2)


def recherche_dome(d, dome,nom):
    file = File()
    file.enfiler(d)
    
    if int(nom[-1]) == 1:
        alpha = 3*pas  # exemple : dôme incliné de 30°*
        cos_a, sin_a = np.cos(alpha), np.sin(alpha)
        robot, h = trouver_max_hauteur_dome(lsr, alpha)
        d.borne_hauteur_inf = 0
        d.borne_hauteur_sup = h
    else:
        alpha = (len(lst)-3)*pas  # exemple : dôme incliné de 30°
        cos_a, sin_a = np.cos(alpha), np.sin(alpha)
        robot, h = trouver_max_hauteur_dome_inf(lsr, alpha)
        d.borne_hauteur_inf = 0
        d.borne_hauteur_sup = h
    if d in dome:
        dome.remove(d)

    while len(dome) > 0 and not file.est_vide():
        elt = file.premier()
        if elt.suivant != [] : rp, tp = elt.suivant[-1].coord[1], elt.suivant[-1].coord[0]
        else:rp, tp = elt.coord[1], elt.coord[0]
        if int(nom[-1]) == 1:
            xp, yp = rp * np.cos(tp), rp * np.sin(tp)
            hauteur = (xp - 1) * np.cos(alpha) + yp * np.sin(alpha)
        else: 
            xp, yp = rp * np.cos(tp), rp * np.sin(tp)
            hauteur = (xp - 1) * np.cos(alpha) - yp * np.sin(alpha)

        
        mini = np.Infinity
        minval = None
        
        
        for e in dome:
            r, t = e.coord[1], e.coord[0]
            x, y = r * np.cos(t), r * np.sin(t)
            if int(nom[-1]) == 1:
                theta, rayon = e.coord
                x, y = rayon * np.cos(theta), rayon * np.sin(theta)
                hauteur = (x - 1) * cos_a + y * sin_a  # projection sur vecteur direction dôme
                
                #hauteur = (x - 1) * np.cos(alpha) + y * np.sin(alpha)
                proj = (x - 1) * cos_a + y * sin_a
                #proj = x * np.cos(alpha) + (y - 1) * np.sin(alpha)
                proj_ref = (xp - 1) * cos_a + yp * sin_a
                #proj_ref = xp * np.cos(alpha) + (yp - 1) * np.sin(alpha)
            else: 
                theta, rayon = e.coord
                x, y = rayon * np.cos(theta), rayon * np.sin(theta)
                hauteur = -( (x - 1) * cos_a + y * sin_a )  
                
                proj = -((x - 1) * cos_a + y * sin_a)
                
                proj_ref = -((xp - 1) * cos_a + yp * sin_a)
                
        
            if elt.borne_hauteur_inf <= hauteur <= elt.borne_hauteur_sup:
                dist = proj - proj_ref
                dist = np.linalg.norm([x - xp, y - yp])
                if 0 <= abs(dist) <= mini:  # réveil vers la droite uniquement
                    mini = dist
                    minval = e

        #print(mini,elt.borne_hauteur_inf,elt.borne_hauteur_sup)
        if minval is None:
            file.defiler()
        else:
            elt.aller(minval)
            minval.reveil = True
            minval.borne_hauteur_sup = (elt.borne_hauteur_inf + elt.borne_hauteur_sup) / 2
            minval.borne_hauteur_inf = elt.borne_hauteur_inf
            elt.borne_hauteur_inf = minval.borne_hauteur_sup

            file.enfiler(minval)
            dome.remove(minval)

ls = [i for i in lsr if i.groupe == 'cone 3']

"""recherche_dome(lst[0], dic['dome 1'], 'dome 1')
recherche_dome(lst[0], dic['dome 2'], 'dome 2')"""

ls1 = [i for i in lsr if i.groupe == 'cone 3']
ls2 = [i for i in lsr if i.groupe == 'cone 3']
coordsdejavu = []
for i, el in enumerate(dic.keys()):
    if el[0]=='c':
        print(el,lst[i].coord, i)
        recherche_cone(lst[i], dic[el], el)
        coordsdejavu.append(lst[i].coord)
    elif el[0] == 'd':
        print(el,lst[i].coord, i)
        recherche_dome(lst[i], dic[el], el)
        coordsdejavu.append(lst[i].coord)
        

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
lsk = [i for i in lsr if i.reveil == True]
plt.plot(np.cos(robot.coord[0])*robot.coord[1],np.sin(robot.coord[0])*robot.coord[1],'ro')
plt.plot(np.cos(robot2.coord[0])*robot2.coord[1],np.sin(robot2.coord[0])*robot2.coord[1],'mo')

plt.show()
    
    
    
    
    
    
    
    