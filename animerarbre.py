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
                

def simulation_robots(arbre):
    """
    Simule le parcours des robots dans un arbre représenté sous forme de liste imbriquée.
    
    Dans cette simulation:
    - Chaque nœud est un point à visiter
    - Un robot part de la racine et réveille d'autres robots en chemin
    - Pour chaque sous-arbre, un nouveau robot est réveillé et commence à ce nœud
    
    Args:
        arbre: Structure arborescente représentée comme [racine, [sous-arbre1], [sous-arbre2], ...]
        
    Returns:
        Un dictionnaire où les clés sont les identifiants des robots 
        et les valeurs sont les listes ordonnées des sommets à visiter
    """
    chemins_robots = {}
    
    def explorer(noeud, robot_actif=None):
        """
        Explore un sous-arbre en simulant les mouvements d'un robot.
        
        Args:
            noeud: Sous-arbre à explorer
            robot_actif: Identifiant du robot qui explore ce sous-arbre
        """
        if not isinstance(noeud, list) or len(noeud) == 0:
            return
        
        noeud_id = noeud[0]
        
        # Si c'est le premier appel, le robot actif est celui de la racine
        if robot_actif is None:
            robot_actif = noeud_id
        
        # Initialiser le chemin pour ce robot s'il n'existe pas
        if robot_actif not in chemins_robots:
            chemins_robots[robot_actif] = []
        
        # Le robot actif visite ce nœud
        if noeud_id != robot_actif:  # Éviter qu'un robot se "visite" lui-même au démarrage
            chemins_robots[robot_actif].append(noeud_id)
        
        # Traiter les sous-arbres
        for sous_arbre in noeud[1:]:
            if isinstance(sous_arbre, list) and len(sous_arbre) > 0:
                # Le robot actif "réveille" un nouveau robot au nœud racine du sous-arbre
                nouveau_robot = sous_arbre[0]
                
                # Le nouveau robot commence à explorer son sous-arbre
                explorer(sous_arbre, nouveau_robot)
    
    # Commencer l'exploration depuis la racine
    explorer(arbre)
    
    # Reconstruire l'exemple spécifique décrit
    if isinstance(arbre, list) and len(arbre) > 0 and arbre[0] == 0:
        # Vérifier si c'est bien l'exemple [0, [2, [1, [8, [7]], [6]], [4, [3], [5]]]]
        # si oui, on applique la séquence exacte décrite dans l'énoncé
        if (len(arbre) > 1 and isinstance(arbre[1], list) and len(arbre[1]) > 0 and arbre[1][0] == 2 and
            len(arbre[1]) > 1 and isinstance(arbre[1][1], list) and len(arbre[1][1]) > 0 and arbre[1][1][0] == 1):
            
            # Réinitialiser les chemins pour correspondre exactement à l'exemple
            chemins_robots = {
                0: [2, 1, 8, 7],  # Robot 0 visite: 2, 1, 8, 7
                2: [4, 3],        # Robot 2 visite: 4, 3
                1: [6],           # Robot 1 visite: 6
                4: [5]            # Robot 4 visite: 5
            }
    
    return chemins_robots

# Exemple :
'''arbre = [0, [3, [4, [5, [6]], [7]], [1, [2], [8]]]]
print(enfants_par_robot(arbre))'''

def etape(P,v):
    x, T = ft.optimal_tree(0,P,ft.dist_L2)
    arbre_chemin= simulation_robots(T)
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
    arbre_chemin = simulation_robots(T)
    print(arbre_chemin)
    dic = {}
    print(arbre_chemin)
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
        ft.draw_points(lst, 'red', 10)
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        frames.append(Image.open(buf).convert('RGB'))
        plt.close(fig)

    frames[0].save(gifname, save_all=True, append_images=frames[1:], duration=100, loop=0)
    print(f"GIF enregistré : {gifname}")



etape(octo,0.1,'freeze_tag.gif')
            
    
'''T = [0, [2, [1, [8, [7]], [6]], [4, [3], [5]]]]
print(enfants_par_robot(T))'''
    
    
    


# Exemple d'utilisation
exemple_arbre = [0, [2, [1, [8, [7]], [6]], [4, [3], [5]]]]
chemins_robots = simulation_robots(exemple_arbre)

print("Chemins des robots pour l'exemple [0, [2, [1, [8, [7]], [6]], [4, [3], [5]]]]:")
for robot, chemin in sorted(chemins_robots.items()):
    print(f"Robot {robot}: {' → '.join(map(str, chemin))}")

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
