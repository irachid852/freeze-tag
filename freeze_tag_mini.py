"""

            An Optimal Algorithm for the Freeze-Tag Problem in 2D

    
    The algorithm is based on a talk I gave in November 2022 at the ANR meeting in Chasseneuil-du-Poitou, France. Slides are available on my webpage, or directly here https://dept-info.labri.fr/~gavoille/gavoille_publis.html#talk:Gav22. It is a tiny part (so called "_mini") of a larger set of programs I've developped during Spring 2022, and that has been used in the DISC'24 paper and also in its full version (https://dept-info.labri.fr/~gavoille/gavoille_publis.html#BCGH24). It uses dynamic programing, and its time complexity is 3^n, up to some polynomials. It is significantly faster than the naive (n/e)^(n+O(1)) like algorithms, testing all possible wake-up trees on n sleeping robots. The number A(n) of trees to be tested is n! times the number of non-planar (unordered) rooted binary trees on n nodes (each non leaf having 1 or 2 children). The number A(n) does not seem to be known from OEIS, but is clearly no more than C(n) = binom(2n,n)/(n+1), the n-th Catalan number.

    It works with any distance function, not only with the L2 norm. Because of the dynamic programing nature of the program, we use numbers of global variables. This is a bad programming technique, sorry for that, but it speeds up a lot the running time. In practice, on my laptop (MacBook Pro 2017), it takes about 1'35" for a set of 15 asleep robots, and 5'30" for 16 asleep robots, i.e., for 17 robots. Note that 17! = 355,687,428,096,000 ~ 4.1 days at 1 GHz. This confirms that the algorithm is significantly faster than n! or n^n. 
    
    Here is a short description of the (recursive) algorithm OPT(), for waking up any point set from any point r, and for any distance function dist(), which is not necessarily metric. We use OPT(r,A,b), which returns the optimal wake-up time for waking up the subset A of points from b∊{1,2} robots woken up and placed at r.

    For any point r, set A not containing r, and b∊{1,2}:

        OPT(r,A,b):
        | OPT(r,A,b) = 0, if |A| = 0
        | OPT(r,A,1) = min_{u∊A} { dist(r,u) + OPT(u,A\{u},2) }, if |A| > 0
        | OPT(r,A,2) = min_{A1,A2} max{ OPT(r,A1,1), OPT(r,A2,1) }

            where {A1,A2} is a partition of A. If dist() is metric, then we
            can force |A1|=0 => |A2|=1 (or the opposite) so that the nodes
            with only one child in the wake-up tree necessarily lead to a leaf.
            A distinction must then be made between the cases |A|=0, |A|=1, |A|>=2.

    The program is based on this formulation with some optimizations. It has been programmed under the Python 3.7.3 64-bits version.

    -- Cyril Gavoille, 09/2024.
"""

from math import pi, atan2, sin, cos, floor, sqrt
from functools import lru_cache
from sys import setrecursionlimit
from time import process_time
from random import random, randrange, seed, vonmisesvariate

import os
import matplotlib as mpl
import matplotlib.pyplot as plt

# remove marks/keys in front of legends (not 100% accurate)
mpl.rcParams['legend.markerscale'] = 0
mpl.rcParams['legend.handlelength'] = 0
mpl.rcParams['legend.handleheight'] = 0
mpl.rcParams['legend.frameon'] = False

# set max recursion depth for @lru_cache()
setrecursionlimit(10*8) # don't enlarge too much, otherwise a 'segmentation fault' will occur


"""
DEFINITIONS & CONVENTIONS

What we call a ** rooted tree T ** is a simply a non empty list [ r, T_1, ..., T_k ] where r is the root of T and T_i is the i-th child of T, that is itself a rooted tree. So, for instance, [ r ] is a tree with only one node (the smallest possible tree), and [ r, [u], [v] ] is a tree with only two children (in this case, with the two leaves u and v). Note that a rooted tree has to contain at least one node, there is no special representation for an empty tree. A ** rooted tree T on a list P of points ** is simply a rooted tree whose vertices are indices of P, so integer in [0,|P|).
"""

##########################
#
# Some distance functions
#
##########################


def dist_L2(p,q):
    """
    Return the standard Euclidean distance (L2 norm) between points p and q, a point being a couple or a list with two elements. We may alternatively use math.dist(p,q) which does seem to be fater.
    """
    return sqrt( (p[0]-q[0])**2 + (p[1]-q[1])**2 )

def dist_L1(p,q):
    """
    Like dist_L2(p,q), but for the L1 norm.
    """
    return abs(p[0]-q[0]) + abs(p[1]-q[1])

def dist_Linf(p,q):
    """
    Like dist_L2(p,q), but for the L∞ norm.
    """
    return max(abs(p[0]-q[0]), abs(p[1]-q[1]))

def dist_Lp(p,q):
    """
    Like dist_L2(p,q), but for L_p norms, where p>=1 is given by the constant LP_NORM.
    """
    return ( abs(p[0]-q[0])**LP_NORM + abs(p[1]-q[1])**LP_NORM )**(1/LP_NORM)

def dist_polygon(p,q):
    """
    Distance defined by a regular k-gon oriented so that the first of the k isosceles triangles of the k-gon has a horizontal side, where k = NGON. So for k = 4, this is dist_L1(), not dist_Linf(). For k = 6, this is the hexagonal distance. In order to define a valid norm, NGON must be an even number.

    Explanation. First we calculate the norm of the difference vector v = p-q, which we express in polar coordinates rho(v) and arg(v). We reduce to the case where v falls within the first isosceles triangle by calculating arg(v)%(2𝜋/k). This vector is then projected onto the perpendicular bisector of the triangle using cos(t-arg(v)). Normalise the distance with the height of the triangle, which is cos(𝜋/k).
    """
    dx, dy = p[0]-q[0], p[1]-q[1] # coordinates of vector p-q
    t = pi/NGON # half-angle of a triangle
    return sqrt(dx*dx+dy*dy) * cos(t - (atan2(dy,dx)%(2*t))) / cos(t)

# @cache() is actually faster than @lru_cache() from Python 3.9
@lru_cache(maxsize=None)
def dist(i,j):
    """
    Compute the distance between POINTS[i] and POINTS[j] according to the constant DIST, the distance function. We cache this function, because its intensively use.
    """
    return DIST( POINTS[i], POINTS[j] )


#############################################################
#
# Global variables to avoid parameter passing whenever they
# are constant thru recursive calls. This becomes necessary
# with the use of @lru_cache().
#
#############################################################


DIST = dist_L2   # default distance function, used by dist()
POINTS = []      # set of points, used by dist()
ROOT = 0         # index in POINTS of the source, i.e., the awake robot
UPPER_BOUND = -1 # upper bound used in branching algorithms
SEED = None      # used to retrieve the seed for random generator
PROCESS_TIME = 0 # elapsed time in seconds for display in draw_all()
NGON = 6         # number of sides for dist_polygon() distance function
OPTSYM = False   # variant for opt() allowing symmetry optimization
LP_NORM = 2      # constant p for dist_Lp() distance function


######################################
#
# eccentricity(), diameter(), depth()
# 
######################################


def eccentricity(r=None, P=None, d=None):
    """
    Compute the excentricity of point P[r] in P according to the distance function d(). The algorithm has complexity O(|P|).

    Global variables modified: POINTS and DIST for the later use of dist()
    """
    global POINTS, DIST # important for dist()
    if r == None: r = ROOT # default value, does not work otherwise
    if P == None: P = POINTS # default value, does not work otherwise
    if d == None: d = DIST # default value, does not work otherwise
    POINTS, DIST = P, d # set values for dist()

    x = 0 # value to return
    for v in range(len(P)): # for all points of P
        x = max(x,dist(r,v))
    return x

def diameter(P=None, d=None):
    """
    Compute the diameter of P according to the distance function d(). We do not suppose that d is symmetric. The algorithm has complexity O(|P|^2).

    Global variables modified: POINTS and DIST because eccentricity()
    """
    if P == None: P = POINTS # default value, does not work otherwise
    if d == None: d = DIST # default value, does not work otherwise

    x = 0 # value to return
    for r in range(len(P)): # for all points of P
        x = max(x,eccentricity(r,P,d))
    return x

def depth(T, P=None, d=None):
    """
    Compute the depth of a rooted tree T on point set P according to the distance function d(). The algorithm has complexity O(|P|).

    Global variables modified: POINTS and DIST for the use of dist()
    """
    global POINTS, DIST # important for dist()
    if P == None: P = POINTS # default value, does not work otherwise
    if d == None: d = DIST # default value, does not work otherwise
    POINTS, DIST = P, d # set values for dist()

    if len(P) == 0 or len(T) <= 1: return 0
    r = T[0] # root of T
    if r < 0 or r >= len(P): return 0 # point outside of P

    x = 0 # value to return
    for i in range(1,len(T)): # for all subtree of T
        f = T[i] # i-th child of T
        x = max(x, dist(r,f[0]) + depth(f,P,d))
    return x


###########################
#
# Some point distributions
#
###########################


def init_seed(n=None):
    """
    Set to n (or not, if n=None) for random generating functions. The reason is to allow us to run again random based point sets, to replay a very bad case. Warning! Don't use this function too many times, in a for loop for instance, since a short seed (4 digits) is used. This is more practical to notice it by hand.
    """
    global SEED
    SEED = randrange(10000) if n == None else n
    seed(SEED)

def generate_regular_polygon(n, r=1, c=(0,0)):
    """
    Generates a list of n points placed regularly on a circle of radius r and centred at c, i.e. a regular n-gone. The points are ordered anti-clockwise from (1,0), if n>0.
    """
    P = [] # list to return
    if n > 0: theta = 2*pi/n # used only if n>0
    for i in range(n): # repeat n times
        P += [ (c[0]+r*cos(i*theta), c[1]+r*sin(i*theta)) ]
    return P

def generate_von_mises(n, p=0.5, k=1, f=0):
    """
    Generate a list of n random points in the disc centred on the origin where the distance to the centre is x^p with x is uniform over [0,1] and where the angle follows a circular normal distribution (more precisely a von Mises distribution) of mean angle a_i = i*2𝜋/k where i is a uniform integer over [0,k[. We need k >= 1. The parameter f is related to the standard deviation of the angle from the direction with angle a_i. More precisely, the standard deviation is defined as sigma = (a_{i+1} - a_i)/(2f) = 𝜋/(k*f). So the larger f is, the more angles are concentrated in a_i. We then arrange the standard deviation of the von Mises law so that it corresponds to sigma. To do this, we set the law parameter kappa = 1/sigma^2, which is a very good approximation (see https://dlwhittenbury.github.io/ds-1-sampling-and-visualising-the-von-mises-distribution.html).
    
    If p=0, then the distance is unitary and the points are all on the unit circle. If f=0, then the distribution of the angle will be uniform (infinite variance), whatever the value of k. If f=+∞ (or any value of f<0), then the variance will be zero (kappa=+∞) and the possible angles will all be a_i. Set k=0 to indicate k=+∞, which will have a similar effect to setting f=0.
    
    Note that p=0.5 does not necessarily give a uniform distribution over the 'flower', i.e. the disc deformed by the k directions, especially when k is small and f is large enough.

    Examples:
    - generate_von_mises(n,0)        -> random pts uniformaly on the circle
    - generate_von_mises(n,0.5)      -> random pts uniformaly on the disk
    - generate_von_mises(n,1,k,-1)   -> random pts uniformaly on a star with k branchs
    - generate_von_mises(n,0.1,k,3)  -> random pts concentrate towards the k directions
    - generate_von_mises(n,-0.1,k,3) -> idem, but points are outside the circle

    Note that points are not necessarily ordered according to the angle at the origin. If you need this ordering, you need to proceed as follows:

        P = generate_von_mises(...)
        P.sort(key=lambda A:atan2(*A))

    If you need points that looks uniform on unit disc (or circle) w.r.t. some distance function d() that is not dist_L2(), you can proceed as follows:
    
        P = generate_von_mises(n,0.5) # replace 0.5 by 0 for the unit circle
        distortion(P,dist_L2,d)

    """
    P = [] # list to return
    if f < 0: kappa = float("inf") # a uniform distribution (deviation ∞)
    else: kappa = (k*f/pi)**2 # standard deviation in 𝜋/(kf)
    if k == 0: k,kappa = 1,float("inf")

    for _ in range(n): # repeat n times
        i = randrange(k) # one of the k random directions
        a = vonmisesvariate(i*2*pi/k,kappa) # random angle around the chosen direction
        x = random()**p # random distance uniform in [0,1]
        P += [ (x*cos(a), x*sin(a)) ] # add the new point at the end of P

    return P

def generate_convex(n):
    """
    Generates a list of n random points in convex position. The algorithm is in O(n*log(n)) due to sorting.

    PRINCIPLE.
    
        Start with a list of n random points in [0,1[², then calculate the difference (i.e. the vector) between two consecutive points. The sum of the n vectors is zero. The vectors are sorted by angle, then the points of the convex envelope are drawn in ascending order.
    """
    P = []
    for _ in range(n): # n random points in [0,1[²
        P += [(random(), random())]
    
    for i in range(n-1): # compute vector differences
        P[i] = (P[i][0] - P[(i+1)%n][0], P[i][1] - P[(i+1)%n][1])

    P.sort(key=lambda A: atan2(*A)) # sort according angles

    for i in range(1,n): # final positions
        P[i] = (P[i][0] + P[i-1][0], P[i][1] + P[i-1][1])

    return P

def normalize_bc(P=None, d=None):
    """
    Translate and scale a given list of points P by adding their barycentre at the end of P and such that the eccentricity from this barycentre is 1. ROOT is then modified so that it corresponds to this barycentre. The points of P are first recentred so that the barycentre is placed at the origin (0,0), then the coordinates are scaled so that the eccentricity of this point is 1. Relative distances are preserved if the distance function d() is linear, i.e., if d(λA+C, λB+C) = d(A,B) for all points A,B,C and real λ>0. This is the case for all norms, in particular L_p for all p>=1.

    Global variables modified: POINTS and DIST because of eccentricity(), ROOT
    """
    if P == None: P = POINTS # default value, does not work otherwise
    if d == None: d = DIST # default value, does not work otherwise
    dist.cache_clear() # important for eccentricity()

    # calculates the barycentre (sx,sy)
    n = len(P)
    if n == 0: return # nothing to do in this case
    x = y = 0 # barycentre
    for i in range(n): # for all points in P
        x,y = x+P[i][0], y+P[i][1]
    x,y = x/n, y/n
    global ROOT
    ROOT = n # modifies the root
    P += [(x,y)] # add the barycentre at the end of the list
    ex = eccentricity(ROOT,P,d) # calculates the eccentricity of the barycentre
    if ex == 0: ex = 1 # in this case, we put all the points at (0,0)
    for i in range(n+1): # also transforms the barycentre
        P[i] = ((P[i][0]-x)/ex, (P[i][1]-y)/ex)
    return

def distortion(P, d1, d2):
        """
        Transform each point A of P into a new point A' so that d2(O,A') = d1(O,A) and O,A,A' are aligned, where O = (0,0) is the origin. In other words, the points of P are scaled with respect to the origin so that the distances are the same according to d1() before the transformation and d2() after the transformation. For example, if P is a set of points on the L2 circle, then after distortion(P,dist_L2,dist_L1), all the points of P will be on L1 circle (i.e. on a square). The point A' is set to O, if distances d1(O,A) or d2(O,A) are not positive.
        """
        O = (0,0)
        for i in range(len(P)):
            A = P[i]
            t1,t2 = d1(O,A), d2(O,A)
            r = 0 if t1 <= 0 or t2 <= 0 else t1/t2
            P[i] = (r*A[0], r*A[1])

def sample_circle(d=DIST, n=256, radius=1):
    """
    Return a sample of 4n points on a circle w.r.t. the distance function d(). It is centered at the origin and has radius r. This is a 4n-gon approximation.
    """
    P = [] 
    for i in range(1,n+1): # generate P with 4n points on the L1-circle
        x = radius*i/n
        P += [(x,radius-x),(radius-x,-x),(-x,x-radius),(x-radius,x)]

    P.sort(key=lambda A:atan2(*A)) # sort points according to angle
    distortion(P,dist_L1,d) # points are now on the d() circle of radius r
    return P

def half_perimeter(d=DIST, n=256):
    """
    Return an approximation (actually a minorant) of half of the perimeter of the unit disk w.r.t. distance function d(). The number n defines the accurency (the number of sides used un the inscribed polygon, namely 4n). The accurency is about 10^(-2*log_10(n)). In other words, the number of correct digits in this approximation is about 2*log_10(n).

    Ex:
    
        print(half_perimeter(dist_L2),pi)
        > 3.1415871144373253 3.141592653589793

        print(half_perimeter(dist_L1),pi)
        > 4.0

        NGON = 6
        print(half_perimeter(ft.dist_polygon))
        > 3.0        

    """
    P = sample_circle(d,n,n) # enlarge the radius of the disk to avoid computation errors
    s = 0 
    for i in range(0,len(P)): # compute the sum of polygon side length 
        s += d(P[i-1],P[i]) # starts with d(P[-1],P[0]) = d(P[0],P[n-1])
    return s/(2*n) # divide by 2 because it's half the perimeter


############################################
#
# Drawing functions: disc, points, trees, ...
#
############################################


def draw_disc(d=DIST, n=256, radius=1, circle='black', fill='powderblue', lw = 1, alpha=0.3):
    """
    Draw a circle (or a disc) in the current plt object w.r.t. the distance function d(). It is centered at the origin and has radius r. This is a 4n-gon approximation. We actually sample 4n points regularily spaced in the L1 circle, so at uniform L1 distance. These 4n points are linked by segments with color *circle* and linewidth lw. If *circle* = None, then the color of the circle, i.e., the border of the disc, is set of color given by *fill*. If *fill*='None' (with the quotes), the disc is not fill, otherwise it is filled with color *fill* and with transparency given by *alpha* in [0,1].

    Ex:

    def f(p,q):
        dx, dy = abs(p[0]-q[0]), abs(p[1]-q[1])
        cx, cy = 2.3, 1.2
        return ( dx**cx + dy**cy )**(1/(cx*dx+cy*dy))
    
    plt.clf()
    plt.axis('equal')
    draw_disc(f)
    plt.show()

    """
    P = sample_circle(d,n,radius) # sample 4n points of the circle of radius r
    if circle == None: circle = fill
    x = [ P[i][0] for i in range(len(P)) ] # list of x-coordinates
    y = [ P[i][1] for i in range(len(P)) ] # list of y-coordinates
    plt.fill(x, y, facecolor = fill, edgecolor=circle, linewidth=lw, alpha = alpha)
    
def draw_points(P, c='gray', s=5, m='o'):
    """
    Draw a set of points P, each with colour c, size s, and shape m, in the current plt drawing. For a simple display of the points, a following extra plt.show() or draw_all() is required.
    """
    for (x,y) in P:
        plt.plot([x], [y], color = c, marker = m, ms = s, ls = 'none')

def draw_tree(T, P=None, c_edge='blue', w=0.005, arc=False, pts=True, c_pts='red', s=5, long=False, c_long='green'):
    """
    Recursively draw a rooted tree T on a list of points P.

     c_edge = color of the edges
          w = width of the edges
        arc = arc from the parent to the child (True) or unoriented edge (False)
        pts = draw the points of the tree (True) or not (False)
      c_pts = color of the points (if pts=True)
          s = size of the points (if pts=True)
       long = draw edges/arcs of the longest branch (according to DIST)
     c_long = color of the longest branches (if long=True)

    Global variable modified (if long=True): POINTS (for dist())
    """
    if long:
        global POINTS
        if P == None: P = POINTS
        POINTS = P
    if len(P) == 0 or len(T) == 0: return # nothing to do in this case
    r = T[0] # the root of T
    if r<0 or r>len(P): return # index out of range of P
    r = P[r] # r = coordinate of the root
    if pts: # draw the root r
        plt.plot([r[0]], [r[1]], color = c_pts, marker = 'o', ms = s, ls = 'none')
    hl, hw = (0.1, 0.05) if arc else (0,0) # length and width of the arcs
    if long: ex = depth(T) # eccentricity of r

    for i in range(1,len(T)): # for each subtree

        # recursively draw subtree T[i]
        draw_tree(T[i],P,c_edge,w,arc,pts,c_pts,s,long,c_long)

        # draw the arc between root r and its subtree T[i]
        f = T[i][0] # root of the child
        if f<0 or f>len(P): return # index outside P
        c = c_long if long and ex == dist(T[0],f) + depth(T[i]) else c_edge # if long branch
        f = P[f] # coordinates of f
        lss, col = (':','orange') if T[0] == ROOT else ('-',c) # dash line for the first edge
        plt.arrow(r[0], r[1], f[0]-r[0], f[1]-r[1], color = col, width = w, head_width = hw, head_length = hl, length_includes_head = True, ls = lss, overhang = .2)

def draw_all(title=None, x=None, T=None, P=None, d=None, save=None, xaxis=[], yaxis=[], disc=True):
    """
    Draw and display a time solution x with its rooted tree T on a set of points P, title being just a comment for the legend. It also calculates and displays the depth of T according to the distance function d(). A simple call to draw_all(), with no arguments, displays the points of POINTS, its eccentricity and its diameter. Setting ROOT=None prevents the root from being distinguished from the other points. To avoid displaying the 'time' line, set PROCESS_TIME = 0. To avoid displaying the 'seed' line, set SEED = -1. Warning! PROCESS_TIME is reset to zero by this function. A unit disc (w.r.t. d() function) is displayed iff *disc* is True.

    Global variables modified: POINTS, DIST (via eccentricity()), PROCESS_TIME
    """
    global PROCESS_TIME
    if P == None: P = POINTS # default value, does not work otherwise
    if d == None: d = DIST # default value, does not work otherwise
    dist.cache_clear() # needed for the correct use of dist(), eccentricity(), and diameter()

    n = len(P) # number of points
    s = max(1, 6 - floor(n/50)) # variable point size
    w = s/1000 # variable edge width
    arc = True if n < 60 else False # orientation or not of the edges
    r = ROOT if T == None else T[0] # root, possibly of the tree if it exists

    plt.clf() # clear the current plt drawing
    plt.axis('equal') # aspect ratio: to have circle with real circle shape
    if len(xaxis) == 0 and len(yaxis) == 0: plt.axis('off') # remove axis
    for z in xaxis: plt.axvline(z,zorder=0)
    for z in yaxis: plt.axhline(z,zorder=0)
    
    # NB. Objects are displayed one on top of the other, so one draws the root of T last to see it, even if there are a lot of points. Having said that, when the grid axis are drawn first, they are displayed on top of the points, and I don't know why. Perhaps, one should set zorder in a more clever way?

    if disc: draw_disc(d) # draw the unit disc
    draw_points(P,'gray',s) # draw the points in gray and with size s
    if T != None: # if the tree does exist
        draw_tree(T,P,'blue',w,arc,True,'red',s,True) # draw the tree with arcs to children
    if ROOT != None: # if the root does exist
        draw_points([P[r]],'orange',1.1*s,'o') # draw the root

    # compute the legend L to display
    L = []; t = ""
    if title != None: t += title # add title
    if x != None: t += f" = {x:.3f}" # add the value x
    if t != "": L += [t]
    if T != None: L += [f"depth = {depth(T,P,d):.3f}"] # add the depth
    L += [f"#pts = {n-1}+1"] # add the number of points
    L += [f"diam = {diameter(P,d):.3f}"] # add the diameter information
    if ROOT != None: L += [f"ecc = {eccentricity(r,P,d):.3f}"] # add the eccentricity information
    if PROCESS_TIME > 0: # add the current process time information
        if 0 <= s < 10:     L += [f"time = {PROCESS_TIME:.3f}s"]
        if 10 <= s < 100:   L += [f"time = {PROCESS_TIME:.2f}s"]
        if 100 <= s < 1000: L += [f"time = {PROCESS_TIME:.1f}s"]
        if 1000 <= s:       L += [f"time = {PROCESS_TIME:.0f}s"]
        PROCESS_TIME = 0 # reset to avoid to display the same time again and again
    if SEED >= 0: L += [f"seed = {SEED}"] # add the SEED information
    plt.figlegend(L) # display the legend L

    # save the figure in current directory
    if save != None:
        file = os.getcwd() + '/' + save
        plt.savefig(file, format='svg', dpi=1200) # to put before plt.show()
        print(f"File {file} saved.")

    plt.show() # show drawing


#####################################################
#
#  Algorithm optimal_tree(r,P,d,s):
#
#    r = root (index in P of the awake robot)
#    P = point set
#    d = distance function
#    s = symmetry optimization if True
#
#  It uses symmetry(r,A) and opt(r,A,b) sub-routines.
#####################################################


@lru_cache(maxsize=None)
def symmetry(r,A):
    """
    Construct a subset V of A without certain symmetries, (r,A) being an intance of the function opt(r,A,b). More precisely, let us denote D(i,j) the distance vector between i and all the points of A\{i,j}, ordered according to the increasing indices of A. We can calculate D(i,j) = [ dist(i,u) for u in A if u!=i and u!=j ].
    
    The subset V returned (actually a tuple) has the property that there are no i,j in V such that [ dist(r,i), D(i,j) ] = [ dist(r,j), D(j,i) ]. This is because the optimal trees for r -> i -> A\{i,j} and r -> j -> A\{i,j} will be identical. So either i or j can be deleted. This is linked to the fact that the distances between the points of A\{i,j} are the same for D(i,j) and D(j,i). Unfortunately, very few (r,A) instances have such symmetries. Even uniform points on a circle do not have these symmetries because the distances between points of index [1,2,3,4] for example, will not be the same from i or j, even if these 4 points are symmetrical around the origin. More symmetries internal to A\{i,j} should be taken into account.
    """
    
    n = len(A)
    V = list(A)[:] # list to construct, copy of A
    D = [(dist(r,u),u) for u in A] # distance vector between r and A, according to dist() function
    D.sort(key=lambda d:d[0]) # decreasing distance sorting
    B = [] # list of blocks
    d0 = D[0][0]+1

    for i in range(n): # determine blocks of D where symmetries has to be seek
        if D[i][0] < d0:
            B += [i] # start/end of a block
            d0 = D[i][0] # new distance
    B += [n] # last block

    for i in range(1,len(B)): # each bloc of D
        L = [d[1] for d in D[B[i-1]:B[i]]] # L = list of vertices equidistant from r
        m = len(L)
        # for each pair {u,v} of L, optionally removes v from L
        for i in range(m):
            u = L[i]
            for j in range(i+1,m):
                v = L[j]
                Du = [ dist(u,w) for w in A if w!=u and w!=v ]
                Dv = [ dist(v,w) for w in A if w!=u and w!=v ]
                if Du == Dv: V.remove(v) # v is useless in V

    return tuple(V)

@lru_cache(maxsize=None)
def opt(r,A,b):
    """
    Returns the optimal wake-up time and tree for waking up the subset A (subset of POINTS indices) from 1 or 2 robots woken up and positioned at r (an index of POINTS). 
    If b=True, only 1 robot is awakened (and r has a single son in the tree), otherwise 2 are awakened (and r has two children). 
    The subset A must neither contain r nor be empty. We assume the metric case, i.e. DIST satisfies the triangle inequality.
    
    Be careful! A is here a tuple, not a list, so that it is 'hashable', which is important for @lru_cache(). 
    One property (which is not used) is that all calls to opt(r,A,b) are such that the indices of A are arranged in ascending order. 
    It can be shown that the complexity of opt(r,A,b), for a subset A of POINTS, is about O( |A|*|POINTS| * 2^|A| * binom(|POINTS|,|A|), up to some polynomials in |POINTS|.

    STRATEGY

        CASE 1 (b=True): choose the best target v in A, then wake up all the points of A\{v} with two robots placed at v by applying opt(v,A\{v},False and thus create a root tree r with a son, v, and its subtree.

        CASE 2 (b=False): cut A into two subsets, A1 and A2, each with at least one element, and apply the strategy opt(r,A1,True) and opt(r,A2,True) to create a tree rooted at r with two children. 
        Of course, all the partitions {A1,A2} must be tested and the best tree selected. 
        In the non-metric case, we would also have to test the case where a single robot moves to wake up the whole of A (A1=A and A2={}). 
        This is unnecessary in the metric case.
    
    For the complexity, that is exponetial in n, we will concentrate on the most important term: the exponent base of n. Let n=|POINTS|, a=|A|, and x = a/n. 
    We want to show that the time complexity C(a,n) = 2^a * binom(n,a) satisfies C(a,n) = O(3^n). If a = o(n), then clearly C(a,n) < 2^n. 
    It is well known that, for any constant ratio x = 0..1, binom(n,xn) ~ h(x)^n / √(2𝜋x(1-x)n) where h(x) = (1/x)^x * (1/(1-x))^(1-x). 
    As a result, C(xn,n) = 2^(xn) * binom(n,xn) ≃ (2^x * h(x))^n. This is maximum whenever x=2/3. 
    In short, we can bound the complexity O( a*n * 2^n * binom(n,a) ) by O( n^{3/2} * 3^n ) (don't forget the 1/√n in binom(n,xn)). 
    Note that a log(n) term is probably missing accessing data in the cache (memorisation), although it is probably possible to implement it in O(1) time, at least in a damped way.

    NB. Care must be taken when managing lists and not to generate too many of them. Subsets could be coded as integers, perhaps to increase speed and make more efficient use of the cache. In any case, the presence of tuple(A)/list(A) in recursive calls is necessary to avoid the 'TypeError: unhashable type: “list”' error (you can't cache function calls that contain lists, but it's OK for tuples).
    """

    n = len(A)
    if n == 1: return dist(r,A[0]), [r, list(A)] # it's over: r has a child that is a leaf

    # optimization assuming symmetric dist(), but there is no significant gain (often, this is worst ...)
    # if n == 2: 
    #    if dist(r,A[0]) < dist(r,A[1]): return dist(r,A[0])+dist(A[0],A[1]), [r, [ A[0], [A[1]] ]]
    #    return dist(r,A[1])+dist(A[1],A[0]), [r, [ A[1], [A[0]] ]]
    
    xmin = UPPER_BOUND # a (strict) majorant over the result

    if b:
    
    ###########################
    # CASE 1: one awake robot #
    ###########################

        if OPTSYM and n>3: # search for symmetry, if |A| is large enough
            V = symmetry(r,A)
        else: # no search
            V = A

        # find the best point v in V to wake up, NB. Here |V|>=2
        for v in V: # for each possible v in V
            B = tuple(u for u in V if u != v) # V = V\{v}
            x,T = opt(v,B,False) # computation of the solution for (v,V), NB. Here |V|>=1
            x += dist(r,v) # add time r->v
            if x < xmin: # we have found a strictly better solution starting with v
                xmin,Tmin = x,T # keep the best solution

        return xmin, [r, Tmin] # tree with one child

    ############################
    # CASE 2: two awake robots #
    ############################

    # In this case, we need to divide A into two subsets, A1 and A2, each with at least 1 point, which is possible because here |A|>=2. To do this, we represent a subset of A by a binary word w of |A| bits. The 1 bits of w are the elements of A1, the 0 bits those not of A1 (thus of A2). We avoid w=000...00 and w=111...11 to guarantee |A1|,|A2|>=1. Since the {A1,A2} and {A2,A1} partitions produce an equivalent tree, we avoid the possibility of having A1 and then bar(A1), its complementary binary word. We can therefore assume that for A1, we go from w=000...01 to 011...11. Finally, we will guarantee that |A1|>=|A2| in order to cut faster (A1 being recursively launched before A2) because a priori it is for large sets that the tree has the best depth.

    for w in range(1,2**(n-1)): # w = 000..01 to 011..11, |w|=n

        # construct A1 and A2 according to w
        A1 = []; A2 = [] # do not put A1=A2=[] ...
        for i in range(n):
            if w & 1: # test the last bit of w
                A1 += [A[i]]
            else: A2 += [A[i]]
            w >>= 1 # remove the last bit from w
            # NB. This does not modify the w in the loop 'for w in range(...)'.

        if len(A1)<len(A2): A1,A2 = A2,A1 # NB. Here |A1|>=|A2|>=1

        x1,T1 = opt(r,tuple(A1),True) # 1 robot wakes up A1
        if x1 >= xmin: continue # no need to continue, xmin is the best tree
        x2,T2 = opt(r,tuple(A2),True) # 1 robot wakes up A2
        if x2 >= xmin: continue # no need to continue, xmin is the best tree
        # we've found a better tree
        xmin,T1min,T2min = max(x1,x2),T1,T2 # the larger of the two branches

    # Here T1min and T2min have the same root r. They must be merged.
    # Example: if T1min=[r,T1'] and T2min=[r,T2'] (one child each)
    # then we should return [r, T1', T2'] (two children for r)

    return xmin, [r, T1min[1], T2min[1]]

def optimal_tree(r=None, P=None, d=None, s=False):
    """
    Returns the optimal wake-up time and tree for a set of points P from an awake robot r, which must be an index of P. Therefore |P| >= 1 and r in [0,|P|[. The list P is not modified. We assume the metric case, i.e., d() (or DIST) satisfies the triangle inequality, but this is easy to modify to work without this asymption. If s=True, then a symmetry search is performed in CASE 1 of opt(r,P,d).
    
    The complexity of the algorithm is majored by the call to opt(r,A,b), i.e. by O( n^{3/2} * 3^n ) where n = |P|. In practice, for n = 15+1 points, this takes 1'35" and for n = 16+1 points, this takes 5'30". Note that 17! = 355,687,428,096,000 ~ 4.1 days at 1 GHz). This confirms that the algorithm is significantly faster than n! or n^n. 

    Global variables modified: ROOT, POINTS, DIST, OPTSYM, UPPER_BOUND, PROCESS_TIME.
    """
    global ROOT, POINTS, DIST, OPTSYM, UPPER_BOUND, PROCESS_TIME
    if r == None: r = ROOT    # default value, does not work otherwise
    if P == None: P = POINTS  # default value, does not work otherwise
    if d == None: d = DIST    # default value, does not work otherwise
    ROOT = r    # set global variable for opt()
    POINTS = P  # set global variable for opt()
    DIST = d    # set global variable for opt()
    OPTSYM = s  # set global variable for opt()
    UPPER_BOUND = 100*eccentricity(r,P,d) # majorant (strict) on the best solution
    opt.cache_clear()   # important! otherwise results might be incorrect
    dist.cache_clear()      # important! otherwise results might be incorrect
    symmetry.cache_clear()  # important! otherwise results might be incorrect
    
    PROCESS_TIME = process_time() # for the elapse time

    if len(P) == 1: return 0, [r] # no robot to wake up, point r is already wake up

    A = list(range(len(P))) # A = [0,|P|[ = list of indices of P
    del A[r] # A = [0,r[ u ]r,|P|[ = list of indices of P except r
    x,T = opt(r,tuple(A),True) # here r must wake up A, NB. |A|>=1

    PROCESS_TIME = process_time() - PROCESS_TIME # duration

    return x,T
