import numpy as np
import random as rd
import matplotlib.pyplot as plt
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def diamant_carre(n):
    m = 2 ** n + 1
    matrice = np.zeros((m,m))
    # valeur aléaoires pour les 4 coins 
    matrice[0,0] = rd.uniform(-m,m)
    matrice[0,m-1] = rd.uniform(-m,m)
    matrice[m-1,m-1] = rd.uniform(-m,m)
    matrice[m-1,0] = rd.uniform(-m,m)
    
    pas = m - 1
    while pas > 1:
        pas2 = pas // 2
        # phase du diamant
        for x in range(pas2 , m, pas ):
            for y in range(pas2 , m, pas ):
                moy = (matrice[x - pas2 , y - pas2] + matrice[x - pas2 , y + pas2] + matrice[x + pas2 , y + pas2] + matrice[x + pas2 , y - pas2]) / 4
                matrice[x,y] = moy + rd.uniform(-pas2,pas2)
        #phase du carre
        for x in range(0,m,pas2):
            if x % pas == 0:
                decalage = pas2
            else:
                decalage = 0
            
            for y in range(decalage , m,pas):
                somme = 0
                i = 0
                if x >= pas2:
                    somme += matrice[x - pas2, y]
                    i += 1
                    
                if x + pas2 < m:
                    somme += matrice[x + pas2, y]
                    i +=1
                
                if y >= pas2 :
                    somme += matrice[x , y - pas2]
                    i +=1
                
                if y + pas2 < m:
                    somme += matrice[x , y + pas2]
                    i += 1
                
                matrice[x , y] = somme / i + rd.uniform(-pas2, pas2)
        
        pas = pas2
    return matrice
    
    
def map2D(n):
    matrice = diamant_carre(n)
    m = 2 ** n + 1
    img = Image.new('RGB', (m , m),color = (0,0,0) )
    pixel = img.load()
    
    for x in range(m):
        for y in range(m) : 
            if matrice[x,y] < 0 :
                k = - matrice[x,y]
                pixel[x , y] = ( 0, 0 , int(255 - (255 / m ) * k ))
            else : 
                k = matrice[x,y]
                pixel[x,y] = (0 , int((255 / m) * k),0)
                
    img.show()
    
def preMap3D(n):
    matrice = diamant_carre(n)
    m = 2 ** n + 1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X = np.arange(1, m + 1, 1)
    Y = np.arange(1, m + 1, 1)
    X , Y = np.meshgrid( X , Y)
    Z = matrice.reshape(X.shape)

    ax.plot_surface(X, Y, Z)

    plt.show()
     
def alea3D(n):
    m = 2 ** n + 1
    matrice = m * np.random.rand(m,m)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X = np.arange(1, m + 1, 1)
    Y = np.arange(1, m + 1, 1)
    X , Y = np.meshgrid( X , Y)
    Z = matrice.reshape(X.shape)

    ax.plot_surface(X, Y, Z)

    plt.show()
    
def alea3Dcouleur(n):
    m = 2 ** n + 1
    matrice = m * np.random.rand(m,m)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X = np.arange(1, m + 1, 1)
    Y = np.arange(1, m + 1, 1)
    X , Y = np.meshgrid( X , Y)
    Z = matrice.reshape(X.shape)

    ax.plot_surface(X, Y, Z, cmap=cm.gist_earth)

    plt.show()

def alea2D(n):
    m = 2 ** n + 1
    matrice = m * np.random.rand(m,m)
    img = Image.new('RGB', (m , m),color = (0,0,0) )
    pixel = img.load()
    
    for x in range(m):
        for y in range(m) : 
            if matrice[x,y] < 0 :
                k = - matrice[x,y]
                pixel[x , y] = ( 0, 0 , int(255 - (255 / m ) * k ))
            else : 
                k = matrice[x,y]
                pixel[x,y] = (0 , int((255 / m) * k),0)
                
    img.show()
    
def Map3D(n):
    matrice = diamant_carre(n)
    m = 2 ** n + 1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X = np.arange(1, m + 1, 1)
    Y = np.arange(1, m + 1, 1)
    X , Y = np.meshgrid( X , Y)
    Z = matrice.reshape(X.shape)

    ax.plot_surface(X, Y, Z, cmap=cm.hot , linewidth = 0)

    plt.show()
    
def Map2D3D(n):
    Map3D(n)
    map2D(n)
