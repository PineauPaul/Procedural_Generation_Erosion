import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import random as rd
from scipy import interpolate

    
def Height_Map_FFT(f_hurst,n):
    def regulateur(x, y):
        if x == 0 and y == 0:
            return 0.0
        return f_hurst(np.sqrt(x**2 + y**2))
            
    noise = np.fft.fft2(np.random.rand(n,n)) #O(nlog(n))    
    amplitude = np.zeros((n,n))
        
    for i in range(n):
        for j in range(n):            
            amplitude[i, j] = regulateur(i, j)
                
    return np.fft.ifft2(noise * amplitude)

def HMFFT3(n,b):
        
    noise = np.fft.fft2(np.random.normal(size = (n,n))) #O(nlog(n))

    for i in range(n):
        for j in range(n):            
            f = np.sqrt(((i - n/2)/n) ** 2 + ((j-n/2)/n)**2)
            if f < 1/n :
                f= 1/n
            noise[i,j] *= 1/(pow(f,b))
            
    return np.absolute(np.fft.ifft2(noise))

def Height_Map_FFT_gauss(n,sigma):
    def regulateur(x, y):
        return (1/(2 * np.pi * sigma **2)) * np.exp(-(x**2 + y**2)/(2*sigma ** 2))
        
    noise = np.fft.fft2(np.random.normal(size = (n,n))) #O(nlog(n))
    amplitude = np.zeros((n,n))
    for i in range(n):
        for j in range(n):            
            amplitude[i, j] = regulateur(i, j)
            
    return np.fft.ifft2(amplitude * noise)

def FFT(n,hurst):
    plt.close()
    out = Height_Map_FFT(lambda k: k**(-hurst),n)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    X = np.arange(1, n+1, 1)
    Y = np.arange(1, n+1, 1)
    X , Y = np.meshgrid( X , Y)
    Z = (out.real.reshape(X.shape))

    ax.plot_surface(X, Y, Z ,cmap = cm.gist_earth, linewidth = 0)

    plt.show()
    
def Water_Map_func2(n,carte):
    ''' n : taille de la carte '''
    W = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if carte[i,j] < 0 :
                W[i,j] = -carte[i,j]
    return W
                
def Water_Map_func(n):
    return np.zeros((n,n))

def Sediment_Map_func(n):
    return np.zeros((n,n))

def Flux_Map_func(n):
    return np.array([[[0,0,0,0]]*n]*n)

def pluie(Water_map,eau_max):
    carte_pluie = Water_map
    n = len(Water_map)
    for k in range(n):
        for j in range(n):
            carte_pluie[k,j] += (eau_max * np.random.random())
    return carte_pluie


def norme(vect):
    '''norme d'un vecteur x=(a,b)'''
    a,b = vect
    return np.sqrt(a ** 2 + b ** 2)

def carte(taille,hurst):
    '''renvoie une carte de hauteur générée par transformée de Fourier'''
    Z = Height_Map_FFT(lambda k: k**(-hurst),taille)
    Height_map = Z.real
    return Height_map

def plot3d(carte_avant, carte_apres,taille):
    '''affiche la carte avant et apres modification par algorithme d'érosion'''
    X = np.arange(1, taille+1, 1)
    Y = np.arange(1, taille+1, 1)
    X , Y = np.meshgrid( X , Y)
    
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1, projection='3d')
    ax.plot_surface(X, Y, carte_avant ,cmap = cm.gist_earth, linewidth = 0)

    ax = fig.add_subplot(1,2,2, projection='3d')
    ax.plot_surface(X, Y, carte_apres ,cmap = cm.gist_earth, linewidth = 0)
    plt.show()

def erosion(taille,hurst,A,l,eau,dt,Kc,Ks,Kd,Ke,n):
    '''algorithme d'érosion'''
    carte_a = carte(taille,hurst)
    carte_apres = Erosion_t_dt(taille,carte_a,A,l,eau,dt,Kc,Ks,Kd,Ke,n)
    plot3d(carte_a,carte_apres,taille)
    

def Erosion_t_dt(taille,carte,A,l,lxy,eau,dt,Kc,Ks,Kd,Ke,n):#complexité O(n*taille*taille)
    '''algorithme d'érosion calculant la carte de hauteur apres érosion '''
    X = np.arange(1, taille+1, 1)
    Y = np.arange(1, taille+1, 1)
    #l,r,t,b
    Height_map = np.copy(carte)
    Water_map = Water_Map_func2(taille,carte)
    Water_map2 = Water_Map_func(taille)
    Water_map3 = Water_Map_func(taille)
    Sediment_map = Sediment_Map_func(taille) 
    Flux_map = Flux_Map_func(taille)
    Velocity_field = np.array([[[1.0,1.0]]*taille]*taille)
    angle = np.zeros((taille,taille))
    Capacite = np.zeros((taille,taille))
    s1 = Sediment_Map_func(taille)


    debut_barre()
    prog = 0 # pour la barre de progression
    
    for k in range(n):

        
        Water_map2 = pluie(Water_map,eau)
        for x in range(taille):
                for y in range(taille):
                    Flux = np.array([1.0,1.0,1.0,1.0])
                    HW = Height_map[x,y] + Water_map2[x,y]

                    #On calcule les vecteurs flux a partir des cases du haut, du bas , des cotes donc on differencie juste les cas ou une des ces 4 cases n'existe pas
                    
                    if x == 0 :
                        if y == 0 :
                            Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW - Height_map[0,1] - Water_map2[0,1]) / l))
                            Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW - Height_map[1,0] - Water_map2[1,0]) / l))

                        elif y == taille - 1 :
                            Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW - Height_map[0,taille - 2] - Water_map2[0,taille - 2]) / l))
                            Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW - Height_map[1,taille - 1] - Water_map2[1,taille - 1]) / l))

                        else :
                            Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW- Height_map[0,y-1] - Water_map2[0,y-1]) / l))
                            Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW- Height_map[0,y+1] - Water_map2[0,y-1]) / l))
                            Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW- Height_map[1,y] - Water_map2[1,y]) / l))

                                                    
                    elif x == taille - 1 :
                        if y == 0 :
                            Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW - Height_map[x,1] - Water_map2[x,1]) / l))
                            Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,0] - Water_map2[x-1,0]) / l))
                            

                        elif y == taille - 1 :
                            Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW - Height_map[x,y-1] - Water_map2[x,y-1]) / l))
                            Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,y] - Water_map2[x-1,y]) / l))

                        else :
                            Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW - Height_map[x,y-1] - Water_map2[x,y-1]) / l))
                            Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW - Height_map[x,y+1] - Water_map2[x,y+1]) / l))
                            Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,y] - Water_map2[x-1,y])/ l))
                            
                    elif y == 0 :
                        Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW - Height_map[x,y+1] - Water_map2[x,y+1]) / l))
                        Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,y] - Water_map2[x-1,y]) / l))
                        Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW - Height_map[x+1,y] - Water_map2[x+1,y]) / l))
                       
                    elif y == taille - 1:
                        Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW - Height_map[x,y-1] - Water_map2[x,y-1]) / l))
                        Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,y] - Water_map2[x-1,y]) / l))
                        Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW - Height_map[x+1,y] - Water_map2[x+1,y]) / l))
                        
                    else : #Cas général
                        
                        Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW - Height_map[x,y-1] - Water_map2[x,y-1]) / l))
                        Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW - Height_map[x,y+1] - Water_map2[x,y+1]) / l))
                        Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,y] - Water_map2[x-1,y]) / l))
                        Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW - Height_map[x+1,y] - Water_map2[x+1,y]) / l))

                        #Calcul de l'angle local d'inclinaison , utile apres
                        
                        dhx = (Height_map[x,y+1] - Height_map[x,y-1])/2
                        dhy = (Height_map[x-1,y] - Height_map[x+1,y])/2
                        angle[x,y] = np.sqrt(dhx ** 2 + dhy ** 2) / np.sqrt( 1 + dhx ** 2 + dhy ** 2 )

                    #On ne peut pas enlever plus d'eau dans la case qu'il y en a
                    if sum(Flux) != 0:
                        K = min(1,(Water_map2[x,y] * lxy**2 )/(sum(Flux) * dt))
                    else:
                        K = 1
                    Flux = [K*Flux[0],K*Flux[1],K*Flux[2],K*Flux[3]]

                    Flux_map[x,y] = Flux
                    
        #Calcul nouveau niveau d'eau à partir des vecteurs flux
        '''deltaV = dt*(somme flux_entrant - somme flux_sortant)'''
        for x in range(taille):
                for y in range(taille):
                    if x == 0 :
                        if y == 0 :
                            deltaV = dt*(Flux_map[x,y+1,0] + Flux_map[x+1,y,2] - sum(Flux_map[x,y]))
                            deltaWX = Flux_map[x,y,1] - Flux_map[x,y+1,0]
                            deltaWY = Flux_map[x,y,3] - Flux_map[x+1,y,2]
                            
                        elif y == taille - 1 :
                            deltaV = dt*(Flux_map[x,y-1,1] + Flux_map[x+1,y,2] - sum(Flux_map[x,y]))
                            deltaWX = Flux_map[x,y-1,1] - Flux_map[x,y,0]
                            deltaWY = Flux_map[x,y,3] - Flux_map[x+1,y,2]
                            
                        else :
                            deltaV = dt*(Flux_map[x,y-1,1] + Flux_map[x+1,y,2] + Flux_map[x,y+1,0] - sum(Flux_map[x,y]))
                            deltaWX = (Flux_map[x,y-1,1] - Flux_map[x,y,0] + Flux_map[x,y,1] - Flux_map[x,y+1,0]) / 2
                            deltaWY = Flux_map[x,y,3] - Flux_map[x+1,y,2]
                            
                                                    
                    elif x == taille - 1 :
                        if y == 0 :
                            deltaV = dt*(Flux_map[x,y+1,0] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                            deltaWX = Flux_map[x,y,1] - Flux_map[x,y+1,0]
                            deltaWY = Flux_map[x-1,y,3] - Flux_map[x,y,2]

                        elif y == taille - 1 :
                            deltaV = dt*(Flux_map[x,y-1,1] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                            deltaWX = Flux_map[x,y-1,1] - Flux_map[x,y,0]
                            deltaWY = Flux_map[x-1,y,3] - Flux_map[x,y,2]
                            
                        else :
                            deltaV = dt*(Flux_map[x,y-1,1] + Flux_map[x,y+1,0] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                            deltaWX = (Flux_map[x,y-1,1] - Flux_map[x,y,0] + Flux_map[x,y,1] - Flux_map[x,y+1,0]) / 2
                            deltaWY = Flux_map[x-1,y,3] - Flux_map[x,y,2]
                            
                            
                            
                    elif y == 0 :
                            deltaV = dt*(Flux_map[x,y+1,0] + Flux_map[x+1,y,2] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                            deltaWX = Flux_map[x,y,1] - Flux_map[x,y+1,0]
                            deltaWY = (Flux_map[x-1,y,3] - Flux_map[x,y,2] + Flux_map[x,y,3] - Flux_map[x+1,y,2])/2
                        
                       
                    elif y == taille - 1:
                            deltaV = dt * (Flux_map[x,y-1,1] + Flux_map[x+1,y,2] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                            deltaWX = Flux_map[x,y-1,1] - Flux_map[x,y,0]
                            deltaWY = (Flux_map[x-1,y,3] - Flux_map[x,y,2] + Flux_map[x,y,3] - Flux_map[x+1,y,2])/2
                        
                    else :
                            deltaV = dt * (Flux_map[x,y+1,0] +Flux_map[x,y-1,1] + Flux_map[x+1,y,2] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                            deltaWX = (Flux_map[x,y-1,1] - Flux_map[x,y,0] + Flux_map[x,y,1] - Flux_map[x,y+1,0]) / 2
                            deltaWY = (Flux_map[x-1,y,3] - Flux_map[x,y,2] + Flux_map[x,y,3] - Flux_map[x+1,y,2])/2

                    
                    Water_map3[x,y] = Water_map2[x,y] + deltaV/(lxy**2)
                    d = (Water_map[x,y] + Water_map2[x,y]) / 2
                    
                    if (lxy*d)!= 0:
                        u = deltaWX / (lxy * d)
                        v = deltaWY / (lxy * d)
                    else :
                        u = 0
                        v = 0
                    
                    Velocity_field[x,y] = [u,v]
                    Capacite[x,y] = Kc * angle[x,y] * norme(Velocity_field[x,y])
                    
                    if Capacite[x,y] > Sediment_map[x,y] : 
                        Height_map[x,y] -= Ks * (Capacite[x,y] - Sediment_map[x,y])
                        s1[x,y] = Sediment_map[x,y] + Ks * (Capacite[x,y] - Sediment_map[x,y])
                    else:
                        Height_map[x,y] += Kd * (Sediment_map[x,y] - Capacite[x,y])
                        s1[x,y] = Sediment_map[x,y] - Kd * (Sediment_map[x,y] - Capacite[x,y])


        #Equation d'advection --> actualisation de la carte sediment

        f = interpolate.interp2d(X,Y,s1,kind='linear')
        for x in range(taille):
                for y in range(taille):
                    u,v = Velocity_field[x,y][0] ,Velocity_field[x,y][1]
                    a = x - u * dt
                    b = y - v * dt

                    Sediment_map[x,y] = f(a,b)[0]
        Water_map[x,y] = Water_map3[x,y] * (1 - Ke * dt)

        # Affichage de la barre de progression
        if k > prog:
            avance_barre()
            prog += n / 50
                        
    fin_barre()
    return Height_map
                
            

    

def debut_barre():
    """Affiche
╔══════════════════════════════════════════════════╗
║"""
    print(chr(9556) + chr(9552)*50 + chr(9559) + '\n' + chr(9553), end='')



def avance_barre():
    """Affiche un ▒ à la suite"""
    print(chr(9618), end='')



def fin_barre():
    """Affiche
                                                   ║
╚══════════════════════════════════════════════════╝
"""  


def Erosion_t_dt_2(taille,carte,A,l,eau,dt,Kc,Ks,Kd,Ke,n):#complexité O(n*taille*taille)
    '''algorithme d'érosion calculant la carte de hauteur apres érosion '''
    
    #l,r,t,b
    Height_map = carte
    Water_map = Water_Map_func(taille)
    Water_map2 = Water_Map_func(taille)
    Water_map3 = Water_Map_func(taille)
    Sediment_map = Sediment_Map_func(taille) 
    Flux_map = Flux_Map_func(taille)
    Velocity_field = np.array([[[1.0,1.0]]*taille]*taille)
    angle = np.zeros((taille,taille))
    Capacite = np.zeros((taille,taille))
    s1 = Sediment_Map_func(taille)


    debut_barre()
    prog = 0 # pour la barre de progression
    
    for k in range(n):

        
        Water_map2 = pluie(Water_map,eau)
        for x in range(1,taille-1):
                for y in range(1,taille-1):
                    Flux = np.array([1.0,1.0,1.0,1.0])
                    HW = Height_map[x,y] + Water_map2[x,y]

                    #On calcule les vecteurs flux a partir des cases du haut, du bas , des cotes donc on differencie juste les cas ou une des ces 4 cases n'existe 
                        
                    Flux[0] = max(0,Flux[0] + dt*A*(9.81*(HW - Height_map[x,y-1] - Water_map2[x,y-1]) / l))
                    Flux[1] = max(0,Flux[1] + dt*A*(9.81*(HW - Height_map[x,y+1] - Water_map2[x,y+1]) / l))
                    Flux[2] = max(0,Flux[2] + dt*A*(9.81*(HW - Height_map[x-1,y] - Water_map2[x-1,y]) / l))
                    Flux[3] = max(0,Flux[3] + dt*A*(9.81*(HW - Height_map[x+1,y] - Water_map2[x+1,y]) / l))

                    #Calcul de l'angle local d'inclinaison , utile apres
                        
                    dhx = (Height_map[x,y+1] - Height_map[x,y-1])/2
                    dhy = (Height_map[x-1,y] - Height_map[x+1,y])/2
                    angle[x,y] = np.sqrt(dhx ** 2 + dhy ** 2) / np.sqrt( 1 + dhx ** 2 + dhy ** 2 )

                    #On ne peut pas enlever plus d'eau dans la case qu'il y en a
                    if sum(Flux) != 0:
                        K = min(1,(Water_map2[x,y] * l**2 )/(sum(Flux) * dt))
                    else:
                        K = 1
                    Flux = [K*Flux[0],K*Flux[1],K*Flux[2],K*Flux[3]]

                    Flux_map[x,y] = Flux
                    
        #Calcul nouveau niveau d'eau à partir des vecteurs flux
        '''deltaV = dt*(somme flux_entrant - somme flux_sortant)'''
        for x in range(1,taille-1):
                for y in range(1,taille-1):
                    
                    deltaV = dt * (Flux_map[x,y+1,0] +Flux_map[x,y-1,1] + Flux_map[x+1,y,2] + Flux_map[x-1,y,3] - sum(Flux_map[x,y]))
                    deltaWX = (Flux_map[x,y-1,1] - Flux_map[x,y,0] + Flux_map[x,y,1] - Flux_map[x,y+1,0]) / 2
                    deltaWY = (Flux_map[x-1,y,3] - Flux_map[x,y,2] + Flux_map[x,y,3] - Flux_map[x+1,y,2])/2

                    
                    Water_map3[x,y] = Water_map2[x,y] + deltaV/(l**2)
                    d = (Water_map[x,y] + Water_map2[x,y]) / 2
                    
                    if (l*d)!= 0:
                        u = deltaWX / (l * d)
                        v = deltaWY / (l * d)
                    else :
                        u = 0
                        v = 0
                    
                    Velocity_field[x,y] = [u,v]
                    Capacite[x,y] = Kc * angle[x,y] * norme(Velocity_field[x,y])
                    #print(Capacite[x,y] - Sediment_map[x,y])
                    if Capacite[x,y] > Sediment_map[x,y] : 
                        Height_map[x,y] -= Ks * (Capacite[x,y] - Sediment_map[x,y])
                        s1[x,y] = Sediment_map[x,y] + Ks * (Capacite[x,y] - Sediment_map[x,y])   
                    else:
                        Height_map[x,y] += Kd * (Sediment_map[x,y] - Capacite[x,y])
                        s1[x,y] = Sediment_map[x,y] - Kd * (Sediment_map[x,y] - Capacite[x,y])


        #Equation d'advection --> actualisation de la carte sediment
        for x in range(taille):
                for y in range(taille):
                    u,v = Velocity_field[x,y][0] ,Velocity_field[x,y][1]
                    a = int(np.around(x - u * dt))
                    b = int(np.around(y - v * dt))
                    if a <0 :
                        a = 0
                    if b < 0 :
                        b = 0
                    if a >= taille :
                        a = taille - 1
                    if b >= taille:
                        b = taille - 1

                    Sediment_map[x,y] = s1[a,b]

        Water_map[x,y] = Water_map3[x,y] * (1 - Ke * dt)

        # Affichage de la barre de progression
        if k > prog:
            avance_barre()
            prog += n / 50
                        
    fin_barre()
    return Height_map               
                
def somme(mat):
    S = 0
    for x in range(len(mat)):
        for y in range(len(mat)):
            S += abs(mat[x,y])
    return S
                
                

    

    


                

            
                        

                
                
                    
                
                    
    
    
    
