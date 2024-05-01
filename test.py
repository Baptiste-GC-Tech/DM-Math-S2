import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


import numpy as np

import numpy as np

def solve_(ode_func, t_span, initial_conditions, t_eval=None, method='RK45'):
    """
    Résout numériquement une équation différentielle ordinaire (EDO) en utilisant la méthode de Runge-Kutta d'ordre 4 (RK45).

    Paramètres :
    - ode_func : La fonction décrivant le système d'EDO.
    - t_span : La plage de valeurs sur laquelle résoudre les équations.
    - initial_conditions : Les conditions initiales du système.
    - t_eval : Les points de temps où évaluer la solution (optionnel).
    - method : La méthode de résolution (par défaut, 'RK45').

    Retourne :
    - t_eval : Les valeurs de temps évaluées.
    - y_results : Les résultats de la solution aux points de temps correspondants.

    La méthode de Runge-Kutta d'ordre 4 (RK45) est une méthode numérique pour résoudre les EDO. Elle consiste à approximer la solution à chaque étape de temps en utilisant une combinaison de pentes calculées à différents points de la plage de temps. 
    Cette méthode est d'ordre 4, ce qui signifie qu'elle atteint une précision quadratique par rapport à la taille de l'intervalle de temps, ce qui en fait une méthode assez précise pour de nombreuses applications.
    """

    # Déballage des arguments
    t0, tf = t_span
    y0 = initial_conditions
    
    # Vérification de la méthode de résolution
    if method != 'RK45':
        raise ValueError("Méthode de résolution non supportée. Seule la méthode 'RK45' est implémentée.")
    
    # Préparation de l'évaluation du temps
    if t_eval is None:
        t_eval = np.linspace(t0, tf, 100)
    
    # Initialisation des résultats
    y_results = np.zeros((len(t_eval), len(y0)))
    y_results[0] = y0
    
    # Boucle principale de résolution
    for i in range(1, len(t_eval)):
        t = t_eval[i]
        y_prev = y_results[i - 1]
        
        # Calcul de la pente en utilisant la fonction fournie
        slope = ode_func(t, y_prev)
        
        # Mise à jour des valeurs en utilisant la méthode de RK45
        k1 = slope
        k2 = ode_func(t + 0.5 * (t_eval[i] - t_eval[i - 1]), y_prev + 0.5 * (t_eval[i] - t_eval[i - 1]) * k1)
        k3 = ode_func(t + 0.5 * (t_eval[i] - t_eval[i - 1]), y_prev + 0.5 * (t_eval[i] - t_eval[i - 1]) * k2)
        k4 = ode_func(t + (t_eval[i] - t_eval[i - 1]), y_prev + (t_eval[i] - t_eval[i - 1]) * k3)
        y_results[i] = y_prev + (t_eval[i] - t_eval[i - 1]) * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    return t_eval, y_results




def trajectoire_nautilus(duree, dt):

    # Paramètres physiques
    mass = 10000000  # Masse du Nautilus
    g = 9.81  # Accélération gravitationnelle e
    rho_eau = 1000  # Densité de l'eau
    volume = 100  # Volume du Nautilus 
    rho_objet = 500  # Densité de l'objet
    rayon = 2  # Rayon du Nautilus en m
    coefficient_traînée = 0.5  # Coefficient de traînée
    surface_frontale = np.pi * rayon ** 2  # Surface avant du Nautilus

    # Conditions initiales
    position = np.array([0, 0]) 
    vitesse = np.array([0, 0])  

    def f(t, y):
        position_x,position_y, vitesse_x, vitesse_y = y

        # Calcul de la densité effective de l'eau en fonction de la profondeur
        if position_y <= 0:
            rho_effectif = rho_eau
        else:
            rho_effectif = rho_eau + 0.1 * position_y  # Variation linéaire de la densité avec la profondeur

        # Calcul des forces
        poids = -mass * g * np.array([0, 1])  # Poids dirigé vers le bas
        poussée_archimède = rho_effectif * volume * g * np.array([0, 1])  # Poussée d'Archimède dirigée vers le haut
        frottements = -0.5 * rho_effectif * np.linalg.norm([vitesse_x, vitesse_y]) * np.array([vitesse_x, vitesse_y]) * surface_frontale * coefficient_traînée  # Frottements avec l'eau
        perturbation_force = np.random.normal(loc=0, scale=10, size=2)  # Perturbation aléatoire des frottements
        frottements += perturbation_force

        # Calcul de l'accélération
        acceleration_x = (poids[0] + poussée_archimède[0] + frottements[0]) / mass
        acceleration_y = (poids[1] + poussée_archimède[1] + frottements[1]) / mass

        return np.array([vitesse_x, vitesse_y, acceleration_x, acceleration_y])

    # Définition des instants de temps
    t_span = (0, duree)
    t_eval = np.arange(0, duree, dt)

    # Résolution des équations différentielles avec solve_ivp
    t_eval, solution = solve_(f, t_span, [position[0], position[1], vitesse[0], vitesse[1]], t_eval=t_eval)

    # Extraction des résultats
    positions_x = solution[:,0]
    positions_y = solution[:,1]

    # Tracer la trajectoire

    plt.plot(positions_x, positions_y)
    plt.xlabel('Position en x (km)')
    plt.ylabel('Position en y (m)')
    plt.title('Trajectoire du Nautilus')
    plt.grid(True)
    plt.show()

# Appel de la fonction avec une durée de 100 secondes et un pas de temps de 0.1 seconde
trajectoire_nautilus(1000, 0.1)
