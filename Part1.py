import matplotlib.pyplot as plt
import numpy as np
import math

abscisse = [2, 6, 4, 8, 8, 9, 9, 6, 6, 2]
ordonne = [7, 7, 9, 9, 7, 5, 3.2, 3, 2, 3]

abscisseTangeante = [1, 3, 5.5, 6.5, 3, 5, 6.5, 10, 7, 9, 8, 10, 8, 10, 5, 7, 5.6, 7.4, 1.6, 2.2]
ordonneTangeante = [6, 8, 6, 8, 10, 8, 10, 8, 6, 8, 4, 6, 3.2, 3.4, 2, 4, 3, 1, 5, 1]

def tangente(at,ot):
    coefDirecteur = []
    for i in range (0,len(at),2):
        coefDirecteur.append((ot[i+1]-ot[i])/(at[i+1]-at[i]))
    return coefDirecteur
        
def needRotation(abscisse1, abscisse2):
    if (abscisse2 < abscisse1):
        return True
    else:
        return False
    
def rotation(_teta):
    ## a modifier   <-- Condition need to add
    n = 1000
    min = -3
    max = 3
    h= (max-min) /n

    x = []
    y = []

    xprim = []
    yprim = []

    for i in range (n):
        
        x.append(h * i + min)
        y.append((x[i]**2) + 1)

        xprim.append(x[i]*math.cos(_teta) - math.sin(_teta) * y[i])
        yprim.append(x[i]*math.sin(_teta) + math.cos(_teta) * y[i])
   
    plt.plot(xprim, yprim)

    plt.show()

# Calculates the coef dir for a point on a curve.
# 'x' is the abs, and 'fx' is the ord.
# (x, fx) represent the point we are doing this on. Same logic with _prev and _next.
def fDERIV(x, fx, x_prev, fx_prev, x_next, fx_next):  #TODO Make sure that there is no need to handle cases where we can't do a "centered" derivative
    Fi_forward = (fx_next - fx) / (x_next - x)
    Fi_backward = (fx - fx_prev) / (x - x_prev)

    Fi_centered = ((Fi_forward + Fi_backward) / 2)

    return Fi_centered

def teta(x, xi, xiplus1):
    teta = ((x - xi) / (xiplus1 - xi))
    return teta

# Hermite mathematically : H(x) = f(xi)*S1((x - xi) / (xi_plus1 - xi)) + f(xi_plus1)*S2(...) + (xi_plus1 - xi)*f'(x)*S3(...) + (xi_plus1 - xi)*f'(xi)*S4(...)
# 'xLIM_min' and 'xLIM_max' represent the interval we are interpolating. (xLIM_min, fxLIM_min) represents a point. Same logic for _max.
# 'x' varies within limits. Also, the same interval is used to calculate the derivatives in the formula.
# 'Der' represents even further points we choose so that we can use the CENTERED in the formula for more precision.
def hermites(x, xLIM_min, yLIM_min, xLIM_max, yLIM_max, xDer_min, yDer_min, xDer_max, yDer_max):
    # Preparing the tool polynoms' argument value, and calculating results for all 4 of them
    PhiArg = (x - xLIM_min) / (xLIM_max - xLIM_min)
    Phi_1 = (PhiArg - 1)**2 * (2*PhiArg + 1)
    Phi_2 = PhiArg**2 * (-2*PhiArg + 3)
    Phi_3 = (PhiArg - 1)**2 * PhiArg
    Phi_4 = PhiArg**2 * (PhiArg - 1)

    # Calculates and returns the Hermite result for the specific point x
    # This line includes every derivatives and multiplications. Read carefully if it needs debug.
    resFuncPart = (yLIM_min * Phi_1) + (yLIM_max * Phi_2)
    resDerPart1 = ((xLIM_max - xLIM_min) * fDERIV(xLIM_min, yLIM_min, xDer_min, yDer_min, xLIM_max, yLIM_max) * Phi_3)
    resDerPart2 = ((xLIM_max - xLIM_min) * fDERIV(xLIM_max, yLIM_max, xLIM_min, yLIM_min, xDer_max, yDer_max) * Phi_4)

    return resFuncPart + resDerPart1 + resDerPart2



#TODO Algorithm to actually draw the stuff
stock = tangente(abscisseTangeante, ordonneTangeante)

print(fDERIV(ordonne[3], ordonne[4], ordonne[5], abscisse[3], abscisse[4], abscisse[5]))

plt.plot(abscisseTangeante, ordonneTangeante, marker='o', linestyle='-', color='b')
plt.xlabel('Abscisse')
plt.ylabel('Dérivée estimée')
plt.title('Estimation des dérivées aux points donnés')
#plt.show()
