import numpy as np
import scipy.integrate as sc


# -32/16 | (−0.02135x4−0.7434x3+ 0.1465x2+ 7.867x+ 5.721)dx =
# = -0.00427x5 - 0.18585x4 + 0.0488333x3 + 3.9335x2 + 5.721x + C |-32 ~~ 16

def trapezcomp(f, a, b, n):

    # f - funckja do policzenia
    # a- dolna granica
    # b- górna granica
    # n - dokładność 

    # Inicjalizacja
    h = (b - a) / n
    x = a

    In = f(a)
    for k in range(1, n):
        x  = x + h
        In += 2*f(x)

    return (In + f(b))*h*0.5

def romberg(f, a, b, p):

    # f - funckja do policzenia
    # a - dolna granica
    # b - górna granica
    # p - dokładność liczenia romberga

    I = np.zeros((p, p))
    for k in range(0, p):
        # Composite trapezoidal rule for 2^k panels
        I[k, 0] = trapezcomp(f, a, b, 2**k)

        # Romberg rekursja
        for j in range(0, k):
            I[k, j+1] = (4**(j+1) * I[k, j] - I[k-1, j]) / (4**(j+1) - 1)

    return I


# Podpunkt a
def policzwiel(x):
    return -0.00427*x**5 - 0.18585*x**4 + 0.0488333*x**3 + 3.9335*x**2 + 5.721*x

def calka(x):
    return -0.02135*x**4 - 0.7434*x**3 + 0.1465*x**2 + 7.867*x + 5.721    

wynik_anal = policzwiel(16) - policzwiel(-32)    

wynik_romberg = romberg(calka,-32,16,3)


wyniki_gauss, err = sc.quadrature(calka, -32, 16)


print(f'Analitycznie: \t \t \t \t \t \t \t \t \t{wynik_anal}')
print('-----------------------------------------------------------------------------------------------------------')
print('Romberg:', wynik_romberg[0,0], ' '.join(map(str,wynik_romberg[1,0:2])), ' '.join(map(str, wynik_romberg[2])))
print('-----------------------------------------------------------------------------------------------------------')
print(f'Gauss: \t \t \t \t \t \t \t \t \t \t{wynik_anal}')
