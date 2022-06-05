import numpy as np
from numpy import linalg
from scipy import optimize
def fit_circle(x, y):
    # == METHOD 1 ==
    method_1 = 'algebraic'

    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)

    # calculation of the reduced coordinates
    u = x - x_m
    v = y - y_m

    # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv  = sum(u*v)
    Suu  = sum(u**2)
    Svv  = sum(v**2)
    Suuv = sum(u**2 * v)
    Suvv = sum(u * v**2)
    Suuu = sum(u**3)
    Svvv = sum(v**3)

    # Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
    uc, vc = linalg.solve(A, B)

    xc_1 = x_m + uc
    yc_1 = y_m + vc

    # Calcul des distances au centre (xc_1, yc_1)
    Ri_1     = np.sqrt((x-xc_1)**2 + (y-yc_1)**2)
    R_1      = np.mean(Ri_1)
    residu_1 = np.sum((Ri_1-R_1)**2)

    return xc_1, yc_1, R_1

def fit_circle2(x, y):

    #  == METHOD 2 ==
    # Basic usage of optimize.leastsq

    method_2  = "leastsq"

    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((x-xc)**2 + (y-yc)**2)

    def f_2(c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    x_m = np.mean(x)
    y_m = np.mean(x)

    center_estimate = x_m, y_m
    center_2, ier = optimize.leastsq(f_2, center_estimate)

    xc_2, yc_2 = center_2
    Ri_2       = calc_R(xc_2, yc_2)
    R_2        = Ri_2.mean()
    residu_2   = np.sum((Ri_2 - R_2)**2)
    residu2_2  = np.sum((Ri_2**2-R_2**2)**2)
    #ncalls_2   = f_2.ncalls

    return xc_2, yc_2, R_2

def create_circle(xc, yc, r):
    circle_x = []
    circle_y = []
    for i in range(0, 361):
        circle_x.append(xc + r * np.sin(np.radians(i)))
        circle_y.append(yc + r * np.cos(np.radians(i)))
    return np.array(circle_x), np.array(circle_y)

