from math import acos,exp, pi

def prandtl_tip_loss_factor(r, R, B, a, tsr):
    mu = r/R
    return (2/pi) * acos(exp((mu - 1)*B*tsr/(2*(1-a))))
