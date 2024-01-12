from math import acos,exp, pi, sin

def prandtl_tip_loss_factor_approx(r, R_tip, B, a, tsr):
    mu = r/R_tip
    return (2/pi) * acos(exp((mu - 1)*B*tsr/(2*(1-a))))

def prandtl_tip_loss_factor(r, R_tip, B, phi):
        return (2/pi) * acos(exp(-B*(R_tip-r)/(2*r*sin(phi))))
    
def xu_sankar_tip_loss(r, R_tip, B, phi):
    mu = r/R_tip
    if mu < 0.7:
        r = 0.7*R_tip
        return 1 - mu * (1 - prandtl_tip_loss_factor(r, R_tip, B, phi))/0.7
    else:
        return (prandtl_tip_loss_factor(r, R_tip, B, phi)**0.85 + 0.5)/2

def prandtl_hub_loss_factor(r, R_hub, B, phi):
    return (2/pi) * acos(exp(-B*(r-R_hub)/(2*r*sin(phi))))