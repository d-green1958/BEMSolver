# caluclate the power coefficient from a power
def C_power(power, rho, U, A):
    return 2 * power / (rho * U**3 * A)
