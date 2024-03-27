
# Using Oye model from aerodyn user manual

class DynamicStallModel:
    """
    A class used to calculate quantities relating to dynamic stall.
    """

    def __init__(self, blade):
        """
        Initialise the DynamicStallModel object with the data stored in the blade object.

        Args:
            blade (BladeGeometry): The blade geometry file with links to the aerofoil
            data through AerofoilLookups objects.
        """

        print("Adding dynamic stall model")

        from bempy.source import BladeGeometry
        self.blade = blade
        self.aerofoil_names = self.blade.aerofoil_type
        self.aerfoil_dict = self.blade.aerofoil_dict
        self.chord_lengths = self.blade.chord_lengths

        # array of angle of attacks
        self.alpha = self.__form_alpha()

        # the lift coefficient at each angle of attack
        self.C_l = self.__form_C_l()

        # the array of AOAs giving 0 lift for each element
        self.alpha_0 = self.blade.alpha_0

        # array of seperation function values for each element
        self.f_ssteady = [0]*blade.number_of_nodes
        self.f_s = [0]*blade.number_of_nodes
        self.f_s_prev = [0]*blade.number_of_nodes
        

        # array of arrays of interpolated gradient values
        self.C_lalpha = self.__form_C_lalpha()

        # the inviscide lift coefficient
        self.C_linv = self.__form_C_linv()

        # the fully seperated life coefficeint
        self.C_lfs = [0]*blade.number_of_nodes

        # scales and constants used to relax steady f_s into dynamic f_s
        self.tau = None
        self.k = 3
        self.dt = None

        # the dynamic lift coefficeint (the output of the class)
        self.C_ldyn = None
        
        

    def __form_alpha(self):
        """
        import the data from the aerofoil look up tables for use later

        Returns:
            array: array of the angles of attack values
        """
        alpha = []
        for ind, name in enumerate(self.aerofoil_names):
            element_alpha = self.aerfoil_dict[name].angle_of_attack
            alpha.append(element_alpha)

        self.alpha = alpha
        return alpha

    def __form_C_l(self):
        """
        import the data from the aerofoil look up tables for use later

        Returns:
            array: array of the lift coefficient values
        """
        C_l = []
        for ind, name in enumerate(self.aerofoil_names):
            element_C_l = self.aerfoil_dict[name].C_lift
            C_l.append(element_C_l)

        self.C_l = C_l
        return C_l

    def __form_C_lalpha(self):
        """
        Takes the aerofoil polars and computes C_lalpha, returning it and 
        passing it so self.C_lalpha
        Returns:
            C_lalpha (array): an array of arrays for the gradient in C_l 
        """
        C_lalpha = []
        for ind, name in enumerate(self.aerofoil_names):
            # first get the steady state lift coefficients
            C_l = self.C_l[ind]
            alpha = self.alpha[ind]
            alpha_0 = self.alpha_0[ind]

            # now we compute the gradient
            from numpy import gradient, interp
            gradient_arr = (gradient(C_l, alpha)).tolist()
                
            # now find value at alpha_0
            C_lalpha.append(interp(alpha_0, alpha, gradient_arr))
            

        self.C_lalpha = C_lalpha
        return C_lalpha

    def __form_C_linv(self):
        """
        Computes the inviscid life coefficint evaluated at each alpha.
        Args:
            alpha_value (float): the value of alpha we are interested in
        Returns:
            C_linv (float): the inviscid lift coefficient
        """
        C_linv = []
        for ind, name in enumerate(self.aerofoil_names):
            alpha = self.alpha[ind]
            alpha_0 = self.alpha_0[ind]
            C_lalpha = self.C_lalpha[ind]

            from numpy import sin, deg2rad
            # C_linv.append(C_lalpha * sin(deg2rad(alpha - alpha_0)))
            C_linv.append(C_lalpha * (alpha - alpha_0))
            # documentation also uses a small angle approximation on the sin

        self.C_linv = C_linv
        return C_linv
    
    def interp_C_linv(self, alpha_value, element_num):
        """
        linearly interpolates the values of the inviscid lift coefficient to some value
        given by alpha_value
        Args:
            alpha_value (float): evaluation point
            element_num (int): the blade element number
        """
        from numpy import interp
        return interp(alpha_value, self.alpha[element_num], self.C_linv[element_num])
    
    def interp_C_l(self, alpha_value, element_num):
        """
        linearly interpolates the values of the lift coefficient to some value
        given by alpha_value
        Args:
            alpha_value (float): evaluation point
            element_num (int): the blade element number
        """
        from numpy import interp
        return interp(alpha_value, self.alpha[element_num], self.C_l[element_num])


    def __form_f_ssteady(self):
        """
        calculate the steady state seperation function value for some given AOA
        at a given blade element.
        Args:
            alpha_value (float): the angle of attack
            element_num (int): the blade element number
            C_l (float): the steady state lift coefficient
        """
        f_ssteady = []
        from numpy import sqrt, isnan
        for element_num, alpha_arr in enumerate(self.alpha):
            f_s_arr = []
            for alpha_value in alpha_arr:
                numer = self.interp_C_l(alpha_value, element_num)
                denom = (alpha_value - self.alpha_0[element_num]) * self.C_lalpha[element_num]
                
                print(denom, numer)
                
                if numer/denom <0:
                    raise ValueError("ASDOASHDOASDUH")
                    # possible errors at extreme angles due to signs of numerator and denom!!!
                
                temp = (2 * sqrt(numer/denom) - 1)**2
                temp = min(temp,1)
                f_s_arr.append(temp)
            
            f_ssteady.append(f_s_arr)
                
        self.f_ssteady = f_ssteady
        return f_ssteady
    
    def calculate_f_ssteady(self, alpha_value, element_num):
        """
        Calculate the steady state seperation function given some bllade element and the 
        angle of attack. Note: this is not done in advance due to issues with divergence
        at angles of attack that are unlikely to occur (extreme angles of attack).

        Args:
            alpha_value (float): The AOA
            element_num (int): Blade element index
            
        Returns:
            f_ssteady (float): the steady state seperation function
        """
        
        C_l = self.interp_C_l(alpha_value, element_num)
        numer = C_l
        denom = (alpha_value - self.alpha_0[element_num]) * self.C_lalpha[element_num]
        
        if denom == 0: # then we have alpha = alpha_0 or C_lalpha = 0
            #  this is likely the case for cylindrical blade types and so we assume f_ssteady = 1
            return 1
        
        temp = numer/denom
        
        f_ssteady = (2*(temp)**0.5 - 1)**2
        from numpy import isnan
        if isnan(f_ssteady):
            err_str = f"f_ssteady has diverged! the numerator is {numer} and the denominator is {denom}"
            raise ValueError(err_str)
        
        f_ssteady = min(f_ssteady, 1)
        self.f_ssteady[element_num] = f_ssteady
        
        return f_ssteady
    
    def calculate_C_lfs(self, alpha_value, element_num):
        """
        Function to calculate the fully seperated lift coefficient
        Args:
            alpha_value (float): AOA value
            element_num (int): blade element index

        Returns:
            float: the fully seperated lift coefficient
        """
        
        C_l = self.interp_C_l(alpha_value, element_num)
        C_lalpha = self.C_lalpha[element_num]
        alpha_0 = self.alpha_0[element_num]
        f_ssteady = self.f_ssteady[element_num]
        
        if f_ssteady != 1:
            numer = C_l - C_lalpha*(alpha_value - alpha_0) * f_ssteady
            denom = 1- f_ssteady
            temp = numer/denom
            
        else:
            temp = C_l/2
            
        
        self.C_lfs[element_num] = temp
        return temp
        
        
    
    def relax_f_ssteady(self, element_num, V_rel):
        """
        Function used to relax f_ssteady into f_s which is then used for calculations
        of the dynamic lift coefficient.
        Args:
            element_num (int): the blade element
            V_rel (float)): the relative velocity at the blade element
        Returns:
            f_s (float): the relaxed seperation function
        """
        f_ssteady = self.f_ssteady[element_num]
        f_s_prev = self.f_s_prev[element_num]
        dt = self.dt
        tau = self.k * self.chord_lengths[element_num] / V_rel
                
        from math import exp
        f_s = f_ssteady + (f_s_prev - f_ssteady) * exp(-dt/tau)
        
        self.f_s[element_num] = f_s
        return f_s
    
    def calculate_C_ldyn(self, alpha_value, element_num):
        """
        function to calculate the dynamic lift coeff. 

        Args:
            alpha_value (_type_): _description_
            element_num (_type_): _description_

        Returns:
            _type_: _description_
        """
        f_s = self.f_s[element_num]
        C_linv = self.interp_C_linv(alpha_value, element_num)
        C_lfs = self.calculate_C_lfs(alpha_value, element_num)
        
        C_ldyn = f_s * C_linv + (1-f_s)*C_lfs
        self.C_ldyn = C_ldyn
        return C_ldyn
    
    
    
    def get_C_ldyn(self, alpha_value, element_num, V_rel, dt):
        """
        Function to get the dynamic lift coefficient given the blade element,
        angle of attack and the relative velocity. This function calls many of the 
        other member functions and is intended to be the only function called outside 
        of the class.

        Args:
            alpha_value (flaot): angle of attack
            element_num (int): blade element index
            V_rel (float): relative velocity at blade element

        Returns:
            float: dynamic lift coeff
        """
        # set the time step
        self.dt = dt
        
        # calculate the steady state seperation function
        f_ssteady = self.calculate_f_ssteady(alpha_value, element_num)
        
        # relax f_ssteady into f_s
        f_s = self.relax_f_ssteady(element_num, V_rel)
        
        # calculate the fully seperated lift coeff
        C_lfs = self.calculate_C_lfs(alpha_value, element_num)
        
        # and finally return the dynamic lift coefficient
        C_ldyn = self.calculate_C_ldyn(alpha_value, element_num)
        return C_ldyn
        
        
    

    
    
        
        
        
        
    
    
    
    

