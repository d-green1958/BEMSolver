
class ConvergenceError(Exception):
    """Exception raised for lack of convergence in a computation.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="Lack of convergence"):
        self.message = message
        super().__init__(self.message)
        
        
        
class InputDataError(Exception):
    """Exception raised if input data appears inccorect.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="Lack of convergence"):
        self.message = message
        super().__init__(self.message)