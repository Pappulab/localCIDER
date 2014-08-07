class Residue:
    """ 
       Residue class which defines the properties for each residue.       
    """
    def __init__(self,name,letterCode3,letterCode1,hydropathy,charge):

        
        self.name = name                      # full name     [glycine]
        self.letterCode3 = letterCode3        # 3 letter code [gly]
        self.letterCode1 = letterCode1        # 1 letter code [g]
        self.hydropathy = hydropathy          # hydropathy score
        self.charge = charge                  # charge
