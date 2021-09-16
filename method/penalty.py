#
# OED (penalty) score function to reduce the duration of the protocol.
#
import numpy as np
import pyoed


class DurationPenalty(pyoed.DesignMeasure):
    """
    A penalty for long ('unnecessary') duration steps and protocols.

    This tries to minimise the duration of the steps. It takes the form:

        d = b * ||duration||_1

    where ||.||_1 is the L1-norm of the vector, and b is the proportionality
    constant. The choice of its size (i.e. the proportionality constant b)
    depends on the design measure that it is added to. Perhaps a reasonable
    choice would be 10% of the best value.
    """
    
    def __init__(self, n, b):
        """
        # n: Number of parameters.
        # b: The proportionality constant of the penalty function.
        """
        super(DurationPenalty, self).__init__()
        self.n = int(n)
        self.b = float(b)

    def n_parameters(self):
        return self.n

    def __call__(self, param):
        # Return the penalty value.

        # Assuming param = [v1, t1, v2, t2, ...]
        e = np.sum(np.abs(param[1::2]))
        return self.b * e
