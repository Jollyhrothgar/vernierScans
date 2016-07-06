#!/usr/bin/env python
import numpy as np
from math import log

def get_cross_section(s):
    """
    s is the sqrt-s beam collision energy. I.e. 510 GeV, 200 GeV. Returns the
    cross-section in mb
    """
    xsec = 35.45 + 0.308*(log(s/28.94))**2. + 42.53*(1./s)**0.458-33.34*(1./s)**0.545
    return xsec

energies = np.arange(200,520,10)

for s in energies:
    print s,get_cross_section(s)
