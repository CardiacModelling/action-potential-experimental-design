#
# Parameter transformations.
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import numpy as np


def log_transform_from_model_param(param):
    # Apply natural log transformation to all parameters
    out = np.copy(param)
    out = np.log(out)
    return out


def log_transform_to_model_param(param):
    # Inverse of log_transform_from_model_param()
    # Apply natural exp transformation to all parameters
    out = np.copy(param)
    out = np.exp(out)
    return out


def donothing(param):
    out = np.copy(param)
    return out
