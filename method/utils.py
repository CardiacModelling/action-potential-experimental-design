#
# Utility functions
#
import pints


def save_protocol(fn, p, mode='vc'):
    # fn: file name to save as.
    # p: parameters to save.
    if mode == 'vc':
        with open(fn, 'w') as f:
            f.write('# Voltage [mV]\tDuration [ms]\n')
            for i in range(len(p) // 2):
                f.write(pints.strfloat(p[2 * i]) \
                        + '\t' \
                        + pints.strfloat(p[2 * i + 1]) \
                        + '\n' \
                        )
    elif mode == 'cc':
        with open(fn, 'w') as f:
            f.write('# Holding duration [ms]\n')
            for i in range(len(p)):
                f.write(pints.strfloat(p[i]) \
                        + '\n' \
                        )
    else:
        raise ValueError('mode must be either `vc` or `cc`.')
