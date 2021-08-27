#
# Utility functions
#
import pints


def save_protocol(fn, p):
    # fn: file name to save as.
    # p: parameters to save.
    with open(fn, 'w') as f:
        f.write('# Voltage [mV]\tDuration [ms]\n')
        for i in range(len(p) // 2):
            f.write(pints.strfloat(p[2 * i]) \
                    + '\t' \
                    + pints.strfloat(p[2 * i + 1]) \
                    + '\n' \
                    )
