import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

base = 'eigenvalues-vc/eigenvalues-lsa-gsa-'
models = ['ohara-2011', 'model-list']
lsas = ['LSA-A', 'LSA-D', 'LSA-E']
gsas = ['GSA-A', 'GSA-D', 'GSA-E']

start = -0.25
increment = 0.5
space = 0.5
current = start

plt.figure(figsize=(10, 4))

'''
for i in range(3):
    x = np.loadtxt(f'{base}rand-{i}.txt')[:, 0]
    for v in x:
        plt.plot([current, current + increment], [v] * 2, c='#7f7f7f')
    current += increment + space
'''

for k in ['ch3', 'groenendaal-2015']:
    x = np.loadtxt(f'{base}{k}.txt')[:, 0]
    for v in x:
        plt.plot([current, current + increment], [v] * 2, c='#7f7f7f')
    current += increment + space


c = ['C0', 'C1']
for i, model in enumerate(models):
    for lsa in lsas:
        x = np.loadtxt(f'{base}{model}-{lsa}.txt')[:, 0]
        for v in x:
            plt.plot([current, current + increment], [v] * 2, c=c[i])
        current += increment + space

plt.yscale('log')
plt.ylabel(r'Eigenvalues of FIM$^{-1}$')
plt.xticks(np.arange(8),
           ['Lei et al.', 'Groenendaal\net al.',
            'LSA A\nSingle', 'LSA D\nSingle', 'LSA E\nSingle',
            'LSA A\nAverage', 'LSA D\nAverage', 'LSA E\nAverage'])
plt.tight_layout()

plt.show()
