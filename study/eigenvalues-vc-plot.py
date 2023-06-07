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

#'''
for i in range(1):
    x = np.loadtxt(f'{base}rand-{i}.txt')[:, 0]
    x = x ** (-1)
    for v in x:
        plt.plot([current, current + increment], [v] * 2, c='#7f7f7f')
    current += increment + space
#'''

for k in ['ch3', 'groenendaal-2015']:
    x = np.loadtxt(f'{base}{k}.txt')[:, 0]
    x = x ** (-1)
    c = '#7f7f7f'
    for v in x:
        plt.plot([current, current + increment], [v] * 2, alpha=.5, c=c)
    current += increment + space

c = ['C0', 'C1']
for i, model in enumerate(models):
    for lsa in lsas:
        x = np.loadtxt(f'{base}{model}-{lsa}.txt')[:, 0]
        x = x ** (-1)
        for v in x:
            plt.plot([current, current + increment], [v] * 2, c=c[0])
        current += increment + space

for i, model in enumerate(models):
    for gsa in gsas:
        x = np.loadtxt(f'{base}{model}-{gsa}.txt')[:, 0]
        x = x ** (-1)
        for v in x:
            plt.plot([current, current + increment], [v] * 2, c=c[1])
        current += increment + space

plt.yscale('log')
#plt.ylabel(r'Eigenvalues of FIM$^{-1}$')
plt.ylabel(r'Eigenvalues of FIM')
#plt.xticks(np.arange(9),
#           ['Random', 'Lei et al.', 'Groenendaal\net al.',
#            'LSA A\nSingle', 'LSA D\nSingle', 'LSA E\nSingle',
#            'LSA A\nAverage', 'LSA D\nAverage', 'LSA E\nAverage'])
plt.xticks(np.arange(15),
           ['Random', 'M', 'N', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
            'J', 'K', 'L'])
plt.tight_layout()

#plt.show()
plt.savefig('fig/eigenvalues-vc-lsa.pdf', format='pdf')
