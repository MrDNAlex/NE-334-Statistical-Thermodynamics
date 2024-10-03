import numpy as np
sample_size=10000
x_uniform=np.random.random(sample_size)
f=open('x_uniform.dat','w')
for x in x_uniform:
	print(x)
	f.write(str(x)+'\n')
f.close()
