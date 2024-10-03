import numpy as np
sample_size=10000

def GenerateUniformDistribution(sampleSize):
	x_uniform=np.random.random(sampleSize)
	f=open('x_uniform.dat','w')
	for x in x_uniform:
		f.write(str(x)+'\n')
	f.close()
	return x_uniform

# Q2 b
for i in range(10):
	distribution = GenerateUniformDistribution(sample_size)
 
 	#Get Average and Display Stats
	average = np.average(distribution)
	difference = abs(average - 0.5)**2
	print(f"Sample Size : {sample_size}, Average : {average:.3f}, Difference : {difference:.3e}")
	sample_size*=2

print("\n\n\n")

#Q2 c
# Modify things to generate numbers between -1 and 3
sample_size=10000
for i in range(10):
	distribution = GenerateUniformDistribution(sample_size)
 
	#Multiply by 4 and subtract 1 to Shift things to (-1 <-> 3)
	distribution*= 4
	distribution-=1
 
	#Get Average and Display Stats
	average = np.average(distribution)
	difference = abs(average - 1)**2
	print(f"Sample Size : {sample_size}, Average : {average:.3f}, Difference : {difference:.3e}")
	sample_size*=2

