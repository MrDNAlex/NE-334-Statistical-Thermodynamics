import numpy as np
sample_size=10000

def EstimatePi (sampleSize):
    # Generate Uniform Distribution
	x=np.random.random(sampleSize)
	y=np.random.random(sampleSize)
 
	# Determine the number of points inside the Circle
	count_circle = 0
	for i in range(sample_size):
		r=x[i]**2+y[i]**2
		if r <= 1:
			count_circle+=1
	
	#Get the Probablility of a Point being inside the Circle
	p = count_circle/sample_size
 
	# Estimate Pi
	pi_estimate = 4*p
	
	#Calculate the Standard Deviation for Binomial
	std_binomial = np.sqrt(p * (1-p))
	
	# Calculate the Standard Error
	standard_error = std_binomial/np.sqrt(sampleSize)
 
	return [pi_estimate, standard_error]

# Q3 b and c
for i in range(15):
    [pi_estimate, standard_error] = EstimatePi(sample_size)
    print(f"Sample Size = {sample_size}, PI Estimate = {pi_estimate}, Standard Error : {standard_error}")
    sample_size*=2

    