import numpy as np
import matplotlib.pyplot as plot
#imput parameters
x_mean=5.
standard_deviation=3.

# Bundled Provided Code into a function, Save as an Image once distribution is made
def GenerateNormalDistribution (sampleSize):
    
    # numpy function to sample the normal distribution
    x_norm=np.random.normal(x_mean,standard_deviation,sampleSize)
    #
    f=open('x_normal.dat','w')
    for x in x_norm:
        f.write(str(x)+'\n')
    f.close()
    # histogram of x
    nbins=100
    (h,x)=np.histogram(x_norm,nbins)
    fhisto=open('x_histo.dat','w')
    for i in range(len(h)):
        fhisto.write(str(x[i])+' '+str(h[i])+' '+'\n')
    fhisto.close()
    
    #Change in |<x> - x-mean|^2
    average = np.average(x_norm)
    variance = np.var(x_norm)
    std = np.sqrt(variance)
    
    # Get the Differences
    meanDiff = abs(average - x_mean)**2
    varianceDiff = abs(variance - standard_deviation**2)**2
    stdDiff = abs(std - standard_deviation)**2
    
    # Textify it
    mean_text = f'Mean = {meanDiff:.3e}'
    std_text = f'STD = {stdDiff:.3e}'
    variance_text = f'Variance = {varianceDiff:.3e}'

    #Save Plot as Image
    plot.hist(x_norm, bins=nbins)
    plot.title(f"Normal Distribution (Sample Size = {sampleSize:.3e})")
    
    # Add Text to Plot
    plot.text(0.95, 0.95, mean_text, transform=plot.gca().transAxes, fontsize=12,
          verticalalignment='top', horizontalalignment='right')
    plot.text(0.95, 0.90, std_text, transform=plot.gca().transAxes, fontsize=12,
          verticalalignment='top', horizontalalignment='right')
    plot.text(0.95, 0.85, variance_text, transform=plot.gca().transAxes, fontsize=12,
          verticalalignment='top', horizontalalignment='right')
    
    # Save Image
    plot.savefig(f"NormalDist_SampleSize_{sampleSize}.png")
    plot.close()
    
# Start with Sample Size of 10000
sample_size = 10000

# Loop 10 times and multiply the Sample Size by 2
for i in range(10):
    GenerateNormalDistribution(sample_size)
    sample_size*=2
    
