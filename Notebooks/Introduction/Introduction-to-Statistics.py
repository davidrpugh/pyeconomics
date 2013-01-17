# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import random
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

# <markdowncell>

# Task 1: Plotting probability distributions using SciPy and Matplotlib
# =====================================================================
# 
# In this task you will learn how to create Python objects representing random variables using the SciPy sub-module [scipy.stats](http://docs.scipy.org/doc/scipy/reference/stats.html).  Once we have learned how to create randoms variable objects representing specific probability distributions, we can use methods defined for such objects to make plots of the [probability density function (pdf)](http://en.wikipedia.org/wiki/Probability_density_function), [cumulative distribution function (cdf)](http://en.wikipedia.org/wiki/Cumulative_distribution_function), and the [survival function (sf)](http://en.wikipedia.org/wiki/Survival_function).  Furthermore, we will also be able to generate random samples of data from our distributions which can be used in simulations.
# 
# Let's get started by looking at the contents of the scipy.stats module:

# <codecell>

stats.

# <markdowncell>

# Although all of the probability distributions are defined individually in the scipy.stats module, they can also be accessed collectively in the sub-sub-module called scipy.stats.distributions.

# <codecell>

stats.distributions.

# <markdowncell>

# Let's look at the [Normal (or Gaussian)](http://en.wikipedia.org/wiki/Normal_distribution) distribution more closely.

# <codecell>

stats.distributions.norm?

# <markdowncell>

# Using the information about syntax provided in the help menus we can create a random variable object representing a Normal distribution for specified values of $\mu$ and $\sigma$ as follows.

# <codecell>

# specify values for mu and sigma
mu = 0
sigma = 1

# create a random variable object representing the normal distribution
normal_rv = stats.distributions.norm(loc=mu, scale=sigma)

# <codecell>

# verify the object type
type(normal_rv)

# <markdowncell>

# Now that we have created a Normal distribution as random variable object, we can access various methods associated with all rv_frozen objects.

# <codecell>

normal_rv.

# <markdowncell>

# We can use some of these methods to create plots of the probability density function (pdf), the cumulative distribution function (cdf), and the survival function (sf)  of our Normal distribution.

# <codecell>

# create a grid of values at which to evaluate pdf and cdf
grid = np.linspace(-4, 4, 1000)

# create new Figure object
fig = plt.figure()

# plot the pdf
ax = fig.add_subplot(121)

ax.plot(grid, normal_rv.pdf(grid), 'r-')
ax.set_ylim(0, 0.5)

# label the axes...
xlab = ax.set_xlabel('X')
xlab.set_family('serif')
ylab = ax.set_ylabel('Probability density function, f(x)')
ylab.set_family('serif')

# ...and add a title
ax.set_title(r'$\mathrm{Normal\ pdf,}\ \mu=0,\ \sigma=1$')

# plot the cdf
ax1 = fig.add_subplot(122)

ax1.plot(grid, normal_rv.cdf(grid), 'b-')

# label the axes...
xlab = ax1.set_xlabel('X')
xlab.set_family('serif')
ylab = ax1.set_ylabel('Cumulative distribution function, F(x)')
ylab.set_family('serif')

# ...and add a title
ax1.set_title(r'$\mathrm{Normal\ cdf,}\ \mu=0,\ \sigma=1$')

# adjust the layout
fig.tight_layout()

plt.savefig('Normal pdf and cdf.png')
plt.show()

# <codecell>

fig = plt.figure()

ax = fig.add_subplot(111)

# plot the survival function
ax.plot(grid, normal_rv.sf(grid), 'g-')
ax.set_xscale('log')
ax.set_yscale('log')

# label the axes...
xlab = ax.set_xlabel('X (log scale)')
xlab.set_family('serif')
ylab = ax.set_ylabel('Survival function, 1 - F(x) (log scale)')
ylab.set_family('serif')

# ...and add a title
ax.set_title(r'$\mathrm{Survival\ Function\ of\ a\ Normal\ distribution\ with}\ \mu=0,\ \sigma=1$')

plt.savefig('Normal sf.png')
plt.show()

# <markdowncell>

# Exercise 1.1
# ------------
# 
# Create a rv_frozen object representing a [Chi-Square](http://en.wikipedia.org/wiki/Chi-squared_distribution) distribution with 3 degrees of freedom and re-generate the above plots.  You should be able to reuse almost all of the code I wrote above to complete this exercise, although you will need to create a new random variable object to represent you chosen distribution.

# <codecell>

# insert your code here!

# <codecell>

# specify degrees of freedom
df = 3

# create a random variable object representing the normal distribution
chi_square_rv = stats.distributions.chi2(df)

# check type of chi_square_rv
type(chi_square_rv)

# <codecell>

# create a grid of values at which to evaluate pdf and cdf
grid = np.linspace(0, 15, 1000)

# create new Figure object
fig = plt.figure()

# plot the pdf
ax = fig.add_subplot(211)

ax.plot(grid, chi_square_rv.pdf(grid), 'r-')
ax.axis('tight')
ax.set_ylim(0, 0.5)

# label the axes...
xlab = ax.set_xlabel('X')
xlab.set_family('serif')
ylab = ax.set_ylabel('f(x)')
ylab.set_family('serif')

# ...and add a title
ax.set_title(r'$\mathrm{Probability\ density\ function\ of\ a\ \chi^2\ r.v.\ with\ 3\ d.f.}$')

# plot the cdf
ax1 = fig.add_subplot(212)

ax1.plot(grid, chi_square_rv.cdf(grid), 'b-')
ax1.axis('tight')

# label the axes...
xlab = ax1.set_xlabel('X')
xlab.set_family('serif')
ylab = ax1.set_ylabel('F(x)')
ylab.set_family('serif')

# ...and add a title
ax1.set_title(r'$\mathrm{Distribution function\ of\ a\ \chi^2\ r.v.\ with\ 3\ d.f.}$')

# adjust the layout
fig.tight_layout()

plt.savefig('exercise 1.1a.png')
plt.show()

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot the survival function
ax.plot(grid, chi_square_rv.sf(grid), 'g-')
ax.set_xscale('log')
ax.set_yscale('log')

# label the axes...
xlab = ax.set_xlabel('X (log scale)')
xlab.set_family('serif')
ylab = ax.set_ylabel('Survival function, 1 - F(x) (log scale)')
ylab.set_family('serif')

# ...and add a title
ax.set_title(r'$\mathrm{Distribution function\ of\ a\ \chi^2\ r.v.\ with\ 3\ d.f.}$')

plt.savefig('exercise 1.1b.png')
plt.show()

# <markdowncell>

# Task 2: Generating random numbers using Python, NumPy and SciPy
# ===============================================================
# 
# An important part of any simulation is the ability to [draw random numbers](http://en.wikipedia.org/wiki/Random_number_generation). We have at least three methods available for generating random numbers using Python:
# 
# 1. We can use the rvs() method of our random variable objects defined using any probability distributions from stats.distributions.
# 2. We can use built-in Python routines from the module random that we imported above.
# 3. We can use NumPy's module numpy.random. 
# 
# The numbers generated are pseudo random in the sense that they are generated deterministically from a seed number, but are distributed in what has statistical similarities to random fashion. Python, NumPy, and SciPy all use a particular algorithm called the [Mersenne Twister](http://en.wikipedia.org/wiki/Mersenne_twister) to generate pseudorandom numbers.
# 
# A minimum requirement for any research (including your MSc theses!) is that it be [reproducible](http://reproducibleresearch.net/index.php/Main_Page). For research that uses random number generators, this means that seed values should almost aways* be set and these values should be reported. The seed is always an integer value. Any procedure to generate random numbers that starts with the same seed will generate exactly the same sequence of numbers each time it is run.
# 
# In addition to making your work more reproducible, manually setting seed values can be useful for debugging purposes. The random number seed can be set using the np.random.seed() command. If this command is not run, NumPy automatically selects a random seed (based on the time) that is different every time a random number generator is run. 
# 
# <p style="font-size:10px">*One does not always need to specify the seed.  For example, when we perform multiple runs of some simulation and average the runs together we will want each run to have be generated from a different sequence of random numbers.</p>

# <codecell>

np.random.seed(293423)

# <markdowncell>

# First, I want to cover how to generate random samples from our rv_frozen objects defined in the previous task.  To generate random samples from rv_frozen objects we use the rvs() method.

# <codecell>

# generates a single draw
normal_rv.rvs()

# <codecell>

# generates T draws
T = 10
normal_rv.rvs((10,))

# <codecell>

# generates an (T, n) array of draws
T = 5
n = 3
shape = (T, n)
normal_rv.rvs(shape)

# <markdowncell>

# In the cell below I draw a large number of observations from normal using the rvs() method, plot a histogram of the sample, and then overlay the probability density function using the pdf() method.

# <codecell>

# generate a large sample
T = 10000
sample = normal_rv.rvs((T,))

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# histogram the sample
n, bins, patches = ax.hist(sample, color='purple', bins=25, normed=True, alpha = 0.75)
ax.set_xlim(-4, 4)
ax.set_ylim(0, 0.5)

# overlay the theoretical pdf
grid = np.linspace(-4, 4, 1000)
ax.plot(grid, normal_rv.pdf(grid), 'r--', label=r'$N(\mu,\sigma)$')

# label the axes
ax.set_xlabel('X')
ax.set_ylabel('Probability density')

# add a title
ax.set_title('We can use the rvs() method to check the pdf() method!')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Normal histogram.png')
plt.show()

# <markdowncell>

# Exercise 2.1
# ------------
# 
# Regenerate the above graph using the rv_frozen object that you created in Exercise 1.4 from the previous task.

# <codecell>

# insert your code here!

# <codecell>

# generate a large sample
T = 10000
sample = chi_square_rv.rvs((T,))

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# histogram the sample
n, bins, patches = ax.hist(sample, color='purple', bins = 25, normed=True, alpha = 0.75)
ax.set_xlim(0, 15)

# overlay the theoretical pdf
grid = np.linspace(0, 15, 1000)
ax.plot(grid, chi_square_rv.pdf(grid), 'r--', label=r'$\chi^2(df)$')

# label the axes
ax.set_xlabel('X')
ax.set_ylabel('Probability density')

# add a title
ax.set_title('We can use the rvs() method to check the pdf() method!')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('exercise 2.1.png')

# <markdowncell>

# In the cell below I have provided code for drawing one sample of length T from a normal (Gaussian) distribution with μ=1.5 and σ=4.0 and plotting the sample path.  Note that I have set the seed value so that you can reproduce my plot <i>exactly</i>.

# <codecell>

# set params
T = 1000
mu = 1.5
sigma = 4.0

# create a new rv object
normal_rv = stats.distributions.norm(loc=mu, scale=sigma)

# draw from a continuous normal (Gaussian) distribution
np.random.seed(int(12345))
sample_path = normal_rv.rvs((T,))

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot the sample path
ax.plot(sample_path, 'r-')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel('Observation')
ylab.set_family('serif')

# set the title
title = ax.set_title('Random sample of N(1.5, 4.0) data')
title.set_family('serif')

plt.savefig('Random sample of Normal data.png')
plt.show()

# <markdowncell>

# Exercise 2.2
# -----------
# 
# Create a new rv_frozen object representing a [Student's t distribution](http://en.wikipedia.org/wiki/Student's_t-distribution) with 1 degree of freedom and recreate the above plot. Be sure to remember to set the seed!

# <codecell>

# insert your code here!

# <codecell>

#Set degrees of freedom
df = 11

# create a t-distribution object
students_t_rv = stats.distributions.t(df)

# <codecell>

# set the seed
np.random.seed(12)

# generate data
T = 1000
sample_path = students_t_rv.rvs((T,))

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot the sample path
ax.plot(sample_path, 'g-')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel('Observation')
ylab.set_family('serif')

# set the title
title = ax.set_title("Random sample of data from Student's t with 1 d.f.")
title.set_family('serif')

plt.savefig('exercise 2.2.png')

# <markdowncell>

# Task 3: Where we learn to simulate flipping a coin...
# =====================================================
# 
# In this task you will simulate flipping a coin.  First we need to create a random variable representing the [Bernoulli distribution](http://en.wikipedia.org/wiki/Bernoulli_distribution) from which coin flips are drawn.

# <codecell>

# create an object representing a fair coin
fair_coin = stats.distributions.bernoulli(0.5)

# <codecell>

# flip a fair coin 10 times
coin_flips = fair_coin.rvs((10,))
print coin_flips

# <codecell>

# count the fraction of heads...
fraction_heads = np.mean(coin_flips)
print fraction_heads

# <markdowncell>

# Exercise 3.1
# ------------
# 
# Modify the above code to simulate 100 flips of a "loaded" coin where the probability of "success" or heads is 0.75.  Increase the number of flips to 1000 and then 10000  Describe what happens to fraction_heads as the number of flips increases.

# <codecell>

# insert your code here!

# <codecell>

# create a loaded coin
loaded_coin = stats.distributions.bernoulli(0.75)

# <codecell>

# flip the loaded coin
T = 1000
loaded_coin_flips = loaded_coin.rvs((T,))

# count the fraction of heads...
fraction_heads = np.mean(loaded_coin_flips)
print fraction_heads

# <markdowncell>

# Now we know that flipping a fair coin will result in heads exactly 50% of the time.&nbsp; Thus after many independent flips, the [Law of of Large Numbers (LLN)](http://en.wikipedia.org/wiki/Law_of_large_numbers) suggests that the number of heads should more or less equal the expected number of heads...or does it!
# 
# A South African mathematician named John Kerrich was visiting Copenhagen in 1940 when Germany invaded Denmark. Kerrich spent the next five years in an interment camp where, to pass the time, he carried out a series of experiments in probability theory...including an experiment where he flipped a coin by hand 10,000 times. After the war Kerrich was released and published the results of many of his experiments. I have copied the table of the coin flipping results reported by Kerrich below. The first two collumns are self explanatory, the third column, <b>Difference</b>, is the difference between the observed number of heads and the expected number of heads.
# 
# <table border="0" align="center">
# <thead>
# <tr align="center">
# <td><span style="text-decoration: underline;"><strong>Tosses</strong></span></td>
# <td><span style="text-decoration: underline;"><strong>Heads</strong></span></td>
# <td><span style="text-decoration: underline;"><strong>Difference</strong></span></td>
# </tr>
# </thead>
# <tbody>
# <tr align="center">
# <td>10</td>
# <td>4</td>
# <td>-1</td>
# </tr>
# <tr align="center">
# <td>20</td>
# <td>10</td>
# <td>0</td>
# </tr>
# <tr align="center">
# <td>30</td>
# <td>17</td>
# <td>2</td>
# </tr>
# <tr align="center">
# <td>40</td>
# <td>21</td>
# <td>1</td>
# </tr>
# <tr align="center">
# <td>50</td>
# <td>25</td>
# <td>0</td>
# </tr>
# <tr align="center">
# <td>60</td>
# <td>29</td>
# <td>-1</td>
# </tr>
# <tr align="center">
# <td>70</td>
# <td>32</td>
# <td>-3</td>
# </tr>
# <tr align="center">
# <td>80</td>
# <td>35</td>
# <td>-5</td>
# </tr>
# <tr align="center">
# <td>90</td>
# <td>40</td>
# <td>-5</td>
# </tr>
# <tr align="center">
# <td>100</td>
# <td>44</td>
# <td>-6</td>
# </tr>
# <tr align="center">
# <td>200</td>
# <td>98</td>
# <td>-2</td>
# </tr>
# <tr align="center">
# <td>300</td>
# <td>146</td>
# <td>-4</td>
# </tr>
# <tr align="center">
# <td>400</td>
# <td>199</td>
# <td>-1</td>
# </tr>
# <tr align="center">
# <td>500</td>
# <td>255</td>
# <td>5</td>
# </tr>
# <tr align="center">
# <td>600</td>
# <td>312</td>
# <td>12</td>
# </tr>
# <tr align="center">
# <td>700</td>
# <td>368</td>
# <td>18</td>
# </tr>
# <tr align="center">
# <td>800</td>
# <td>413</td>
# <td>13</td>
# </tr>
# <tr align="center">
# <td>900</td>
# <td>458</td>
# <td>8</td>
# </tr>
# <tr align="center">
# <td>1000</td>
# <td>502</td>
# <td>2</td>
# </tr>
# <tr align="center">
# <td>2000</td>
# <td>1013</td>
# <td>13</td>
# </tr>
# <tr align="center">
# <td>3000</td>
# <td>1510</td>
# <td>10</td>
# </tr>
# <tr align="center">
# <td>4000</td>
# <td>2029</td>
# <td>29</td>
# </tr>
# <tr align="center">
# <td>5000</td>
# <td>2533</td>
# <td>33</td>
# </tr>
# <tr align="center">
# <td>6000</td>
# <td>3009</td>
# <td>9</td>
# </tr>
# <tr align="center">
# <td>7000</td>
# <td>3516</td>
# <td>16</td>
# </tr>
# <tr align="center">
# <td>8000</td>
# <td>4034</td>
# <td>34</td>
# </tr>
# <tr align="center">
# <td>9000</td>
# <td>4538</td>
# <td>38</td>
# </tr>
# <tr align="center">
# <td>10000</td>
# <td>5067</td>
# <td>67</td>
# </tr>
# </tbody>
# </table>

# <codecell>

# input the Kerrich data
kerrich_data = np.array([[10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000], [4, 10, 17, 21, 25, 29, 32, 35, 40, 44, 98, 146, 199, 255, 312, 368, 413, 458, 502, 1013, 1510, 2029, 2533, 3009, 3516, 4034, 4538, 5067], [-1, 0, 2, 1, 0, -1, -3, -5, -5, -6, -2, -4, -1, 5, 12, 18, 13, 8, 2, 13, 10, 29, 33, 9, 16, 34, 38, 67]], float).T

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot Kerrich's difference
ax.plot(kerrich_data[:,0], kerrich_data[:,2], 'r-', label='Kerrich Data')

# set the x-axis to have log scale
ax.set_xscale('log')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel('Observed heads - expected heads')

# set the title
title = ax.set_title('Divergence between observed heads and expected heads?')
title.set_family('serif')

# add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('Kerrich difference.png')
plt.show()

# <markdowncell>

# WTF!  This above plot seems to imply that the number of observed heads is diverging from the expected number of heads as we increase the number of flips (which is the exact opposite of our intuition)!  Perhaps Kerrich made a mistake.  Fortunately, we can check his results via simulation! 

# <codecell>

# set params
N = 100
T = 10000

# generate an array of integers (each of n columns can be interpreted as a run of length T)
data = fair_coin.rvs((T,N))

# create an array in which to store our sample averages
difference = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        difference[i, j] = 2 * np.sum(data[:i + 1, j]) - (i + 1)

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot each sample path...
for i in range(N):
    ax.plot(np.arange(1, T + 1, 1), difference[:, i], 'k-', alpha=0.05)

# plot Kerrich's observed difference
ax.plot(kerrich_data[:,0], kerrich_data[:,2], 'r-', label='Kerrich Data')

# set the x-axis to have log scale
ax.set_xscale('log')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel('Heads - Tails')

# set the title
title = ax.set_title("Kerrich's result was typical!")
title.set_family('serif')

# add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('Simulation of Differences.png')
plt.show()

# <markdowncell>

# So where's the LLN? The LLN does <i>not</i> say that as T increases the number of heads will be close to the number of tails! What the LLN says instead is that, as T increases, the average number of heads will get closer and closer to the long-run average (in this case, 0.5).  The technical term for this is that the sample average, which we estimate from data, converges to the expected value, which is a parameter.
# 
# Let's run another simulation to verify the LLN...

# <codecell>

# set params
N = 100
T = 10000

# generate an array of integers (each of n columns can be interpreted as a run of length T)
data = fair_coin.rvs((T,N))

# create an array in which to store our sample averages
sample_averages = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        sample_averages[i, j] = np.mean(data[:i + 1, j])

# <codecell>

# determine the fraction of heads in the kerrich data
kerrichs_fraction_heads = kerrich_data[:,1] / kerrich_data[:,0]
print kerrichs_fraction_heads

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot each sample path
for i in range(N):
    ax.plot(np.arange(1, T + 1, 1), sample_averages[:, i], 'k-', alpha=0.05)
ax.set_ylim(0, 1)

# plot Kerrich's fraction of heads
ax.plot(kerrich_data[:,0], kerrichs_fraction_heads, 'r-', label=r'$\hat{\mu}_{k}$')

# set the x-axis to have log scale
ax.set_xscale('log')

# demarcate the true mean
ax.axhline(y=0.5, color='black', linestyle='dashed', label=r'$\mu$')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel(r'$\hat{\mu}$')
ylab.set_rotation('horizontal')
ylab.set_fontsize(15)

# set the title
title = ax.set_title(r'$\mathrm{Average\ number\ of\ heads\ converges\ to\ \mu=0.5!}$')
title.set_family('serif')

# add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('Demonstration of LLN.png')

# <markdowncell>

# Exercise 3.2
# ------------
# 
# Re-run the above simulation and regenerate the plots (leaving out Kerrich's data!) using the loaded coin you created in the previous exercise. 

# <codecell>

# Insert your code here!

# <codecell>

# set params
N = 100
T = 10000
loaded_coin = stats.distributions.bernoulli(0.75)

# generate an array of integers (each of n columns can be interpreted as a run of length T)
data = loaded_coin.rvs((T,N))

# create an array in which to store our sample averages
loaded_difference = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        loaded_difference[i, j] = 2 * np.sum(data[:i + 1, j]) - (i + 1)

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot each sample path...
for i in range(N):
    ax.plot(np.arange(1, T + 1, 1), loaded_difference[:, i], 'k-', alpha=0.05)

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')
ylab = ax.set_ylabel('Heads - Tails')

# set the title
title = ax.set_title('Difference really diverges with a loaded coin!')
title.set_family('serif')

plt.savefig('exercise 3.2a.png')
plt.show()

# <codecell>

# set params
N = 100
T = 10000

# generate an array of integers (each of n columns can be interpreted as a run of length T)
data = loaded_coin.rvs((T,N))

# create an array in which to store our sample averages
sample_averages = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        sample_averages[i, j] = np.mean(data[:i + 1, j])

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot each sample path from obs. 10 forward...
for i in range(N):
    ax.plot(np.arange(1, T + 1, 1), sample_averages[:, i], 'k-', alpha=0.05)
ax.set_ylim(0, 1)

# set the x-axis to have log scale
ax.set_xscale('log')

# demarcate the true mean
ax.axhline(y=0.75, color='black', linestyle='dashed', label=r'$\mu$')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel(r'$\hat{\mu}$')
ylab.set_rotation('horizontal')
ylab.set_fontsize(15)

# set the title
title = ax.set_title(r'$\mathrm{Average\ number\ of\ heads\ converges\ to\ \mu=0.75!}$')
title.set_family('serif')

# add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('exercise 3.2b.png')

# <markdowncell>

# Task 4: Where we learn to simulate dice rolling...
# ==================================================
# 
# Suppose that we wish to simulate rolling a standard six-sided die.  Such a simulation would requires to generate random integers between 1 and 6.  To generate random integers in the range [min, max] use stats.distributions.randint(min, max+1) command.

# <codecell>

# create a new rv_frozen object
dice_roll = stats.distributions.randint(1, 7)

# <codecell>

# a single roll of the dice
dice_roll.rvs()

# <codecell>

# a sequence of rolls...
T = 10
sample_rolls = dice_roll.rvs((T,))
print sample_rolls

# <markdowncell>

# We will frequently what to perform statistical analysis of runs from our simulation.  For example, perhaps we wish to estimate the sample average and standard deviation...

# <codecell>

# what is our sample average?
print sample_rolls.mean()
print sample_rolls.std()

# <codecell>

# set the simulation length
T = 1000

# create an array with a single roll of the die
np.random.seed(12345)
data = dice_roll.rvs((1,))
sample_average = data

# for loops simulates the repeated roll of a die
for i in range(1, T):
    next_roll = dice_roll.rvs((1,))
    data = np.append(data, next_roll)
    sample_average = np.append(sample_average, np.mean(data))

# display some data to make sure simulation run looks correct    
print 'Obs:', data[0:3]
print 'Sample Avg.:', sample_average[0:3]

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# create the plot
ax.plot(sample_average, 'b-')
ax.set_xscale('log') # sets x-axis scale to be logarithmic
ax.set_ylim(0, 6)

# demarcate the true mean
mu = dice_roll.mean()
sigma = dice_roll.std()
ax.axhline(y=mu, color='black', linestyle='dashed', label=r'$\mu$')

# label axes
xlab = ax.set_xlabel('Sample Size')
xlab.set_family('serif')

ylab = ax.set_ylabel(r'$\hat{\mu}$')
ylab.set_rotation('horizontal')
ylab.set_fontsize(15)

# set the title
title = ax.set_title('Another demonstration of the LLN')
title.set_family('serif')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Another LLN demonstration.png')
plt.show()

# <markdowncell>

# Exercise 4.1
# ------------
# 
# Increase the length of the simulation run T, to 1000, 10000, and finally 100000.  What do you expect should happen?  Regenerate the plot to see if you inuition is confirmed. You may wish to un-comment the line of code that sets the x-axis scale of the plot to be logarithmic for T=1000 and 100000.

# <codecell>

# insert code here!

# <markdowncell>

# You will probably have notice that each time you regenerate the above plot (assumming you are not setting the seed!) you get a slightly different result (this is true even when T is fixed!).  It would seem that our estimate of the population mean, μˆ, itself has some distribution.  How might we characterize the shape of this distribution? Often characterizing the shape of the distribution of our estimator is more important that the value of the estimator itself.  Why? Because the distribution of our estimator tells us something about how uncertain we are about the specific value of our estimator! Quantifying uncertainty in point estimates is hugely important in econometrics!
# 
# How does our distribution change as we get more data?  Now that we can generate 2D arrays of random integers we can generate a simulation to help us answer this question!  The basic strategy is simply to repeat the above experiment (i.e., rolling a dice T times and re-calculating the sample average after each additional role) some large number, say n, times.
# 
# The most straightforward way to do this simulation (that I could think of at least!) uses a double for loop...

# <codecell>

# set params
N = 1000
T = 1000

# generate an array of integers (each of n columns can be interpreted as a run of length T)
data = dice_roll.rvs((T,N))

# create an array in which to store our sample averages
sample_averages = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        sample_averages[i, j] = np.mean(data[:i + 1, j])

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# calculate the standard errors
std_error = np.zeros(T)
for j in range(T):
    sigma = np.std([1,2,3,4,5,6])
    std_error[j] = sigma / np.sqrt(j + 1)

# plot the sample paths
for i in range(N):
    ax.plot(sample_averages[:, i], 'b-', alpha=0.01)
ax.set_ylim(0, 7)

# set the x-axis to have log scale
ax.set_xscale('log')

# add some lines...can you figure out what I am plotting with these commands?
ax.axhline(y=mu, color='black', linestyle='dashed', label=r'$\mu$')
ax.plot(mu + 1.96 * std_error, 'orange', ls='--', label='95% C.I.')
ax.plot(mu - 1.96 *std_error, 'orange', ls='--')

# label axes
xlab = ax.set_xlabel('Sample Size')
xlab.set_family('serif')

ylab = ax.set_ylabel(r'$\hat{\mu}$')
ylab.set_rotation('horizontal')
ylab.set_fontsize(15)

# set the title
title = ax.set_title('Another LLN demonstration')
title.set_family('serif')

# add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('A further LLN demonstration.png')
plt.show()

# <markdowncell>

# Imagine gathering data by taking a vertical slice through the above plot.  This would be a marginal distribution and would look something like...

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# add the histogram of the first row of sample averages. Why?
n, bins, patches = ax.hist(sample_averages[0,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

# don't forget to label axes!
ax.set_xlabel(r'$\hat{\mu}$')
ylabel = ax.set_ylabel('Density')
ylabel.set_family('serif')

# add a title!
title = ax.set_title('Marginal distribution for T=1 for N=1000 trials')
title.set_family('serif')

plt.savefig('Initial roll of dice has uniform distribution.png')
plt.show()

# <markdowncell>

# Now let's plot a histogram of the marginal distribution for $T=2$.

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# add the histogram of the first row of sample averages. Why?
n, bins, patches = ax.hist(sample_averages[1,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

# CLT says that distribution should be Gaussian for large enough T
x = np.linspace(0, 7, 1000)
normal = stats.norm(loc=mu, scale=sigma / np.sqrt(2))
ax.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ax.set_xlabel(r'$\hat{\mu}$')
ylabel = ax.set_ylabel('Density')
ylabel.set_family('serif')

# add a title!
title = ax.set_title('Marginal distribution for T=2 for N=1000 trials')
title.set_family('serif')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Getting closer to normal.png')
plt.show()

# <codecell>

# create new Figure object
fig = plt.figure()

# create the first subplot
ax1 = fig.add_subplot(321)

# add the histogram of the first row of sample averages.
n, bins, patches = ax1.hist(sample_averages[0,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

# don't forget to label axes!
ylabel = ax1.set_ylabel('Density')
ylabel.set_family('serif')

# add a title!
title = ax1.set_title('T=1')
title.set_family('serif')

# create the second subplot
ax2 = fig.add_subplot(322)

n, bins, patches = ax2.hist(sample_averages[49,:], bins=25, align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(49))
ax2.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ylabel = ax2.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax2.set_title('T=50')

# create the third subplot
ax3 = fig.add_subplot(323)

n, bins, patches = ax3.hist(sample_averages[1,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(2))
ax3.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ylabel = ax3.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax3.set_title('T=2')

# create the fourth subplot
ax4 = fig.add_subplot(324, sharex=ax1)

n, bins, patches = ax4.hist(sample_averages[499,:], bins=25, align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(500))
ax4.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ylabel = ax4.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax4.set_title('T=500')

# create the fifth subplot
ax5 = fig.add_subplot(325, sharex=ax1)

n, bins, patches = ax5.hist(sample_averages[2,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(3))
ax5.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ax5.set_xlabel(r'$\hat{\mu}$')
ylabel = ax5.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax5.set_title('T=3')

# create the sixth subplot
ax6 = fig.add_subplot(326, sharex=ax1)

n, bins, patches = ax6.hist(sample_averages[999,:], bins=25, align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(1000))
ax6.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ax6.set_xlabel(r'$\hat{\mu}$')
ylabel = ax6.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax6.set_title('T=1000')

# set the layout
plt.tight_layout()

# add a title for the entire set of plots!
plt.figtext(0.5, 0.95, 'CLT in Action!', ha='center', weight='bold', size='large')

plt.savefig('The CLT in action.png')
plt.show()

# <markdowncell>

# Task 5: Relationships between probability distributions 
# ========================================================
# 
# In this task you will learn about some interesting relationship between probability distributions via simulation, and look at an example in which the LLN fails in a spectacular fashion.

# <codecell>

# insert your code here

# <markdowncell>

# The sum of squares of $k$ standard Normal random variables is a Chi-square random variable with $k$ degrees of freedom. Technically, $||N_{i=1,\dots,k}(0,1)||^2 = \chi^2(k)$

# <codecell>

# create a standard normal distribution object
std_normal_rv = stats.distributions.norm()

# create a chi-square distribution object with k degrees of freedom
k = 5
chi_square_rv = stats.distributions.chi2(k)

# <codecell>

# set the length of our random samples...
T = 1000

# generate an array of random draws from a standard normal
std_normal_samples = np.random.randn(T, k)

# <codecell>

# generate chi_square samples by taking the sum of squared elements each row
chi_square_k_samples = np.sum(std_normal_samples**2, axis=1)

# <codecell>

# create a grid for pretty plotting
grid = np.linspace(0, 15, 1000)

# create a new figure object
fig = plt.figure()

# create the first subplot
ax = fig.add_subplot(111)
        
# histogram of manufactured chi-squares
np.random.seed(12345)
n, bins, patches = ax.hist(chi_square_k_samples, bins=50, normed=True, color='blue', \
    label=r'$\sum^{k}_{i=1}N^2_{i}(0,1)$')
ax.plot(grid, chi_square_rv.pdf(grid), 'r--', label=r'$\chi^2(k)$')
ax.set_xlim(0, 15)

# don't forget to label Axes!
ax.set_xlabel('X')
ax.set_ylabel('Density')

# set the title
ax.set_title(r'$\mathrm{The\ sum\ of\ squares\ of\ k\ N(0, 1)\ random\ variables\ is\ a\ \chi^2(k)}$!')
 
# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Relationship between Normal and Chi-square.png')
plt.show()

# <markdowncell>

# [F-distribution](http://en.wikipedia.org/wiki/F-distribution) with degrees of freedom $(df_1, df_2)$ is the ratio of two $\chi^2(k)$ random variables with degrees of freedom $df_1$ and $df_2$ repsectively.

# <codecell>

# specify the degrees of freedom
df1 = 100
df2 = 100

# create an F distribution object
F_rv = stats.distributions.f(df1, df2)

# create two Chi-square objects
chi_square_n = stats.distributions.chi2(df1)
chi_square_d = stats.distributions.chi2(df2)

# <codecell>

# specify a length for random samples
T = 1000

# create a grid for pretty plotting
grid = np.linspace(0, 3, 1000)

# create a new figure object
fig = plt.figure()

# create the first subplot
ax1 = fig.add_subplot(111)
        
# histogram of random draw from ratio of chi-squares
np.random.seed(12345)
manufactured_F_rv = (chi_square_n.rvs((T,)) / df1) / (chi_square_d.rvs((T,)) / df2)
n2, bins2, patches2 = ax1.hist(F_rv.rvs((T,)), bins=25, normed=True, color='blue', \
    label=r'$\frac{\left(\frac{\chi^2(df_1)}{df_1}\right)}{\left(\frac{\chi^2(df_2)}{df_2}\right)}$')
ax1.plot(grid, F_rv.pdf(grid), 'r--', label=r'$F(df_1, df_2)$')

# don't forget to label Axes!
ax1.set_xlabel('X')
ax1.set_ylabel('Density')

# set the title
ax1.set_title('The ratio of two Chi-squares is an F!')
 
# add a legend
ax1.legend(loc='best', frameon=False)

plt.savefig('Relationship between Chis-square and F.png')
plt.show()

# <markdowncell>

# The ratio of two standard Normal random variables, has a [Cauchy](http://en.wikipedia.org/wiki/Cauchy_distribution) distribution.  The Cauchy distribution has the curious property that all of its moments (i.e., mean, variance, skewness, etc) are infinite! Which means, among other things, that the LLN does not hold... 

# <codecell>

# create a standard Normal distribution object
normal_rv = stats.distributions.norm()

# create a Cauchy distribution object
cauchy_rv = stats.distributions.cauchy(0, 1)

# <codecell>

# Set the length of our data array
T = 10000

# create an array of random draws from our standard normal
normal_data = normal_rv.rvs((T, 2))

# <codecell>

# create a grid for pretty plotting
grid = np.linspace(-10, 10, 1000)

# create a new figure object
fig = plt.figure()

# create the first subplot
ax1 = fig.add_subplot(111)
        
# histogram of random draw from ratio of chi-squares
np.random.seed(12345)
manufactured_Cauchy_rv = normal_data[:,0] / normal_data[:,1]

# need to have many more bins than samples to get good histogram
n2, bins2, patches2 = ax1.hist(manufactured_Cauchy_rv, bins=5000, normed=True, color='blue', \
    label=r'$\frac{N(0,1)}{N(0,1)}$')
ax1.plot(grid, cauchy_rv.pdf(grid), 'r--', label=r'$Cauchy(0, 1)$')
ax1.set_xlim(-10,10)

# don't forget to label Axes!
ax1.set_xlabel('X')
ax1.set_ylabel('Density')

# set the title
ax1.set_title('The ratio of two standard Normals is a Cauchy!')
 
# add a legend
ax1.legend(loc='best', frameon=False)

plt.savefig('Relationship between Normal distribtion and Cauchy.png')
plt.show()

# <codecell>

# set length of sample
np.random.seed(10)
T = 10000

# generate some data
cauchy_data = cauchy_rv.rvs(T)

# create a vector in which to store sample averages
sample_average = np.zeros(T)

# for loops simulates the repeated roll of a die
for i in range(0, T):
    sample_average[i] = np.mean(cauchy_data[:i + 1])

# <codecell>

# create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# create the plot
ax.plot(sample_average, 'g-')
ax.set_xscale('log') # sets x-axis scale to be logarithmic

# demarcate the true location
mu = 0
ax.axhline(y=mu, color='black', linestyle='dashed', label=r'$\mu$')

# label axes
xlab = ax.set_xlabel('Sample Size')
xlab.set_family('serif')

ylab = ax.set_ylabel(r'$\hat{\mu}$')
ylab.set_rotation('horizontal')
ylab.set_fontsize(15)

# set the title
title = ax.set_title('Epic Failure of the LLN')
title.set_family('serif')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Failure of LLN.png')
plt.show()

# <codecell>


