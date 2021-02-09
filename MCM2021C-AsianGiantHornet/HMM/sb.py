#!/usr/bin/env python
# coding: utf-8

# In[ ]:


y = np.asarray(observ)
s =  2803
start_probability = np.zeros([n])
start_probability[s] = 1
#start_probability = np.matrix(start_probability).reshape([-1,1])
#start_probability = start_probability.reshape([-1,1])

nu = start_probability
Q,g,l = HMMbaumwelch(y, nu, tol=1e-4, maxIt = 100)


# In[ ]:


my_data = from numpy import genfromtxt
my_data = genfromtxt('',delimiter = ',')


# In[ ]:


transition_probability = my_data
emission_probability = np.zeros([n,2])
p = 0.9
for i in range(n):
    emission_probability[i,0] = p
    emission_probability[i,1] = 1-p


# In[ ]:


def HMMsample(nu, Q, g, n):
  '''
sample a trajectory from a hidden markov chain
in:   nu = initial distribution as vector of size k
      Q = transition matrix of size k
      n = positive integer
out:  (x,y) = sample trajectory of size n of a HMM defined by (nu, Q, g):
      x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
      y = observations such that the conditionnal distribution of y[k]
      given x[k] is g(x[k], :)
  '''
  x = np.zeros(n)
  y = np.zeros(n)
  x[0] = randm(nu)
  y[0] = randm(g[x[0], :])
  for j in range(1,n):
    x[j] = randm(Q[x[j-1], :])
    y[j] = randm(g[x[j], :])
  return (x,y)

def HMMfilter(y, nu, Q, g):
  '''
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(g.shape[1])
      nu = initial distribution as vector of size k
      Q = transition matrix of size k x k
      g = emission matrix with k rows
out:  phi = filter: P(x[t]=x | y[0:t]=y[0:t]) for 0<=x<k and 0<=t<n
      c(t) = conditional likelihood: P(Y[t] = y[t]| Y[0:t-1]=y[0:t-1])
  '''  
  n = y.size
  k = Q.shape[0]
  phi = np.zeros([k, n])
  c = np.zeros(n)
  
  Z = nu * g[:, y[0]]
  c[0] = sum(Z)
  phi[:,0] = Z/c[0]
  
  for t in range(1,n):
    Z = np.dot(phi[:, t-1], Q) * g[:, y[t]];
    c[t] = sum(Z)
    phi[:, t] = Z / c[t]
  
  return (phi,c)


def HMMsmoother(y, Q, g, c):
  '''
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(Q.shape[0])
      Q = transition matrix of size k x k
      g = emission matrix with k rows
      c = conditional likelihoods, computed by HMMfilter
out:  beta = smoothing factors: P(y[t+1:n]=y[t+1:n] | X[t]=x) / P(Y[t+1:n]=y[t+1:n] | Y[0:t]=y[1:t]) for 0<=x<k and 1<=t<n
      permits to compute the posterior distribution of the hidden states 
      P(X[t]=x | Y[0:n]=y[0:n])  as post = phi * beta
  '''
  n = y.size
  beta = np.ones([Q.shape[0], n])
  for t in range(n-2, -1, -1):
    beta[:,t] = np.dot(Q, g[:, y[t+1]] * beta[:, t+1]) / c[t+1]
  return beta

def HMMbaumwelch(y, nu, tol=1e-4, maxIt = 100):
  '''
compute maximum likehood estimate using Expectation-Maximization
iterations
in:   y = vector of observations 
      nu = initial distribution of the hidden chain
      tol = tolerance for the stopping criterion
      maxIt = maximal number of iterations
out:  Q = estimate of the transition matrix of the hidden markov process
      g = estimated probabilities of transition: g(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
      l = log-likelihood of y for parameters Q and g  
  '''  
  #global myfilter, mysmoother # should be either HMMfilter/HMMsmoother, or HMM_C.HMMfilter/HMM_C.HMMsmoother
  k = nu.size; r = 1+max(y); n = y.size
  Y = np.zeros([n, r]); 
  Y[range(n), np.int_(y)] = 1;
  Q = transition_probability
  g = emission_probability
  for j in range(k):
    Q[j,:] = Q[j,:] / sum(Q[j,:])
    g[j,:] = g[j,:] / sum(g[j,:])
  it = 0; oldQ = Q; oldg = g + tol +1
  
  while (sum(sum(abs(oldQ[:]-Q[:]))) + sum(sum(abs(oldg-g))) > tol) & (it<maxIt):
    it+=1
    (phi, c) = HMMfilter(y, nu, Q, g)
    beta = HMMsmoother(y, Q, g, c)
    post = phi * beta
    N = Q * (np.dot(phi[:, 0:-1], np.transpose(beta[:, 1:]*g[:, np.int_(y[1:])]/np.tile(c[1:], [k, 1]))))
    M = np.dot(post, Y)

    oldQ = Q.copy(); oldg = g.copy()
    for j in range(k):
      Q[j,:] = N[j,:] / sum(N[j,:])
      g[j,:] = M[j,:] / sum(M[j,:])
  l = sum(np.log(c))
  return (Q, g, l)

