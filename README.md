# Computational parametric mapping of functional neuroimaging data

## Abstract

The analysis of functional neuroimaging data and its integration with computational models of cognition typically depends on fitting models to behaviour and regressing the latent variables of those models onto the neural data. Here we develop a technique that directly fits computational cognitive models to neuroimaging data. This technique performs a voxel-wise comparison of cognitive models and visualises the latent parameters of these models over neuroanatomical space. This allows cognitive spaces to be topographically mapped analogous to the mapping of sensory spaces such as retinotopic space. Furthermore, fitting models directly to neural data rather relaxes the assumption that the cognitive variables of interest must be expressed in behavior. Our approach builds on a probabilistic framework designed for mapping population receptive fields (pRFs) in sensory cortices. The advance we offer is modest but useful. First, we generalise the pRF modelling framework from the sensory to the computational cognitive domain. Second, we substantially lower the computation time per voxel, making it feasible for modelling more extensive neural systems and over longer time series. Finally, we offer simulations for model comparison and parameter recovery, demonstrating the robustness and specificity of this technique. We assert that this approach, which we generically call computational parametric mapping, can play a role in the discovery of topographic principles underlying the neural encoding of cognitive computational variables. 

## Introduction

**Model-based neuroimaging.** Deploying computational models of the mind, brain, and behaviour to analyse functional neuroimaging data is now a powerful and widely used strategy in the cognitive and computational neurosciences [ref find a review]. Since the advent of this strategy, the methods used have typically depended on fitting computational models of cognition to behavioural data first and then regressing the latent variables of these fitted models onto functional neuroimaging data. One of the first examples of this came from reinforcement learning, in which a temporal difference learning algorithm was first fit to choice behaviour, resulting in the estimation of a learning rate parameter determining how quickly rewards are learnt as a function of experience [nature, ben seymour et al.]. This fitted model yielded latent variables such as reward prediction errors that was entered as predictor variables into a first-level multiple linear regression model. Computing the main effect of this predictor variable thus revealed brain regions associated with the behaviorally-fitted reward prediction error variable.
   
**Relative merits of behavior-centric models.** This approach has the advantage of deploying regression techniques predominant in functional neuroimaging and for which efficient algorithms exist for analysing large volumes of fMRI data. Furthermore, fitting the computational model to behaviour first ensures that the neuronal populations modelled are potentially behaviorally relevant. The downside of this behavior-centric fitting approach is that it necessarily restricts the mapping of computational variables to those associated with observable behaviour. The brain is a massively parallel organ, the elements of which may not necessarily map to ongoing behaviour in any simple way. Thus, only mapping those latent cognitive variables that fit observable behaviour substantively narrows the scope of neuronal populations that can be modelled and topographically mapped.

**Generative modelling as an alternative.** A different approach has been pioneered in sensory neuroimaging. Kumar and Penny [2014] developed a Bayesian method that allows the estimation of neural response functions directly from fMRI data. The two building blocks of the method are a neural response function  which maps stimuli characteristics to neuronal population responses, and a hemodynamic response function which maps neural responses to fMRI data. The latter is an observation model as it models from latent causes to observable data. The parameters of the two response functions are estimated together using a Bayesian optimisation algorithm known as variational Laplace, widely used in neuroimaging [friston et al. 2007]. This method was initially used to model the receptive field structure of neuronal populations in auditory cortices, and visualise to their topographic mapping [Kumar and Penny, 2014]. 


**Benefits of the generative approach.** There are many benefits to this approach. First, the neural and hemodynamic response functions can be fitted to neural data without needing any behavioural information. This widens the scope for fitting neural populations which may not be expressed in the behavior recorded, or whose association to the behavior may not be adequatedly modelled. Second, by virtue of being a Bayesian approach, the models can easily incorporate prior knowledge and generalise to hierarchical models. Third, by producing a model evidence score that can be used to formally compare models evaluating competing neural response functions, balancing the accuracy against complexity. Finally, by including a hemodynamic model for every voxel, the method can accomodate varations across brain regions and subjects, including non-linearities associated with hemodynamic saturation effects. 


**Population receptive field toolbox.** The method of [Kumar and Penny 2014] was subsequently generalised to model visual (population) receptive fields [Zeidman et al. 2018] and released as the BayespRF toolbox available for the Statistical Parametric Mapping package [ref]. Henceforth, we refer to this method as the pRF method. 

**Our aims.** In this technical note, we are interested in developing the pRF method and adapting it for voxel-wise estimation of cognitive computational models. In other words, we wish to generalise receptive field modelling into cognitive field modelling. In the same way a receptive field model maps the parameters of a receptive field model onto the brain, we wish to generalise this approach to any computational cognitive model. We begin by describing the generative models of the existing pRF framework and we explain how this generalises to computational parametric mapping. We illustrate this via a comparison between sensory and cognitive models, and then step through the application of the approach in detail for one specific reward model. We discuss an approach to speed up the computation time, and demonstrate robustness via parameter and model recovery simulations.


## Generative models

**Response functions.** The pRF framework allows voxel-wise fitting of a non-linear model response to fMRI data. The model is a composition of  two seperate response functions, namely the neuronal response function (NRF) and the hemodynamic response function (HRF) (Fig. 1). The NRF is based on a non-linear parametric form which is free to vary depending on the application. It map from experimental inputs $x$ such as stimulus or behavioral events through to a neuronal response $z$, parameterised by unknown neuronal parameters $\theta_n$. The HRF maps from this neuronal response $z$ to the predicted fMRI activity $g$ via an extended balloon model with unknown parameters $\theta_h$. The whole model (NRF and HRF) is fitted to fMRI data $y$ in order to estimate the unknown parameters. Of interest to the cognitive computational neuroscientist will be the estimated of NRF parameters $\theta_n$, as well as model comparisons between different NRFs. 

|![](https://i.imgur.com/rKjDPtD.jpg)|
| :-: |
|**Fig. 1 Schematic of generative model.** The generative model maps from experimental input, to predicted neuronal input via the NRF, and from neuronal input to predicted fMRI signal via the HRF. The parameters of the model are fit according to how well the generative model predicts the observed fMRI data.| 

**Receptive field model.** We begin by formally outlining a receptive field model equivalent to that specified in the pRF framework. This provides a foundation from which we can draw out the correspondence to more abstract cognitive models in the next section. In Zeidman and colleagues [ref] the model was set up to model retinotopic receptive fields. *[Note that we set up the notation for the visual inputs in a way that allows for generalisation to a wider class of compuational models]*. A receptive field mapping approach maps the topographic distribution of parameters of receptive field models onto the brain. It may also map the topographic distribution of neuronal populations that may be better explained by one receptive field model compared to another. 

**Receptive field model parameters** If we start with a simple gaussian receptive field model that takes three parameters, two for the pixel coordinates $\mu^x$ and $\mu^y$ and one for the variance $\mu^\sigma$. These three parameters define a parameter space for this receptive field model $\mu=\{ \mu^x, ∈ \mathbb{R}, \mu^y ∈ \mathbb{R},\mu^\sigma ∈ \mathbb{R^+}\}$. 

**Experimental inputs for receptive field model.** The experimental inputs to the model are simply the luminance values $x$ and their positions $k$ such that that $x(i,t)$ is the luminance value of a pixel located at position $k(i)$ at timestep $t$.

**Population response space for a receptive field model.** We now simulate the response of neurons to this experimental input, where the population of neurons span the parameter space of the model. This yields a matrix representing a population response over time. For a grid of discrete parameter values $\mu$ this yields a matrix $s$ such that $s(i,t)$ is the response of a neuron with parameters $\{\mu^x(i),\mu^y(i),\mu^\sigma(i)\}$ at the $t^{th}$ timestep.

**Population field.** Having set up the population of hypothetical neurons with different receptive field properties, we can now define a population field, which is a (normalised) gaussian weighted average of these neuronal responses, whose parameters are optimised to best account for the voxelwise response:

$z(t)=\beta\sum_i s(i,t) \cdot N(k(i):\theta_n),$  

$N(k(i):\theta_n)=N(k(i):\mu,\Sigma)=\frac{1}{\sqrt{(2\pi)^3|\Sigma|}} \exp\left(-\frac{1}{2}(k(i)-\mu)^T{\Sigma}^{-1}(k_i-\mu)\right).$

The population field is similar to a receptive field however here it is more abstractly defined as the weighting of a population response space to best explain a voxel response. In other words it models what mixture of the hypothetical neuronal population is necessary to beest explain the voxels response. 

**Cognitive field model.** We now follow an equivalent setup for a cognitive model, here illustrated for a simple reward experiment. 


We begin again with the input $x(t)$, which for our exemplar reward experiment is determined simply by the reward values and their associated sensory cues.We explicitly organise the time-series of inputs into a cognitive space, effectively this is the parameter space of a cognitive model. The cognitive space for our reward learning model has two dimensions, one for a risk aversion parameter, and another for a learning rate parameter. Each location in this cognitvie space is represented by its coordinates $k_i=\{x ∈ \mathbb{R},y ∈ \mathbb{R}\}$. Expressed in this cogntive space, the cognitive input can be represented as $s$ such defined such that $s(i,t)$ is the reward prediction error variable for an agent with the parameter values $k_i$ at timestep $t$. The neuronal response $z$ at time step $t$ is then modelled simply as a gaussian weighted average of these cognitive inputs, which is linearly scaled by a parameter $\beta$. The same equations apply as in the receptive field model:   

$z(t)=\beta\sum_i s(i,t) \cdot N(k_i:\theta_n),$ 

$N(k_i:\theta_n)=N(k_i:\mu,\Sigma)=\frac{1}{\sqrt{(2\pi)^2|\Sigma|}} \exp\left(-\frac{1}{2}(k_i-\mu)^T{\Sigma}^{-1}(k_i-\mu)\right).$

What was the gaussian receptive field is now in this context treated as a cognitive field, $N(k_i:\theta_n)$. The neuronal response function is thus again parameterised by $\theta_n$ comprised of a vector $\mu$ which specifies the location of the cogntive field and its covariance matrix $\Sigma$ which specifies the width and rotation of the cognitive field. 

**Comparing the receptive fields and cognitive fields.** To help compare the receptive field models and cognitive models we juxtapose the two processes in Fig. 2. We make explicit the mapping from an arbitrary set of inputs $s$ to a coherent space, in the form of sensory or cognitive spaces. This way of visualising the process makes explicit the equivalence between the two domains. Note that the HRF for both domains is identical. 
| ![test](https://i.imgur.com/lX6C1qg.jpg) | 
| :-: |
|**Fig. 2. Comparing sensory and cognitive field models.** The upper row represents the generic inputs, response functions, parameters and outputs. Middle two rows illustrate the receptive field and cognitive field models, respectively. |

### Model fitting 

**Introducing variational Bayes.** The posterior distribution of free parameters $\theta = (\theta_n,\theta_h)$ is found by fitting the synthetic BOLD signals to the measured BOLD signals for each voxel using variational Bayes as an approach to approximate Bayesian inference. Variational Bayes provides an efficient means for parameter estimation and model selection. Here we summarise this process.

**Summary of variational Bayes.** Let $m$ be a generative model of data $y$, defined by the likelihood function $p(y|\theta,m)$ and the prior density $p(\theta|m)$ on model parameters $\theta$. In brief, VB is an iterative scheme that indirectly optimizes an approximation to both the model evidence $p(y∣m)$ and the posterior density $p(\theta∣y,m)$. The key is to decompose the log model evidence into:

$log \, p(y|m)=F(q)+D_{KL}(q(\theta)||p(\theta|y,m))$

$F(q)=log \, p(y|m)-D_{KL}(q(\theta)||p(\theta|y,m)),$

where $D_{KL}$ is the Kullback–Leibler divergence, $q$ is the approximate distribution and $F$ is the variational free energy. By maximizing $F(q)$ we are minimizing the divergence between the posterior distribution $p(\theta|y,m)$ and the approximate distribution $q(\theta)$. This maximisation is done by way of gradient ascent. Given two concurrent models $m_1$ and $m_2$, the factor of their minimized free energy $F_1^*$ and $F_2^*$ approximates the Bayes Factor between models :

$B_{12}=\frac{p(y|m_1)}{p(y|m_2)}\approx exp(F_1^*-F_2^*).$


This allows us to perform Bayesian model comparison and selection by choosing the model with the strongest model evidence.

**Variational Laplace.** The variational Laplace (VL) method uses simplifying assumptions on the approximate distribution $q(\theta)$, which render the variational Bayes algorithm computationally and statistically efficient [cite]. One of these assumptions, the Laplace approximation, requires the model's priors to be Gaussian. As such, we will be expressing the cognitive field locations $\mu$ and widths $\Sigma$ as latent Gaussian parameters. Here, we represent our uncertainty in the parameters by constraining them to uniforms distributions on a plausible range of real numbers. Thus, for each parameter, location $\mu$ and width $\sigma$ are expressed in terms of latent parameters $l_{\mu}$ and $l_{\sigma}$ with normal priors:

$\mu = (\mu_{max}-\mu_{min})*ncdf(l_{\mu})+\mu_{min}$
$\sigma = (\sigma_{max}-\sigma_{min})*ncdf(l_{\sigma})+\sigma_{min},$

where $ncdf$ is the cumulative density function for the univariate normal distribution, and $\mu_{min},\mu_{max}, \sigma_{min},\sigma_{max}$ are the boundaries where the parameters can take their values. $l_{\mu}$ and $l_{\sigma}$ are the latent parameters with normal priors. The described transformations turn the normal distribution of a latent parameter into a uniform distribution of the parameter across its support:

$l_{\mu}\sim\mathcal{N}(0, 1) \implies \mu\sim\mathcal{U}(\mu_{min}, \mu_{max})$

$l_{\sigma}\sim\mathcal{N}(0, 1) \implies \sigma\sim\mathcal{U}(\sigma_{min}, \sigma_{max}).$

| ![test](https://i.imgur.com/HQzYphd.png) | 
| :-: |
|**Figure 3** : Transformation of latent parameter $l_\mu$ distribution to the parameter $\mu$ distribution. The prior latent distribution is a unit Gaussian (mean 0, standard deviation 1). The posterior latent distribution is Gaussian with a mean 0.7 and a standard deviation of 0.1. The posterior distribution is not uniform, but approximatively Gaussian with a mean of ~2, when the latent parameter standard deviation is low enough.|

The scaling parameter $\beta$ is constrained to be positive using a normal latent parameter $l_{\beta}$ and through the transformation $\beta=exp(l_{\beta})$ .

| ![test](https://i.imgur.com/1arOXk3.png) | 
| :-: |
|**Fig. 3** : Transformation of latent parameter $l_\beta$ distribution to the parameter $\beta$ distribution. The prior latent distribution is a unit Gaussian. The posterior latent distribution is Gaussian with mean 0.7 and standard deviation 0.1. The posterior distribution is not uniform, but approximatively Gaussian with a mean of ~2, when the latent parameter standard deviation is low enough.|

### Use of a Precomputed Meshgrid

During the fitting, using the VL algorithm, the computing of $s_{\theta_N}$ is performed a large number of times for different sets of parameters. This process can be computationally extensive, especially in the context of models of reward. To avoid excessive computation, the values of $s_{\theta_N}$ over the experimental session is precomputed on a discretised subset of the model parameter space. 

A finite support $[p_{min},p_{max}]$ and a bin size $r$ are used to define a discrete subset of values for each parameter $p$. The choice of the bin size depends on available computing power and the desired precision. The choice of the boundaries depends on a reasonable hypothesis given existing data. A meshgrid $\Theta$ over the neuronal parameter space is obtained, over which $s_{\theta_N}$ is computed at every point. 

| ![test](https://i.imgur.com/Rt3i1XE.png) | 
| :-: |
|**Fig. 4 Meshgrid parametrization for a bidimensionnal parameter space** |


The precomputed meshgrid $\Theta$ has to be chosen such that the centre of the population receptive field (pRF) falls inside its bounds. As a rule of thumb, we defined it such that at least 95% of prior density fall inside the meshgrid area, and such that at least a grid point is contained around the centre of the receptive field.  This gives for each parameter $p$:

- $[p_{min},p_{max}]  \supset[\mu_{min}-2\sigma_{max},\mu_{max}+2\sigma_{max}]$ 
- $r \leq 2\sigma_{min}$



### Observation model

It is not possible to observe neuronal activity directly. The neuronal signal is inferred from a secondary phenomenon that also needs to be modelled. In the case of fMRI, it is the hemodynamic response.

The default hemodynamic model in the BayespRF toolbox is an extended balloon model. Most parameters in this model are fixed, based on empirical measurements, leaving three parameters to be estimated:  the transit time $\tau$, the rate of signal decay $\kappa$, and the ratio of intra-to extra-vascular signal $\epsilon$. They are also expressed as the transformation of latent parameters similar to $\beta$. We will be using this model for what follows.


## Cognitive Field Models of Reward

To illustrate the framework above, we will apply it using a reinforcement-learning algorithm as the cognitive model. 

### Experiment setting:

We consider the following experiment involving a passive reward learning task [2]: subjects were endowed with a wealth of 1000DKK and shown a series of fractal images, each having a deterministic effect on their wealth. The subject is tasked with learning the effect of each fractal via passive observation of the changes in wealth. A session consists of 150 trials with the events in a single trial consisting of:
- <ins>beginning of the trial</ins>, a cue for the start of the trial
- <ins>fractal onset</ins>, a fractal picture
- <ins>wealth variation</ins>, a visible update of the subject's wealth

| ![test](https://i.imgur.com/wNPSZIw.png) | 
| :-: |
|**Figure 5** : Evolution of a subject wealth over time during an experiment session. |!

### Temporal difference learning:

The learning process of a subject can be modelled by a temporal difference (TD) algorithm. 

| ![test](https://i.imgur.com/1VKhebU.png) | 
| :-: |
|**Figure 4** : Computing the reward prediction error with the TD learning algorithm |

Learning in this context is driven by a progressive evolution of expected utility $U$ over trials, caused by the sequence of fractals and their effect on endowed wealth. We are especially interested in this algorithm’s reward prediction error (RPE) signal, as it has been empirically associated with phasic dopaminergic activity of the midbrain. Faced with this paradigm, our model learns an associative strength vector $W$, representing the 'value function' in this context.

The utility here represents the subjective value perceived by a subject to the wealth variation he perceives during the experimental session. This can be modelled as a transition from an amount of money $C_t$ to an amount $C_{t+1}$ by an isoelastic utility function $iso$.

$U_{t+1}=iso(C_{t+1},\eta)-iso(C_{t},\eta)$

with $iso(c,\eta)=\frac{c^{1-\eta}}{1-\eta}$ for $\eta \neq1$ 
and $iso(c,1)=ln(c)$

$\eta$ is termed the aversion degree.

Here we define a standard TD algorithm and use it as the neuronal model, with the reward prediction error as the output signal. This model takes into account the value of rewards being delivered at each time step:

$W_{t+1}= W_{t}+\alpha E(t)Z_t$
$E(t)=U_{t+1}+\gamma W_tX_{t+1}-W_tX_t$
$Z_{t+1}=\gamma \lambda Z_t+X_t$

Where $W_t$ is the associative strength vector, $\alpha$ is the learning rate, $E$ is the reward prediction error, $Z_t$ is the eligibility-trace vector which acts as a sort of backwards memory for the TD algorithm [cite], $U_t$ is the utility, $X_t$ is the complete serial compound (CSC) feature vector which represents the stimulus [cite], $\gamma$ and $\lambda$ are decay rates.

| ![test](https://i.imgur.com/u854FPB.png) | 
| :-: |
|**Figure 4** : Computing the reward prediction error $E(t)$ with the TD learning algorithm |



### Cognitive models:

We are interested in two parameters: the learning rate $\alpha$ and the risk aversion degree $\eta$. We fix all other parameters to standard values:
- Discount factor $\gamma=0.97$
- Eligibility trace decay rate $\lambda=1$

We are also interested in comparing different variants of the proposed model. 

We are interested in comparing different variations of the proposed TD model. For instance, it has been demonstrated that the RPE signals in mice VTA neurons follow a distributional code [cite]. Under this class of learning model, known as distributional reinforcement learning, we have two learning parameters $\alpha^+$ and $\alpha^-$ for positive and negative RPEs respectively. Also, we are interested in knowing if including the aversion degree (the utility function parameter) as a parameter is worth the extra model complexity. 

With these two considerations, we have a total of 4 possible cognitive models outlined in this table:
| ![test](https://i.imgur.com/N2SwIKk.png) | 
| :-: |
|Figure 10 : Cognitive models of reward |


To compare evidence for the data under different models, we compute the free energy. For each voxel, the VL algorithm returns the free energy of each model given the data observed. We use VBA (Variational Bayesian Analysis) toolbox to perform a random-effects RFX analysis using a $N_{voxels} \times N_{models}$ matrix of free energies. This is usually used to compare subjects but can also be applied across a population of voxels within a single subject’s brain. The analysis returns the expected probabilities of each model, which measures how likely it is that the model is more frequent than all other models in the comparison set.


### Model n°3 Specification

We will begin by specifying model 3. The learning rate $\alpha$ is between 0 and 1, so we transform the parameter in one varying on the full real line with a logit transformation , $\tau = logit(\alpha) = log(\frac{\alpha}{1-\alpha})$. This gives us a proper bidimensional learning space. The meshgrid $\Theta$ is defined by setting: $\mu_{min}=-2,\mu_{max}=2,\sigma_{min}=0.1,\sigma_{max}=1,R=2,r=0.2$ for both $\tau$ and $\eta$.

The reward prediction error signal is composed of spikes of various intensity, occurring at stimuli timings. Elsewhere it takes on a null value. The neural response function is then:

$z(t)=\beta\sum_{j=1}^N \delta(t-t_j)\sum_{k \in \Theta} E_k(t)N(k:\boldsymbol\mu,\boldsymbol\Sigma)$

Where $t_j$ is the timing of the j-th stimulus, $E_k(t)$ is the reward prediction error signal for parameters $k=(\tau,\eta)$ at time $t$, $\boldsymbol\mu=(\mu_{\tau},\mu_{\eta})$ is the location of the receptive field and $\boldsymbol\Sigma=diag(\sigma_{\tau}^2,\sigma_{\eta}^2)$   is its covariance matrix. 

| ![test](https://i.imgur.com/qxFi4Y5.png) | 
| :-: |
|**Figure 6** : Prior & Posterior probability density |

In total, we have eight parameters: five neuronal plus three haemodynamic. They are all defined through transformations of latent parameters with normal priors. 



### Parameter Estimation and Recovery

We generate data with the model for a set of parameters and add various noise levels to the BOLD signal before trying to recover the parameters. We generate nine voxels time series by taking the location parameters to be all the pairs $(\mu_\tau,\mu_{\eta})$ of sets $A_{\tau}=\{-1,0,1\}$ and $A_{\eta}=\{-1,0,1\}$ . We set for all voxels $\sigma_\tau=\sigma_\eta=0.5$ and keep default priors of the Bayes pRF toolbox for haemodynamic parameters. The quality of parameter recovery is evaluated for different amplitudes of added Gaussian noise.

| ![test](https://i.imgur.com/kiKCt6j.png) | 
| :-: |
|Figure 6 : Simulated BOLD signal for different noise levels |


By performing Bayesian Parameter Averaging (BPA) over the few voxels latent parameters estimates, we summarise the covariance across voxels. For this, the covariance matrix of latent parameters of each voxel $C_v$ are geometrically averaged to compute the BPA covariance matrix $C=(\sum_{v}C_v^{-1})^{-1}$. From this covariance matrix, we get the correlation between BPA averaged parameters. Correlation between parameters is detrimental to the model evidence and makes inference on them less precise. This issue has been briefly discussed in a previous Bayesian field study. The correlation between our BPA averaged parameters is higher than expected.

| ![test](https://i.imgur.com/Pnc1wnk.png) | 
| :-: |
|Figure 7 : Posterior receptor fields for each simulated voxels, for noise variance of 0.1 |

To evaluate the accuracy of parameter recovery for different levels of observation noise, we compute a Euclidean distance.

For a given noise level, we measure the quality of the recovery by the distance $||p_{MAP}-p_{true}||$ with $||\cdot||$ the euclidean distance in $\mathbb{R}^{|p|}$, with $p_{MAP}$ as a vector containing the maximum the posteriori (MAP) estimate of parameter $p$ for all voxels and $p_{true}$ the true value of parameter p for all voxel.

| ![test](https://i.imgur.com/BbMTZ5M.png) | 
| :-: |
|Figure 8 : Euclidean distance between estimated and true parameters |

| ![test](https://i.imgur.com/y5xIS7i.png) | 
| :-: |
|Figure 9 : Correlation plot |



### Supermodel and Model Reduction

To perform Bayesian model comparison between the four models, fitting them individually to the data is not necessary. We may only fit model n°4, which we will refer to as the supermodel, then reduce it to one of the other three models when needed. This way, we can run parameter estimation once instead of four times.

We will define a new parameter $\epsilon$ that can be "turned off" to obtain the reduced standard TD model for the distributional models. 
Similarly to the previous transformation of $\alpha$, we define $\epsilon$ such that: $\tau - \varepsilon =logit(\alpha^-)$ and $\tau + \varepsilon =logit(\alpha^+)$.

And so, the supermodel (model 4) has 7 free neuronal parameters $\mu_{\eta}$, $\sigma_{\eta}$, $\mu_{\tau}$, $\sigma_{\tau}$, $\mu_{\varepsilon}$, $\sigma_{\varepsilon}$ and $\beta$. 

To obtain an $\eta$ independent model, we will have to "switch off" the parameters $\mu_{\eta}$ and $\sigma_{\eta}$ by fixing their latent values as such:

- We reduce the distribution of latent parameter $l_\mu^\eta$ to its mean $0$. This fixes parameter $\mu_{\eta}$ value to $0$.
- We reduce the distribution of the latent parameter $l_\sigma^\eta$ to a small value $-10$, this fixes $\sigma_\eta$ is to $\sigma_{min}$. The relation between $\sigma_{min}$ and the bin size of the mesh grid $r$ ensures that all the Gaussian weights of the receptive field will be concentrated on the closest precomputed value of $\sigma_\eta$.


We switch off $\epsilon$ parameters in the same manner to obtain a standard TD model.


### Model Selection

Here we perform parameter recovery on data generated by the supermodel and compare its evidence to the reduced models evidences. The meshgrid $\Theta$ is defined over the full model 3-d parameter space by setting as previously $\mu_{min}=-2$, $\mu_{max}=2$, $\sigma_{min}=0.1$, $\sigma_{max}=1$, $R=2$, $r=0.2$ for $\tau$, $\varepsilon$ and $\eta$.

We generate twenty-seven voxels time series by taking the location parameters to be all the triplets $(\mu_\tau,\mu_\varepsilon,\mu_{\eta})$ of sets $A_{\tau}=\{-1,0,1\}$,  $A_{\varepsilon}=\{-0.5,0,0.5\}$ and  $A_{\eta}=\{-1,0,1\}$. For all voxels $\sigma_\tau=\sigma_\varepsilon=\sigma_\eta=0.5$. Similarly to before, haemodynamic parameters keep their default priors set by the toolbox.
eta=0.
The quality of parameter recovery is evaluated for different amplitudes of added Gaussian noise.

| ![test](https://i.imgur.com/nWeRiuo.jpg) | 
| :-: |
|Figure 11 : Posterior receptor fields for noise variance 0.1 |

The VBA toolbox provides a function for group Bayesian model comparison. With the reduction of both mean and width parameters, reduced model evidence is lower than the full model. These poor performances are expected: the reduced model cannot adequately fit the data generated from a full model with a lower parameter precision because the reduced model has too high precision.

We use the protected exceedance probability to compare models over all voxels. It measures how likely it is that a model is more frequent than all the others among all the voxels: $EP_i=P(r_i>r_j|y),j\neq i$ with $r_i$ the frequency of model $i$.

When recovering the parameters from data generated by the previous model n°3 $(\tau,\eta)$, the reduced model with $\epsilon$ fixed at $0$ is chosen by the Bayesian model selection, as expected.

| ![test](https://i.imgur.com/S0UJYGX.png) | 
| :-: |
|Figure 11 : Model comparison, data generated by model n°3 |

### Computation cost and optimization

One major inconvenience of using population receptive field models is the computational cost required to do parameter estimations for a single voxel. Moreover, a single fMRI scan contains many voxels that exceed the number of cores available at our cluster by three orders of magnitude. In this context, parallelisation over CPU cores is not sufficient. 

#### The $spm\_int$ optimization

By timing various parts of the estimation of a single voxel, we find that the execution of the SPM integration function, $spm\_int$, takes up 80% of the computation time. This numerical integration generates a BOLD signal from a neuronal signal $z(t)$ and a hemodynamic model.

A neuronal signal $z(t)$ a BOLD signal is calculated by solving the differential equations representing the hemodynamics of each timestep $t$. In the case of reward models, the neuronal signal only takes on values at stimuli timings (4% of the signal takes on non zero values). The $spm\_int$ function unnecessarily calculates the same solution for the differential equations for every timestep where the neuronal signal equals 0.

We modify this function to store the solution for $z(t)=0$ and only solves the differential equations for non zero values. By testing this on reward models, we get a speed-up of x$13$ compared to the normal execution of $spm\_int$. Of course, this speed-up is dependent on the sparsity of the neuronal signal. To demonstrate the usefulness of this simple modification, we simulate neuronal inputs with varying sparsity and plot the speed up factor:

| ![test](https://i.imgur.com/SBLUoGY.png) | 
| :-: |
|Figure 11 : Speed up gained against the proportion of non-zero values. For 4% non zero values, the modified $spm\_int$ is x$13$ times faster than normal execution. |

This graph is specific to the Balloon extended hemodynamic model. For more complex models (models with larger state vectors), the speedup will be lower. However, the inverse relationship between sparsity and speed up holds.

The presented modification has no downsides and always results in the same output as normal execution.

#### BOLD pre-computation 

The previous modification is effective. However, in the case of high-resolution scans, complex models, or non-sparse signals, it may not be sufficient on its own. Here we present another solution to speed up voxel-wise estimation.

We may extend the parameter space to include hemodynamic parameters. This means that we also discretise the hemodynamic parameter space, and for every combination of points, we calculate the associated BOLD signal. The resulting combined parameter space would be of higher dimension.

By learning a receptive field over this combined parameter space, we are essentially learning both the hemodynamic and neuronal parameters together. This requires us to estimate the hemodynamic parameters in terms of locations $\mu$ and widths $\Sigma$. These new $\mu$ and $\Sigma$ parameters follow the same latent transformations as the neuronal receptive field parameters. The final BOLD signal can be expressed as:

$y(t)= \beta\sum_{k\in \Theta} m(k,t)N(k : \boldsymbol\mu,\boldsymbol\Sigma)$

Where $m(k,t)$ is the BOLD signal at time $t$ for parameters $k=(\theta_{1},..,\theta_{N},\theta_{N+1},...,\theta_{N+H})$
Where $\theta_{1},..,\theta_{N}$ are the neuronal parameters and $\theta_{N+1},...,\theta_{N+H}$ are the hemodynamic parameters. $N$ and $M$ are the number of neural and hemodynamic parameters, respectively.

$m(t)$ is pre-computed for every point in the mesh-grid $\Theta$. The same pre-computed grid can be used for all voxels across all subjects. It only needs to be calculated once for every experimental course of events. When the overhead computation is completed, and for the case of reward models, we may estimate a single voxel in less than 2 seconds instead of 4 minutes with the normal execution.

This method has two drawbacks. First, it is less precise since we are discretising the hemodynamic parameters. Second, the BOLD pre-computation may become computationally intractable for high dimensional parameter spaces. Here the choice of the hemodynamic model and the parameter grid resolution becomes an important consideration. For instance, if we need to estimate a complex neuronal model with many parameters, we may want to use a simple hemodynamic model such as the canonical HRF or reduce the resolution of a specific parameter. As a rule of thumb, the choice of the resolutions and the hemodynamic model should be such that:

$\prod_{i=1}^{N+H} r_i << n_{voxels}*n_{subjects}$

Where $r_i$ is the resolution of the $i$th parameter.

By dividing this product over the number of available cores, we may also roughly estimate its time to precompute all the BOLD signals.

## Discussion

### Summary

### Limitations

### Perspectives


## References
1.	Kumar, S. & Penny, W. Estimating neural response functions from fMRI. Front. Neuroinformatics 8, (2014).
2.	Friston, K., Mattout, J., Trujillo-Barreto, N., Ashburner, J. & Penny, W. Variational free energy and the Laplace approximation. NeuroImage 34, 220–234 (2007).
3.	Zeidman, P., Silson, E. H., Schwarzkopf, D. S., Baker, C. I. & Penny, W. Bayesian population receptive field modelling. NeuroImage 180, 173–187 (2018).
4.	Meder, D. et al. Ergodicity-breaking reveals time optimal economic behavior in humans. ArXiv190604652 Econ Q-Fin (2019).
5.	Schultz, W. Predictive Reward Signal of Dopamine Neurons. J. Neurophysiol. 80, 1–27 (1998).
6.	Ludvig, E. A., Sutton, R. S. & Kehoe, E. J. Evaluating the TD model of classical conditioning. Learn. Behav. 40, 305–319 (2012).
7.	A distributional code for value in dopamine-based reinforcement learning | Nature.
8.	Ludvig, E. A., Sutton, R. S. & Kehoe, E. J. Stimulus Representation and the Timing of Reward-Prediction Errors in Models of the Dopamine System. Neural Comput. 20, 3034–3054 (2008).
9.	Stephan, K.E., Weiskopf, N., Drysdale, P.M., Robinson, P.A., Friston, K.J., 2007.Comparing Hemodynamic Models with DCM.


# Cognitive field mapping of fMRI data 

## Abstract

The analysis of functional neuroimaging data and its integration with computational models of cognition typically depends on fitting models to behaviour and regressing the latent variables of those models onto the neural data. Here we develop a technique that directly fits computational cognitive models to the neuroimaging data, affording model comparison and the visualisation of latent model parameters over neuroanatomical space. Our approach builds on a probabilistic framework designed for mapping population receptive fields (pRFs) in sensory cortices. The advance we offer is modest but useful. First, we adapt the modelling framework from the sensory to the cognitive domain, providing a primer in how this can be achieved generically for cognitive models. Second, we lower the computation time per voxel to a level that makes the technique practically feasible for larger regions of interest and volumes of data. Third, we offer extensive simulations for model comparison and parameter recovery to demonstrate the robustness and specificity of this technique. This approach, which we generically call cognitive field mapping, thus demonstrates its potential in answering important questions about the topographic principles underlying the encoding of cognitive computational variables. 

## Introduction

Deploying computational models of the mind, brain, and behaviour to analyse neuroimaging data is now a powerful and widely used strategy in the cognitive and computational neurosciences. Since the advent of this strategy, the methods used have typically depended on fitting computational models to behavioural data first and then regressing the latent variables of these fitted models onto the neural data. A different approach has however been pioneered in sensory neuroimaging. Kumar and Penny (2014) developed a Bayesian method that allows the estimation of neural response functions from fMRI data. The two building blocks of the method are a neural response function (NRF) which maps stimuli characteristics to neuronal population responses, and a hemodynamic response function (HRF) which maps neural responses to fMRI data. The method was initially used for modelling the receptive field structure of populations of neurons in auditory cortices and their topographic mapping (Kumar and Penny, 2014). This technique was then subsequently generalised to model visual receptive fields (Zeidman ref) and released as a toolbox known as BayespRF (Zeidman ref) available for the Statistical Parametric Mapping package (ref). Henceforth we refer to this method in the sensory domain as the pRF method. 

Here, we are interested in developing the pRFmethod and adapting it to estimate cognitive computational models flexibly. To do so we will describe an analogy between sensory and cognitive domains. A sensory space can be represented via its dimensions, where for a retinotopic space it may have two dimensions corresponding to its x and y coordinates. The receptive field is a distribution of weights within this space, that describe the relative weighting necessary to map from the sensory input to the neuronal response. This spatial distribution of weights may compactly be approximated via a function such as a gaussian or difference of gaussian, taking parameters that control its location and dispersion within sensory space. Likewise a cognitive space can be represented via its dimension, where for a reward model it may have two dimensions corresponding to its learning rate and its risk tolerance parameter. The cognitive field is then a spatial distribution of weights within this cognitive space that describe the relative weighting necessary to map from cognitive inputs to the neuronal response. Likewise the spatial distribution of weights may also be compactly approximated by a gaussian, taking parameters that control its location and dispersion within cognitive space.     

Effectively via this approach, one is ultimately modelling the mean and variances of cognitive models in accounting for the neural data. This method can be generalised for any application that requires fitting a large number of generative models. And via different observational models one can apply this technique to modelling other modalities of data such as calcium imaging. In this paper we stick to fMRI data, and demonstrate the technique by applying it to reward prediction error models, arguably the most ubiquitous computational model in cognitive computational neuroimaging.  

## Formal definition of a Bayesian Receptive Field Model

Let $s$ be a cognitive model with $n$ free parameters $\theta_N=p_{1},...,p_n\in\mathbb{R}$. We consider that these parameters have support on the real number line. An experimental session is comprised of $K$ stimulus onset events happening at times $T = {t_1,...,t_K}$, with events characteristics $C={c_1,...,c_K}$. For a given set of free parameters $\theta_N$ and a given experimental course of events $(T,C)$ as inputs, the cognitive model outputs a temporal signal $s_{\theta_N}(t)$. 

The neuronal response function is obtained by using a receptive field (a Gaussian function) over the parameter space of the cognitive model. From $\theta_N=p_{1},...,p_n$ we specify the receptive field location $\boldsymbol\mu=(\mu_1,...,\mu_n)$ and width $\boldsymbol\Sigma=diag(\sigma_1^2,...,\sigma_n^2)$. We use a precomputed subset of the parameter space $\Theta$ over which the receptive field is learned. The neuronal response function is expressed in the following equation:

$z(t)=\beta\sum_{k\in \Theta} s(k,t)N(k : \boldsymbol\mu,\boldsymbol\Sigma)$
with $N(k:\boldsymbol\mu,\boldsymbol\Sigma)=\frac{1}{\sqrt{(2\pi)^n|\boldsymbol\Sigma|}} \exp\left(-\frac{1}{2}(k-\boldsymbol\mu)^T{\boldsymbol\Sigma}^{-1}(k-\boldsymbol\mu)\right)$


The obtained neuronal signal $z(t)$ is then fed into a hemodynamic model with parameters $\theta_H$ that generates a synthetic BOLD signal. 

The posterior distribution of free parameters $\theta = (\theta_N,\theta_H)$ is then found by fitting the synthetic BOLD signals to the measured BOLD signals for each fMRI voxels using a variational Bayesian approach. Using a variational Bayesian approach (such as variational Laplace) is beneficial as it allows us to quantitatively compare cognitive models.

### Variational Bayes and Model Comparison

Let $s$ be a generative model of data $y$, defined by the likelihood function $p(y|\theta,s)$ and the prior density $p(\theta|s)$ on model parameters $\theta$. In brief, variational Bayes is an iterative scheme that indirectly optimizes an approximation to both the model evidence $p(y∣s)$ and the posterior density $p(\theta∣y,s)$. The key trick is to decompose the log model evidence into:

$log \, p(y|s)=F(q)+D_{KL}(q(\theta)||p(\theta|y,s))$

Where $D_{KL}$ is the Kullback–Leibler divergence, $q$ is the approximate distribution and $F$ is the variational free energy:
$F(q)=log \, p(y|s)-D_{KL}(q(\theta)||p(\theta|y,m))$

By maximizing $F(q)$ we are essentially minimizing the divergence between the posterior distribution $p(\theta|y,s)$ and the approximate distribution $q(\theta)$. This maximisation is done by way of gradient ascent.



Given two concurrent models $s_1$ and $s_2$, the factor of their minimized free energy $F_1^*$ and $F_2^*$ approximates the Bayes Factor between models :

$B_{12}=\frac{p(y|s_1)}{p(y|s_2)}\approx exp(F_1^*-F_2^*)$

This allows us to perform Bayesian model comparison and selection by choosing the model with the strongest evidence.

### Variational Laplace and neural Response function priors

The variational Laplace (VL) method uses simplifying assumptions on the approximate distribution $q(\theta)$, which render the variational Bayes algorithm computationally and statistically efficient [cite]. One of these assumptions, the Laplace approximation, requires the model's prior to follow a Gaussian distribution. As such we will be expressing the receptive field locations $\mu$ and widths $\Sigma$ as latent Gaussian parameters.

We wish the priors to be as uninformative as possible, so we maximize their entropy by constraining them to uniforms distributions on a range of real numbers. Thus, for each parameter, location $\mu$ and width $\sigma$ are expressed in terms of latent parameters $l_{\mu}$ and $l_{\sigma}$ with normal priors:

$\mu = (\mu_{max}-\mu_{min})*ncdf(l_{\mu})+\mu_{min}$
$\sigma = (\sigma_{max}-\sigma_{min})*ncdf(l_{\sigma})+\sigma_{min}$
Where $ncdf$ is the cumulative density function for the univariate normal distribution, and $\mu_{min},\mu_{max}, \sigma_{min},\sigma_{max}$ are the boundaries where the parameters can take their values. $l_{\mu}$ and $l_{\sigma}$ are the latent parameters with normal priors.

The described transformations turn the normal distribution of a latent parameter into a uniform distribution of the parameter across its support range.

$l_{\mu}\sim\mathcal{N}(0, 1) \implies \mu\sim\mathcal{U}(\mu_{min}, \mu_{max})$


$l_{\sigma}\sim\mathcal{N}(0, 1) \implies \sigma\sim\mathcal{U}(\sigma_{min}, \sigma_{max})$

| ![test](https://i.imgur.com/HQzYphd.png) | 
| :-: |
|**Figure 1** : Transformation of latent parameter $l_\mu$ distribution to the parameter $\mu$ distribution. The prior latent distribution is gaussian with a mean 0 and standard deviation of 1. The posterior latent distribution is gaussian with a mean 0.7 and a standard deviation of 0.1. The posterior distribution is not uniform, but approximatively gaussian around a value, here 2, when the latent parameter standard deviation is low enough.|


The scaling parameter $\beta$ is constrained to be positive using a normal latent parameter $l_{\beta}$ and through the transformation $\beta=exp(l_{\beta})$ .

| ![test](https://i.imgur.com/1arOXk3.png) | 
| :-: |
|**Figure 2** : Transformation of latent parameter $l_\beta$ distribution to the parameter $\beta$ distribution. The prior latent distribution is gaussian with mean 0 and a standard deviation 1. The posterior latent distribution is gaussian with mean 0.7 and standard deviation 0.1. The posterior distribution is not uniform, but approximatively gaussian around a value, here 2, when the latent parameter standard deviation is low enough.|




### Use of a Precomputed Meshgrid



During the fitting, using the VL algorithm, the computing of $s_{\theta_N}$ is performed a large number of times for different sets of parameters. This process can be computationally extensive, especially in the context of models of reward. To avoid excessive computation, the values of $s_{\theta_N}$ over the experimental session is precomputed on a discretised subset of the model parameter space. 

A finite support $[p_{min},p_{max}]$ and a bin size $r$ are used to define a discrete subset of values for each parameter $p$. The choice of the bin size depends on available computing power and the desired precision. The choice of the boundaries depends on a reasonable hypothesis given existing data. A meshgrid $\Theta$ over the neuronal parameter space is obtained, over which $s_{\theta_N}$ is computed at every point. 

| ![test](https://i.imgur.com/Rt3i1XE.png) | 
| :-: |
|**Figure 3** : Meshgrid parametrization for a bidimensionnal parameter space |


The precomputed meshgrid $\Theta$ has to be chosen such that the centre of the population receptive field (pRF) falls inside its bounds. As a rule of thumb, we defined it such that at least 95% of prior density fall inside the meshgrid area, and such that at least a grid point is contained around the centre of the receptive field.  This gives for each parameter $p$:

- $[p_{min},p_{max}]  \supset[\mu_{min}-2\sigma_{max},\mu_{max}+2\sigma_{max}]$ 
- $r \leq 2\sigma_{min}$



### Observation model

It's not possible to observe neuronal activity directly. The neuronal signal is inferred from a secondary phenomenon which also needs to be modelled, in the case of fMRI it's the hemodynamic response.

The default Hemodynamic model in the BayespRF toolbox, developed by Kumar and Penny [cite], is an extended Balloon Model. Most parameters in this model are fixed, based on empirical measurements, leaving three parameters to be estimated:  the transit time $\tau$, the rate of signal decay $\kappa$, and the ratio of intra-to extra-vascular signal $\epsilon$. They are also expressed as the transformation of latent parameters similar to $\beta$. We will be using this model for what follows.


## Bayesian Receptive Field Models of Reward


Here we apply the previous formal development for the purpose of modelling the reward system of the brain. We will showcase Bayesian field modelling on simulated data. We will perform parameter recovery, model comparison and model recovery on the simulated data. We will be using a modified temporal difference learning algorithm as the cognitive model. This is a particular class of reinforcement learning algorithm that currently best explain conditional learning behaviour and especially neuronal activity in dopaminergic areas of the brain. The phasic activity of dopamine neurons of the midbrain (VTA & SN) has been linked to the reward prediction error of reinforcement learning algorithms [cite]. 


### Experiment setting:

We consider the following experimental setting: a passive learning task from a broader experiment done at the DRCMR in 2017 [2]. A subject is given a hypothetical wealth of 1000DKK and is then shown a set of fractal pictures. Each fractal picture has a fixed effect on the subject's wealth, which the subject is asked to learn by passive observation.  A session consists of 150 trials with the events in a single trial consisting of:
- <ins>beginning of the trial</ins>, a cue for the start of the trial
- <ins>fractal onset</ins>, a fractal picture
- <ins>wealth variation</ins>, a visible update of the subject's wealth

| ![test](https://i.imgur.com/wNPSZIw.png) | 
| :-: |
|**Figure 5** : Evolution of a subject wealth over time during an experiment session. |!

### Temporal difference learning:

The learning process of a subject can be modelled by a temporal difference (TD) algorithm. 

| ![test](https://i.imgur.com/1VKhebU.png) | 
| :-: |
|**Figure 4** : Computing the reward prediction error with the TD learning algorithm |


Learning in this context is seen as a progressive evolution of expected utility $U$ over trials, caused by temporal relations between stimuli. We are especially interested in the reward prediction error (RPE) signal of this algorithm, as it has been directly linked with dopaminergic activity. Through this experiment, our model is essentially learning an associative strength vector $W$ which represents the 'value function' in this context.

The utility here represents the subjective value given by a subject to the wealth variation he perceives during the experimental session. This can be modelled as a transition from an amount of money $C_t$ to an amount $C_{t+1}$ by an isoelastic utility function $iso$.

$U_{t+1}=iso(C_{t+1},\eta)-iso(C_{t},\eta)$

with $iso(c,\eta)=\frac{c^{1-\eta}}{1-\eta}$ for $\eta \neq1$ 
and $iso(c,1)=ln(c)$

$\eta$ is termed the aversion degree.

Here we define a standard TD algorithm and use it as the neuronal model, with the reward prediction error as the output signal. This model takes into account the value of rewards being delivered at each time step:

$W_{t+1}= W_{t}+\alpha E(t)Z_t$
$E(t)=U_{t+1}+\gamma W_tX_{t+1}-W_tX_t$
$Z_{t+1}=\gamma \lambda Z_t+X_t$

Where $W_t$ is the associative strength vector, $\alpha$ is the learning rate, $E$ is the reward prediction error, $Z_t$ is the eligibility-trace vector which acts as a sort of backwards memory for the TD algorithm [cite], $U_t$ is the utility, $X_t$ is the complete serial compound (CSC) feature vector which represents the stimulus [cite], $\gamma$ and $\lambda$ are decay rates.

| ![test](https://i.imgur.com/u854FPB.png) | 
| :-: |
|**Figure 4** : Computing the reward prediction error $E(t)$ with the TD learning algorithm |



### Cognitive models:

We are interested in two parameters: the learning rate $\alpha$ and the risk aversion degree $\eta$. We fix all other parameters to standard values:
- Discount factor $\gamma=0.97$
- Eligibility trace decay rate $\lambda=1$

We are also interested in comparing different flavours of the proposed model. 


We are interested in comparing different variations of the proposed TD model. For instance, it has been demonstrated that the RPE signals in mice VTA neurons actually follow a distributional code [cite]. Under this class of learning model, known as distributional reinforcement learning, we have two learning parameters $\alpha^+$ and $\alpha^-$ for positive and negative RPEs respectively. Also, we are interested in knowing if including the aversion degree (the utility function parameter) as a parameter is worth the extra model complexity. 

With these two considerations, we have a total of 4 possible cognitive models outlined in this table:
| ![test](https://i.imgur.com/N2SwIKk.png) | 
| :-: |
|Figure 10 : Cognitive models of reward |


To compare evidence for the data under different models, we use the free energy. For each voxel, the VL algorithm returns the free energy of each model given the data observed. We use VBA (Variational Bayesian Analysis) toolbox to perform a random-effects RFX analysis using a $N_{voxels} \times N_{models}$ matrix of free energies. This is usually used for comparison across subjects but can also be applied across a population of voxels within a single subject’s brain. The analysis returns the expected probabilities of each model which measures how likely it is that the model is more frequent than all other models in the comparison set.


### Model n°3 Specification

We will begin by specifying model 3. The learning rate $\alpha$ is between 0 and 1, so we transform the parameter in one varying on the full real line with a logit transformation , $\tau = logit(\alpha) = log(\frac{\alpha}{1-\alpha})$. This gives us a proper bidimensional learning space. The meshgrid $\Theta$ is defined by setting: $\mu_{min}=-2,\mu_{max}=2,\sigma_{min}=0.1,\sigma_{max}=1,R=2,r=0.2$ for both $\tau$ and $\eta$.

The reward prediction error signal is composed of spikes of various intensity, occurring at stimuli timings. Elsewhere it takes on a null value. The neural response function is then:

$z(t)=\beta\sum_{j=1}^N \delta(t-t_j)\sum_{k \in \Theta} E_k(t)N(k:\boldsymbol\mu,\boldsymbol\Sigma)$

Where $t_j$ is the timing of the j-th stimulus, $E_k(t)$ is the reward prediction error signal for parameters $k=(\tau,\eta)$ at time $t$, $\boldsymbol\mu=(\mu_{\tau},\mu_{\eta})$ is the location of the receptive field and $\boldsymbol\Sigma=diag(\sigma_{\tau}^2,\sigma_{\eta}^2)$   is its covariance matrix. 

| ![test](https://i.imgur.com/qxFi4Y5.png) | 
| :-: |
|**Figure 6** : Prior & Posterior probability density |

In total, we have 8 parameters: five neuronal plus three haemodynamic. They are all defined through transformations of latent parameters with normal priors. 



### Parameter Estimation and Recovery

We generate data with the model for a set of parameters and add various levels of noise to the BOLD signal before trying to recover the parameters. We generate nine voxels time series by taking the location parameters to be all the pairs $(\mu_\tau,\mu_{\eta})$ of sets $A_{\tau}=\{-1,0,1\}$ and $A_{\eta}=\{-1,0,1\}$ . We set for all voxels $\sigma_\tau=\sigma_\eta=0.5$ and keep default priors of the Bayes pRF toolbox for haemodynamic parameters. The quality of parameter recovery is evaluated for different amplitudes of added Gaussian noise.

| ![test](https://i.imgur.com/kiKCt6j.png) | 
| :-: |
|Figure 6 : Simulated BOLD signal for different noise level |


By performing Bayesian Parameter Averaging (BPA) over the few voxels latent parameters estimates, we obtain a summary of the covariance across voxels. For this, the covariance matrix of latent parameters of each voxel $C_v$ are geometrically averaged to compute the BPA covariance matrix $C=(\sum_{v}C_v^{-1})^{-1}$. From this covariance matrix, we get the correlation between BPA averaged parameters. Correlation between parameters is detrimental to the model evidence and makes inference on them less precise. This issue has been briefly discussed in a previous Bayesian field study. The correlation between our BPA averaged parameters is higher than expected.

| ![test](https://i.imgur.com/Pnc1wnk.png) | 
| :-: |
|Figure 7 : Posterior receptor fields for each simulated voxels, for noise variance of 0.1 |

To evaluate the accuracy of parameter recovery for different levels of observation noise we use a euclidean distance.

For a given noise level, we measure the quality of the recovery by the distance $||p_{MAP}-p_{true}||$ with $||\cdot||$ the euclidean distance in $\mathbb{R}^{|p|}$, with $p_{MAP}$ as a vector containing the maximum the posteriori (MAP) estimate of parameter $p$ for all voxels and $p_{true}$ the true value of parameter p for all voxel.

| ![test](https://i.imgur.com/BbMTZ5M.png) | 
| :-: |
|Figure 8 : Euclidean distance between estimated and true parameters |

| ![test](https://i.imgur.com/y5xIS7i.png) | 
| :-: |
|Figure 9 : Correlation plot |



### Supermodel and Model Reduction

To perform Bayesian model comparison between the four models, fitting them individually to the data is not necessary. We may only fit model n°4, which we will refer to as the supermodel, then reduce it to one of the other three models when needed. This way, we're able to run parameter estimation once instead of four times.

For the distributional models, we will define a new parameter $\epsilon$ that can be "turned off" to obtain the reduced standard TD model. 
Similarly to the previous transformation of $\alpha$, we define $\epsilon$ such that: $\tau - \varepsilon =logit(\alpha^-)$ and $\tau + \varepsilon =logit(\alpha^+)$.

And so, the supermodel (model 4) has 7 free neuronal parameters $\mu_{\eta}$, $\sigma_{\eta}$, $\mu_{\tau}$, $\sigma_{\tau}$, $\mu_{\varepsilon}$, $\sigma_{\varepsilon}$ and $\beta$. 

To obtain an $\eta$ independent model, we will have to "switch off" the parameters $\mu_{\eta}$ and $\sigma_{\eta}$ by fixing their latent values as such:

- We reduce the distribution of latent parameter $l_\mu^\eta$ to its mean $0$. This fixes parameter $\mu_{\eta}$ value to $0$.
- We reduce the distribution of the latent parameter $l_\sigma^\eta$ to a small value $-10$, this fixes $\sigma_\eta$ is to $\sigma_{min}$. The relation between $\sigma_{min}$ and the bin size of the mesh grid $r$ ensures that all the Gaussian weights of the receptive field will be concentrated on the closest precomputed value of $\sigma_\eta$.


We switch off $\epsilon$ parameters in the same manner to obtain a standard TD model.


### Model Selection

Here we perform parameter recovery on data generated by the supermodel and compare its evidence to the reduced models evidences. The meshgrid $\Theta$ is defined over the full model 3-d parameter space by setting as previously $\mu_{min}=-2$, $\mu_{max}=2$, $\sigma_{min}=0.1$, $\sigma_{max}=1$, $R=2$, $r=0.2$ for $\tau$, $\varepsilon$ and $\eta$.

We generate twenty-seven voxels time series by taking the location parameters to be all the triplets $(\mu_\tau,\mu_\varepsilon,\mu_{\eta})$ of sets $A_{\tau}=\{-1,0,1\}$,  $A_{\varepsilon}=\{-0.5,0,0.5\}$ and  $A_{\eta}=\{-1,0,1\}$. For all voxels $\sigma_\tau=\sigma_\varepsilon=\sigma_\eta=0.5$. Similarly to before, haemodynamic parameters keep their default priors set by the toolbox.

The quality of parameter recovery is evaluated for different amplitudes of added Gaussian noise.

| ![test](https://i.imgur.com/nWeRiuo.jpg) | 
| :-: |
|Figure 11 : Posterior receptor fields for noise variance 0.1 |

The VBA toolbox provides a function for group Bayesian model comparison. With the reduction of both mean and width parameters, reduced model evidence is lower compared to the full model. These poor performances are excepted: by reducing a parameter's width, the reduced model can't accurately fit the data generated from the full model with a greater parameter width.

We use the protected exceedance probability to compare models over all voxels. It measures how likely it is that a model is more frequent than all the others among all the voxels: $EP_i=P(r_i>r_j|y),j\neq i$ with $r_i$ the frequency of model $i$.

When recovering the parameters from data generated by the previous model n°3 $(\tau,\eta)$, the reduced model with $\epsilon$ fixed at $0$ is chosen by the Bayesian model selection, as expected.

| ![test](https://i.imgur.com/S0UJYGX.png) | 
| :-: |
|Figure 11 : Model comparison, data generated by model n°3 |

### Computation cost and optimization

One major inconvenience of using population receptive field models is the computational cost required to do parameter estimations for a single voxel. Moreover, a single fMRI scan contains a large number of voxels that exceed the number of cores available at our cluster by three orders of magnitude. In this context, parallelization over CPU cores is not sufficient. 

#### The $spm\_int$ optimization

By timing various parts of the estimation of a single voxel, we find that the execution of the SPM integration function, $spm\_int$, takes up 80% of the computation time. This numerical integration is essentially what generates a BOLD signal from a neuronal signal $z(t)$ and a hemodynamic model.

For a neuronal signal $z(t)$ a BOLD signal is calculated by solving the differential equations representing the hemodynamics for each timestep $t$. In the case of models of reward, the neuronal signal only takes on values at stimuli timings (4% of the signal takes on non zero values). The $spm\_int$ function unnecessarily calculates the same solution for the differential equations for every timestep where the neuronal signal is equal to 0.

We modify this function such that it stores the solution for $z(t)=0$ and only solves the differential equations for non zero values. By testing this on models of reward, we get a speed-up of x$13$ compared to the normal execution of $spm\_int$. Of course, this speed-up is dependent on the sparsity of the neuronal signal. To demonstrate the usefulness of this simple modification we simulate neuronal inputs with varying sparsity and plot the speed up factor:

| ![test](https://i.imgur.com/SBLUoGY.png) | 
| :-: |
|Figure 11 : Speed up gained against ercentage of non-zero values. For 4% non zero values, the modified $spm\_int$ is x$13$ times faster than normal execution. |

This graph is specific to the Balloon extended hemodynamic model. For more complex models (models with larger state vectors) the speed up will be lower, however the inverse relationship between sparsity and speed up holds.

The presented modification has no downsides and always results in the same exact output as normal execution.

#### BOLD pre-computation 

The previous modification is effective, however, in the case of high-resolution scans, complex models, or non-sparse signals it may not be sufficient. Here we present another solution to speed up voxel-wise estimation.

We may extend the parameter space to include hemodynamic parameters. This means that we also discretize the hemodynamic parameter space and for every combination of points we calculate the associated BOLD signal. The resulting combined parameter space would be of higher dimension.

By learning a receptive field over this combined parameter space, we are essentially learning both the hemodynamic and neuronal parameters together. This requires us to estimate the hemodynamic parameters in terms of locations $\mu$ and widths $\Sigma$. These new $\mu$ and $\Sigma$ parameters follow the same latent transformations as the neuronal receptive field parameters. The final BOLD signal can be expressed as:

$y(t)= \beta\sum_{k\in \Theta} m(k,t)N(k : \boldsymbol\mu,\boldsymbol\Sigma)$

Where $m(k,t)$ is the BOLD signal at time $t$ for parameters $k=(\theta_{1},..,\theta_{N},\theta_{N+1},...,\theta_{N+H})$
Where $\theta_{1},..,\theta_{N}$ are the neuronal parameters and $\theta_{N+1},...,\theta_{N+H}$ are the hemodynamic parameters. $N$ and $M$ are the number of neural and hemodynamic parameters, respectively.

$m(t)$ is pre-computed for every point in the mesh-grid $\Theta$. The same pre-computed grid can be used for all voxels across all subjects. It only needs to be calculated once for every experimental course of events. When the overhead computation is completed, and for the case of models of reward, we may estimate a single voxel in less than 2 seconds as opposed to 4 minutes with the normal execution.

This method has two drawbacks, first, it is less precise since we are discretizing the hemodynamic parameters, second, the BOLD pre-computation may become computationally intractable for high dimensional parameter spaces. Here the choice of the hemodynamic model and the parameter grid resolution becomes an important consideration. For instance, if we need to estimate a complex neuronal model with many parameters we may want to use a simple hemodynamic model such as the canonical HRF, or we may reduce the resolution of a certain parameter. As a rule of thumb, the choice of the resolutions and the hemodynamic model should be such that:

$\prod_{i=1}^{N+H} r_i << n_{voxels}*n_{subjects}$

Where $r_i$ is the resolution of the $i$th parameter.

By dividing this product over the number of available cores we may also roughly estimate the time it would take to precompute all the BOLD signals.













## Discussion

### Summary

### Limitations




### Perspectives


## References
1.	Kumar, S. & Penny, W. Estimating neural response functions from fMRI. Front. Neuroinformatics 8, (2014).
2.	Friston, K., Mattout, J., Trujillo-Barreto, N., Ashburner, J. & Penny, W. Variational free energy and the Laplace approximation. NeuroImage 34, 220–234 (2007).
3.	
, P., Silson, E. H., Schwarzkopf, D. S., Baker, C. I. & Penny, W. Bayesian population receptive field modelling. NeuroImage 180, 173–187 (2018).
4.	Meder, D. et al. Ergodicity-breaking reveals time optimal economic behavior in humans. ArXiv190604652 Econ Q-Fin (2019).
5.	Schultz, W. Predictive Reward Signal of Dopamine Neurons. J. Neurophysiol. 80, 1–27 (1998).
6.	Ludvig, E. A., Sutton, R. S. & Kehoe, E. J. Evaluating the TD model of classical conditioning. Learn. Behav. 40, 305–319 (2012).
7.	A distributional code for value in dopamine-based reinforcement learning | Nature.
8.	Ludvig, E. A., Sutton, R. S. & Kehoe, E. J. Stimulus Representation and the Timing of Reward-Prediction Errors in Models of the Dopamine System. Neural Comput. 20, 3034–3054 (2008).
9.	Stephan, K.E., Weiskopf, N., Drysdale, P.M., Robinson, P.A., Friston, K.J., 2007.Comparing Hemodynamic Models with DCM.