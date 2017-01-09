# Bivariate Survival Association.
![inline fill](imagesForPresentation/semmel03.jpg)![inline fill](imagesForPresentation/220px-Florence_Nightingale_three_quarter_length.jpg)

---

# Types of survival data:

- Yes or No
- Time to event
- Time to event with censoring


---

## Earlier developments

- Fisher, 1931 and Hald, 1949 addressed _censoring_

- Kaplan and Meier, 1958, Product-Limit estimator of the survival function.

- Cox, 1972 introduced proprotional hazard regression model

---

## Defintion of hazard rate

- Cox utilized a concept of _hazard_:
$$\lambda(s|y) =\frac{f(s|y)}{\mathcal{S}(s|y)} = \frac{\partial}{\partial s}\left\{  -log \mathcal{S}(s|y)  \right\}$$

- Basu, 1971 defined a _bivariate failure rate_: 
$$\lambda(s,t) = \frac{f(s,t)}{\mathcal{S}(s,t)} = \frac{\partial}{\partial s}\left\{  -log \mathcal{S}(s,t)  \right\}\frac{\partial}{\partial t}\left\{  -log \mathcal{S}(s,t)  \right\} - \frac{\partial^2}{\partial s\partial t}\left\{  -log \mathcal{S}(s,t)  \right\}$$

---

# Three Main Approaches of Measuring Association


 - Modeling **bivariate hazard ratio**
 - Estimating **bivariate survival surface**
 - Utilizing **marginal cumulative distribution functions**
   - **Probability Scale Residuals** approach

---

## Probability Scale Residuals

- Define _Probability Scale Residuals_(PSR) for uncensored, censored, and _current status data_
- Show that _Speaman_ correlation is equivalent to correlation of PSR
- Define future direction of reseach:
  - Compare performance of PSR with other residuals
  - Designing a test of association under non-informative censoring
  - Developing a test of association for partial PSR  and for conditional in the presence of censoring and for current status data

---

## Modeling bivariate hazard ratio
**Clayton, 1978**:

- Easy to compute
- Expressible as a constant ratio of age-specific rates
- Symmetrical in two variables

$$\theta = \frac{\lambda_s(s_0|t=t_0)}{\lambda_s(s_0|t>t_0)} = \frac{\lambda_t(t_0|s=s_0)}{\lambda_t(t_0|s>s_0)}$$

---

## Modeling bivariate hazard ratio
**Clayton, 1978**:

$$\theta = \frac{\lambda_s(s_0|t=t_0)}{\lambda_s(s_0|t>t_0)} = \frac{\lambda_t(t_0|s=s_0)}{\lambda_t(t_0|s>s_0)}$$

**Interpretation**: how much more likely for a son to have an event at time $$s_0$$ given that his father had an event at time $$t_0$$ compared to a son whose father survived until time $$t_0$$ without the event.

---

## Modeling bivariate hazard ratio
**Clayton, 1978**:

$$f(s,t) = \frac{\theta a'(s)b'(t)}{[1+(\theta - 1)(a(s) + b(t))]^{2+\frac{1}{\theta-1}}},~~~a(0)=0,~~~b(0)=0$$

where $$a(s)$$ and $$b(t)$$ are unknown non-decreasing functions.

**Independence** between fathers and sons would mean $$\theta = 1$$, **positive association** $$\theta > 1$$ and **negative association** $$\theta < 1$$

---

Oakes, 1982


---

Oakes, 1989



---

Fan, Hsu, and Prentice 2000

---




<!-- ---

---

# Sampling method:
### Gibbs + (Metropolis-Hastings)

$$\theta_j^{t} |\alpha^{t-1},\beta^{t-1} \sim Beta(R_j + \alpha^{t-1},~~ N_j - R_j + \beta^{t-1})$$
$$\alpha^{t}| \pmb{\theta}^{t}, \beta^{t-1}  \propto \left[\prod_{j=1}^{N_c}  \frac{\left(\theta_j^{t}\right)^{\alpha-1} \Gamma(\alpha^{t-1} + \beta^{t-1})}{\Gamma(\alpha^{t-1})}\right]  \cdot \left(\alpha^{t-1}\right)^{a-1} e^{-b\alpha^{t-1}}$$
$$\beta^{t}| \pmb{\theta}^{t}, \alpha^{t}  \propto \left[\prod_{j=1}^{N_c}  \frac{(1 - \theta_j^{t})^{\beta-1} \Gamma(\alpha^{t} + \beta^{t-1})}{\Gamma(\beta^{t-1})}\right]  \cdot \left(\beta^{t-1}\right)^{a-1}  e^{-b\beta^{t-1}}$$



---

# Jumping distributions
$$\alpha^*\sim Gamma(0.1\alpha^{t},  0.1)$$
$$\beta^*\sim Gamma(0.1\beta^{t},  0.1)$$

# Acceptance ratio for $$\alpha$$
$$r = \frac{ P(\alpha^*| \pmb{\theta}^{t-1}, \beta^{t-1})/ Gamma(\alpha^*, 0.1\alpha^{t-1},  0.1)}{  P(\alpha^{t-1}| \pmb{\theta}^{t-1}, \beta^{t-1})/ Gamma(\alpha^{t-1}, 0.1\alpha^{*},  0.1) }$$

---

#Sampling problem

$$l(\alpha^{t}| \pmb{\theta}^{t}, \beta^{t-1})\propto $$
$$\propto \sum_{j=1}^{N_c} \left[    (\alpha^{t-1}-1)  log(\theta_j^{t}) + log (\Gamma(\alpha^{t-1} + \beta^{t-1})) - log(\Gamma(\alpha^{t-1}))\right]  +  (a-1)log\left(\alpha^{t-1}\right)  -  b\alpha^{t-1}$$

$$l(\beta^{t}| \pmb{\theta}^{t}, \alpha^{t})\propto$$
$$\propto \sum_{j=1}^{N_c} \left[(\beta^{t-1}-1) log(1 - \theta_j^{t})  + log(\Gamma(\alpha^{t} + \beta^{t-1})) - log(\Gamma(\beta^{t-1}))\right]   +  (a-1)log\left(\beta^{t-1}\right) -b\beta^{t-1}$$

---

# Acceptance ratio modification for $$\alpha$$:

When any $$\theta_j=0$$, instead of:

$$r = \frac{ P(\alpha^*| \pmb{\theta}^{t-1}, \beta^{t-1})/ Gamma(\alpha^*, 0.1\alpha^{t-1},  0.1)}{  P(\alpha^{t-1}| \pmb{\theta}^{t-1}, \beta^{t-1})/ Gamma(\alpha^{t-1}, 0.1\alpha^{*},  0.1) }$$

we compute:

$$r = \frac{ Gamma(\alpha^{t-1}, ~0.1\alpha^{*},  0.1)}{   Gamma(\alpha^{*}, ~0.1\alpha^{t-1},  0.1) }$$

Similar modification is made for $$\beta$$ when $$\theta_j = 1$$

---



## Likelihood estimates

![inline](Likelihood.pdf)

---

# Posterior for $$\theta_j$$ for model II(a)
![inline](PostDistrIIa.pdf)

---

# Posterior for $$\theta_j$$ for model II(b)
![inline](PostDistrIIb.pdf)

---

# Posterior predictive  distributions for model II(a)
![inline](PostPredIIa.pdf)

---

# Posterior predictive  distributions for model II(b)
![inline](PostPredIIb.pdf)

---

## Convergence based on $$\hat{R} = \sqrt{\frac{\hat{var}^+ (\psi|y)}{W}}$$
![inline](covergenceRIIa.pdf)
![inline](covergenceRIIb.pdf)

---

###Trace plots  for model II(a) and  model II(b)
![inline](TraceIIa.pdf)![inline](TraceIIb.pdf)

---

##What does it take to win the OSCAR:
![inline](CheckHypoth.pdf)

---

![left](easternpromiseslede.jpg)

#Thank you:

- **Prof. Leena Choi**: for the interesting course and for your support of this project
- **My classmates**: for your patience, understanding, and support
- **Vigo Mortensen**: for inspiration
 -->
