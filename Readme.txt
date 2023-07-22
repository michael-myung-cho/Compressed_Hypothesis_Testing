Title: Compressed Hypothesis Testing: To Mix or Not to Mix?


A. INTRODUCE
In this paper, we study the problem of determining $k$ anomalous random variables that have different probability distributions from the rest $(n-k)$ random variables. Instead of sampling each individual random variable separately as in the conventional hypothesis testing, we propose to perform hypothesis testing using mixed observations that are functions of multiple random variables. We characterize the error exponents for correctly identifying the $k$ anomalous random variables under fixed time-invariant mixed observations, random time-varying mixed observations, and deterministic time-varying mixed observations. Our error exponent characterization is through newly introduced notions of \emph{inner conditional Chernoff information} and \emph{outer conditional Chernoff information}. It is demonstrated that mixed observations can strictly improve the error exponents of hypothesis testing, over separate observations of individual random variables. We further characterize the optimal sensing vector maximizing the error exponents, which lead to explicit constructions of the optimal mixed observations in special cases of hypothesis testing for Gaussian random variables. These results show that mixed observations of random variables can reduce the number of required samples in hypothesis testing applications.


B. GOALS:
Finding $k$ anomalous random variables out of $n$ random variables


C. CONTENTS:
matlab code 
(a) HT_CLRT: Hypothesis Testing via Compressed Likelihood Ratio Test 
(b) HT_SLRT: Hypothesis Testing via Standard (Singe) Likelihood Ratio Test
(c) HT_MP: Hypothesis Testing via Message Passing method
(d) HT_LASSO: Hypothesis Testing via LASSO type method

D. HOW TO RUN (Example)
(matlab prompt)>> Comp_LRT_MP.m

E. ADDITION
(a) For a sensing matrix, we use "regular Gallager Parity Check" matrixThe code is obtained in the following link: http://www.mathworks.com/matlabcentral/fileexchange/17297-construction-of-regular-gallagher-parity-check/content/genH_regularGallagher.m?requestedDomain=www.mathworks.com


(b) Contact info.: Myung (Michael) Cho 
- Email: myung-cho@uiowa.edu, michael.myung.cho@gmail.com
- Homepage: http://user.engineering.uiowa.edu/~myucho/index.htm


F. REFERENCE PAPER: 
[1] Myung Cho, Weiyu Xu, and Lifeng Lai, "Compressed Hypothesis Testing: To Mix or Not to Mix?"




