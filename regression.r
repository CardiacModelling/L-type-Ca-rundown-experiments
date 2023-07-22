# """
# Perform multi-linear regression for all the run down rates with temperature, 
# thold, and INaCa as the indepedent variables 
# """

# get path
pathwork = getwd() 
pathrel = "output/r_rate_database.csv"
path = paste(pathwork, pathrel, sep = "/") 

# load data
dataset = read.csv(path)
dataset['thold'] = log(dataset['thold'])

# Reorder the levels to set 310 as the reference level
dataset$Temperature = factor(dataset$Temperature, levels = c(310, 298))

# perform multi-linear regression
multi.fit = lm(Run.rate~thold + factor(Temperature) + factor(INaCa), data=dataset)
res = summary(multi.fit)

# print output
print(summary(multi.fit))
print(confint(multi.fit))

# For H1: beta < 0
p_value_less = pt(coef(res)[, 3], multi.fit$df, lower = TRUE)
print(p_value_less)

# For H1: beta > 0
p_value_greater = pt(coef(res)[, 3], multi.fit$df, lower = FALSE)
print(p_value_greater)
#%%%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%

# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.158729 -0.016681  0.003505  0.026190  0.121515 

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.158795   0.019715   8.055 2.16e-13 ***
# thold                  -0.016373   0.006545  -2.502  0.01342 *  
# factor(Temperature)298 -0.062213   0.007123  -8.734 4.16e-15 ***
# factor(INaCa)On        -0.021030   0.007023  -2.994  0.00321 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.04355 on 152 degrees of freedom
# Multiple R-squared:  0.4045,    Adjusted R-squared:  0.3928 
# F-statistic: 34.42 on 3 and 152 DF,  p-value: < 2.2e-16

#                              2.5 %       97.5 %
# (Intercept)             0.11984411  0.197745292
# thold                  -0.02930387 -0.003442791
# factor(Temperature)298 -0.07628548 -0.048139628
# factor(INaCa)On        -0.03490557 -0.007153578



# Null hypothesis: Coefficient is 0
# Std Error: Lower error means better estimation
# t value: Larger t-value (]-2, 2[) indicate evidence against the null hypothesis 
# p-value: a p-value less than .05 rejects the null hypothesis
# Confidence intervals: a confidence interval that does not include 0 is good

# p-values for H1: beta < 0
# (Intercept)            :  1.000000e+00
# thold                  :  6.708997e-03
# factor(Temperature)298 :  2.079688e-15
# factor(INaCa)On        :  1.606375e-03

# p-values for H1: beta > 0
# (Intercept)            :  1.078998e-13
# thold                  :  9.932910e-01 
# factor(Temperature)298 :  1.000000e+00
# factor(INaCa)On        :  9.983936e-01
