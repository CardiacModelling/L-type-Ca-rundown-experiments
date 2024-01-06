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
print("############### Summary of model fit ###############")
print(summary(multi.fit))
print(confint(multi.fit))
print("####################################################")

# For H1: beta < 0
p_value_less = pt(coef(res)[, 3], multi.fit$df, lower = TRUE)
print("## Test for directionality (beta coefficients <0) ##")
print(p_value_less)
print("####################################################")

# Statistical reduction by addition of temperature
initial.model <- lm(Run.rate ~ 1, data = dataset)  
full.model <- lm(Run.rate ~ factor(Temperature), data = dataset)

# Perform F-test
f_test_result <- anova(initial.model, full.model)

print("##### Statistical reduction due to temperature #####")
print(f_test_result)
print("####################################################")

# Statistical reduction by addition of INaCa
initial.model <- lm(Run.rate ~ factor(Temperature), data = dataset)  
full.model <- lm(Run.rate ~ factor(Temperature) + factor(INaCa), data = dataset)

# Perform F-test
f_test_result <- anova(initial.model, full.model)

print("######## Statistical reduction due to INaCa ########")
print(f_test_result)
print("####################################################")

# Statistical reduction by addition of thold
initial.model <- lm(Run.rate ~ factor(Temperature)+ factor(INaCa), data = dataset)  
full.model <- lm(Run.rate ~ factor(Temperature) + factor(INaCa) + thold, data = dataset)

# Perform F-test
f_test_result <- anova(initial.model, full.model)

print("######## Statistical reduction due to thold ########")
print(f_test_result)
print("####################################################")



# For H1: beta > 0
# p_value_greater = pt(coef(res)[, 3], multi.fit$df, lower = FALSE)
# print(p_value_greater)

#%%%%%%%%%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# [1] "############### Summary of model fit ###############"

# Call:
# lm(formula = Run.rate ~ thold + factor(Temperature) + factor(INaCa), 
#     data = dataset)

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
# [1] "####################################################"
# [1] "## Test for directionality (beta coefficients <0) ##"
#            (Intercept)                  thold factor(Temperature)298        factor(INaCa)On 
#           1.000000e+00           6.708997e-03           2.079688e-15           1.606375e-03 
# [1] "####################################################"
# [1] "##### Statistical reduction due to temperature #####"
# Analysis of Variance Table

# Model 1: Run.rate ~ 1
# Model 2: Run.rate ~ factor(Temperature)
#   Res.Df     RSS Df Sum of Sq      F    Pr(>F)    
# 1    155 0.48405                                  
# 2    154 0.31631  1   0.16773 81.661 6.389e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# [1] "####################################################"
# [1] "######## Statistical reduction due to INaCa ########"
# Analysis of Variance Table

# Model 1: Run.rate ~ factor(Temperature)
# Model 2: Run.rate ~ factor(Temperature) + factor(INaCa)
#   Res.Df     RSS Df Sum of Sq      F   Pr(>F)   
# 1    154 0.31631                                
# 2    153 0.30010  1  0.016211 8.2647 0.004619 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# [1] "####################################################"
# [1] "######## Statistical reduction due to thold ########"
# Analysis of Variance Table

# Model 1: Run.rate ~ factor(Temperature) + factor(INaCa)
# Model 2: Run.rate ~ factor(Temperature) + factor(INaCa) + thold
#   Res.Df     RSS Df Sum of Sq      F  Pr(>F)  
# 1    153 0.30010                              
# 2    152 0.28824  1  0.011868 6.2586 0.01342 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# [1] "####################################################"
