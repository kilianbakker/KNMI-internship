testIndices <- c(1:10)
X<-scale(usair[-testIndices,2:length(colnames(usair))])
form <- paste0("y","~",paste0(colnames(usair)[2:length(colnames(usair))], collapse = "+"))
form2 <- paste0("y","~",paste0(paste0("pb(",colnames(usair)[2:length(colnames(usair))],")"), collapse = "+"))
form3 <- paste0("y","~",paste0(paste0("ridge(",colnames(usair)[2:length(colnames(usair))],")"), collapse = "+"))

fitcontrol <- gamlss.control(c.crit = 0.02, n.cyc = 200, mu.step = 1, sigma.step = 1, nu.step = 1, tau.step = 1)
fitcontrol2 <- glim.control(cc = 0.02, cyc = 200, glm.trace = F, bf.cyc = 200, bf.tol = 0.02, bf.trace = F)

m00 <- gamlss(as.formula(form), data = usair[-testIndices,], family = NO(mu.link = "identity", sigma.link = "identity"))
m0 <- gamlss(y~X, data = usair[-testIndices,], family = NO(mu.link = "identity", sigma.link = "identity"))
m0 <- gamlss(as.formula(form2), data = usair[-testIndices,], family = NO(mu.link = "identity", sigma.link = "identity"))
xupdate <- ri(x1, Lp = 1, method = "GAIC", start = 0.05, Lp = 1, kappa = 1e-05, iter = 10000, c.crit = 1e-03, k = 2)
m0 <- gamlss(as.formula(form3), data = usair[-testIndices,], family = NO(mu.link = "identity", sigma.link = "identity"), control = fitcontrol, i.control = fitcontrol2, trace = F)

m0 <- gamlss(y~ri(X, lambda = 1, Lp = 1), data = usair[-testIndices,], family = NO(mu.link = "identity", sigma.link = "identity"), control = fitcontrol, i.control = fitcontrol2, trace = F)
m0$mu.coefficients <- c(m0$mu.coefficients[1], m0$mu.coefSmo[[1]]$coef)
                        #,m0$mu.coefSmo[[2]]$coef,m0$mu.coefSmo[[3]]$coef,m0$mu.coefSmo[[4]]$coef,m0$mu.coefSmo[[5]]$coef,m0$mu.coefSmo[[6]]$coef)
names(m0$mu.coefficients) <- c("(Intercept)", "x1","x2","x3","x4","x5","x6")
predictions2 <- rep(m0$mu.coefficients[1], length(testIndices)) + as.matrix(usair[testIndices,-1])%*%m0$mu.coefficients[-1]
print(RMSE(usair[testIndices,1], predictions2))
predictions1 <- rep(m00$mu.coefficients[1], length(testIndices)) + as.matrix(usair[testIndices,-1])%*%m00$mu.coefficients[-1]

predict(m0, newdata = usair[testIndices,2:length(colnames(usair))], what = c("mu"))

plot(getSmo(m0))

