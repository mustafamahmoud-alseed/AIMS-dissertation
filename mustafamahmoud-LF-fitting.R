###### Mustafa Mahmoud B. G. Alseed #################################################################
#####################################################################################################
#### Building a mathematical model for lymphatic filariasis (LF) ####################################
#####################################################################################################

# Time steps annually

# libraries
library(deSolve)
library(haven)
library(boot)
library(optimx)
library(xlsx)

# Calling the dataset
setwd("/var/autofs/misc/home/mustafamahmoud/Lymphatic_Filariasis")

# Reading in the data
data <- read.xlsx(file.choose(),1)

# Viewing the data (optional)
#View(data)

#convert zonal data to time series
cases  <-ts(data[,11],st=c(2000), end=c(2018),fr=1)

# Plot the data (optional)
plot(cases,col="blue",xlab = "Year",ylab="Number of cases",tybe = "b",xlim =c(2000,2018),lwd = 2.5)
legend("topleft",legend=c("Confirmed"), col=c("blue"), lty=1 , lwd = 2.5, cex=1)#, size=14)



# Setting up function to solve the model
savModel = function(times,init,parms) { 
  # times are the time steps, init are set of initial values and parms a lis of parameters
  
  with(as.list(c(parms,init)), {
    
    # Populations
    Nh = Sh + Eh + Iha + Ihc + Rh 
    Nv = Sv + Ev + Iv
    
    amp = exp(logamp) # amplitude of the seasonality function
   # phi = exp(logphi) # phase angle of the seasonality function
    
    # definition of variables
    
    #Defining seasonal forcing function
    seas <- ( 1 + amp*cos(2*pi*times))*10*exp(-m)   # the term exp(-m) for the interventions

    
    # Forces of infection:
    
    # Infection from mosquitos to humans
    lambdah = (b*qv*Iv)/Nh
    
    
    # Infection from humans to mosquitos
    lambdav = b*qh*(Iha + Ihc)/Nh
    
    
    # Define rate equations
    
    # Human compartments
    dSh= betah + f*Rh - lambdah*Sh*seas - muh*Sh 
    dEh= lambdah*Sh*seas - (1 - p)*x*Eh - muh*Eh - p*x*Eh
    dIha= (1 - p)*x*Eh*seas - z*Iha - (muh + dh)*Iha 
    dIhc= (p*x*Eh - r*k*Ihc - (muh + dh)*Ihc + (1 - r)*k*Ihc)*seas
    dRh= r*k*Ihc*seas + z*Iha - f*Rh - muh*Rh
    
  
    # Mosquito compartments
    dSv= betav - lambdav*Sv*seas - muv*Sv 
    dEv= lambdav*Sv*seas - muv*Ev - epsilon*Ev
    dIv= epsilon*Ev - muv*Iv 
    
    
    # return the rate of change
    list(c(dSh,dEh,dIha,dIhc,dRh,dSv,dEv,dIv))
  } )
}


# Set of Parameter values
parms = c(betah = 25000.00, betav = 1000000, muh = 0.00099, 
          dh = 0.0001, muv = 0.1429, b = 250, qh = 0.01,
          qv = 0.1, x = 0.0238, z = 0.0054,
          p = 0.8, r = 0.9, k = 0.1, f = 0.238, epsilon = 0.0555, m=0.4) ## vary m from 0.04 to 0.4


# Seasonality
logamp   = log(0.6) #amplitude of seasonal function


# Demographic parameters for human population
Nh = 100000000     # population size

# Demographic parameters for mosquito population
Nv = 7000000     # population size 

# Start time amd time steps
startyear= 2000                           # starting year of simulation
tyears  <- 18.999996                      # total years of simulation (19 years)
msteps  <- 1                              # output timestep
tsteps  <- round(tyears/msteps)           # number of time steps
times   <- startyear+seq(0,tyears,msteps) # time vector


# Starting values for solver
init = c(Sh=100000000, Eh = 0, Iha=1, Ihc=1, Rh=0, Sv=7000000, Ev=0, Iv = 1)  

savModel(0,init,parms)

# Solve rate equations
run<-vode(times=times,y=init,func=savModel,parms=parms)

# Viewing the numerical solutions for the compartments (optional)
View(run) 


# Plot of symptomatic cases (for human)   (optional)
plot(times,round(run[,2]+run[,3]+run[,4]+run[,5]+run[,6]),col = "red",xlab = "Year",ylab="Number of cases")
legend("topleft",legend=c("Symptomatic"), col=c("red"), lty=1 , lwd = 2.5, cex=1)


# Susceptible humans 
plot(times,run[,2],col = "green",type="l",xlab="Year",ylab = "Population", lwd=3)
legend("topright",legend=c("S_h"), col=c("green"), lty=1 , lwd = 2.5, cex=1)

# Exposed (latent) humans 
plot(times,run[,3],col = "blue",type="l",xlab="Year",ylab = "Population", lwd=3)
legend("topleft",legend=c("E_h"), col=c("blue"), lty=1 , lwd = 2.5, cex=1)

# Infected-asymptomatic humans 
plot(times,run[,4],col = "black",type="l",xlab="Year",ylab = "Population", lwd=3)
legend("topleft",legend=c("I_ha"), col=c("black"), lty=1 , lwd = 2.5, cex=1)

# Infected-symptomatic humans 
plot(times,run[,5],col = "red",type="l",xlab="Year",ylab = "Population", lwd=3)
legend("topleft",legend=c("I_hc"), col=c("red"), lty=1 , lwd = 2.5, cex=1)

# Recoverd humans 
plot(times,run[,6],col = "purple",type="l",xlab="Year",ylab = "Population", lwd=3)
legend("topleft",legend=c("R_h"), col=c("purple"), lty=1 , lwd = 2.5, cex=1)

# Human population stratification
plot(times,run[,2],lty=1,type="l",xlab="Year", ylab="Population",lwd=3,col="green") 
lines(times,run[,3],lty=2,type="l",lwd=3,col="blue") 
lines(times,run[,4],lty=3,type="l",lwd=3,col="black") 
lines(times,run[,5],lty=4,type="l",lwd=3,col="red") 
lines(times,run[,6],lty=5,type="l",lwd=3,col="purple") 
legend("topright", legend=c("S_h", "E_h","I_ha","I_hc","R_h"), col=c("green", "blue","black","red","purple"), lty=1:5)

# Human population stratification (after 2008)
plot(times,run[,2],lty=1,type="l",xlab="Year", ylab="Population",xlim = c(2008,2018),lwd=3,col="green") 
lines(times,run[,3],lty=2,type="l",lwd=3,col="blue") 
lines(times,run[,4],lty=3,type="l",lwd=3,col="black") 
lines(times,run[,5],lty=4,type="l",lwd=3,col="red") 
lines(times,run[,6],lty=5,type="l",lwd=3,col="purple") 
legend("topright", legend=c("S_h", "E_h","I_ha","I_hc","R_h"), col=c("green", "blue","black","red","purple"), lty=1:5)


# Infectious humans (I_ha + I_hc) 
plot(times,run[,4]+run[,5],col = "red",type="l",xlab="Year",ylab = "Number of cases", lwd=3)
legend("topleft",legend=c("Modelled"), col=c("red"), lty=1 , lwd = 2.5, cex=1)


# Infectious humans (I_ha + I_hc) against the confirmed number of treated people (data)
plot(cases,col="blue",xlim =c(2000,2018),lwd = 2.5,xlab="Year",ylab="Number of cases",type = "l")
lines(times,run[,4]+run[,5],col="red",lwd = 2.5) 
legend("topleft",legend=c("Confirmed","Modelled"), col=c("blue","red"), lty=1 , lwd = 2.5, cex=1)#, size=14)


#####################################################################################################################################
###### The following part for fitting the model onto the data with interventions ####################################################
#####################################################################################################################################

# Use the seasonal force function (1 + amp*cos(2*pi*times))*10*exp(-m) in savModel 
# Repreat all the previous steps 
# Repeat the following process for run[,4], run[,5] and run[,4]+run[,5] 

# Minimise the sum of square errors
sir.ss<-function(data1,parms){
  model<- ode(times=times,y=init, func=savModel, parms=parms)
  error1<-model[,4]-cases 
 # error2<-model[,5]-cases
  sse<-0.7*sum(error1^2)   # Use 0.7*sum(error1^2 + error2^2) for run[4]+run[,5]
  return(sse)
}

sir.ss(data1, parms)

# Fitting process
fit0<-optim(par=parms,fn=sir.ss,data=data1)
run_fit0<- ode(times=times, y=init, func=savModel,parms=fit0$par)

cases_fit    <-ts(run_fit0[,4],st=c(2000), end=c(2018),fr=1)

# Plot fitting results 
plot(times,cases,col="blue",lwd=2,xlab="Year",ylab="Number of cases",type="b")
lines(cases_fit,col="red",lwd=2)
legend("topleft", legend=c("Modelled","Confirmed"), col=c("red","blue"), lty=1:2 , lwd = 2, cex=0.8)  

# Use the optimization method (L-BFGS-B) in fitting
fit0<-optim(par=parms,fn = sir.ss, gr = NULL, method = c("L-BFGS-B"),control = list(), hessian = FALSE,data=data1)
run_fit0<- ode(times=times, y=init, func=savModel,parms=fit0$par)

cases_fit    <-ts(run_fit0[,4],st=c(2000), end=c(2018),fr=1)

# Plot fitting results 
plot(times,cases,col="blue",lwd=2,xlab="Year",ylab="Number of cases",type="b")
lines(cases_fit,col="red",lwd=2)
legend("topleft", legend=c("Modelled","Confirmed"), col=c("red","blue"), lty=1:2 , lwd = 2, cex=0.8) 

