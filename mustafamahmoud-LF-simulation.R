###### Mustafa Mahmoud B. G. Alseed ###############################################################
###################################################################################################
#### LF-simulation of SEIR Model for humans and mosquito vectors ####################################
###################################################################################################

# The script is for simulating the SEIR model for lymphatic filarisis (with and without interventions), and plotting the infected states for the human population with time, considering the forces of infection, initial values and parameter  
# For the numerical solution we used the library "deSolve"

# Import the library "deSolve"
library(deSolve)

# Equations of the simulated model 
seir_LF <- function(times, start, parms)  {
  with(as.list(c(parms, start)), {   
    #X represents the set of initial values,t represents times and parms is a list of parameters
    amp = exp(logamp) # amplitude of the seasonality function
    #phi = exp(logphi) # phase angle of the seasonality function
    
    # For humans
    season <- ( 1 + amp*cos(2*pi*times/180))#*exp(-m)
    Nh=(Sh+Eh+Iha+Ihc+Rh)
    lambdah = (b*qv*Iv)/Nh  # the force of infection for humans 
    
    # For mosquitoes 
    Nv=(Sv+Ev+Iv)
    lambdav = b*qh*(Iha + Ihc)/Nh   # the force of infection for vectors
    
    #sub-System for human population
    dSh= betah + f*Rh - lambdah*Sh*season - muh*Sh 
    dEh= lambdah*Sh*season - (1 - p)*x*Eh - muh*Eh - p*x*Eh
    dIha= (1 - p)*x*Eh*season - z*Iha - (muh + dh)*Iha 
    dIhc= (p*x*Eh - r*k*Ihc - (muh + dh)*Ihc + (1 - r)*k*Ihc)*season
    dRh= r*k*Ihc*season + z*Iha - f*Rh - muh*Rh
    
    #sub-System for mosquito population
    dSv= betav - lambdav*Sv*season - muv*Sv 
    dEv= lambdav*Sv*season - muv*Ev - epsilon*Ev
    dIv= epsilon*Ev - muv*Iv 
    
    output <- c(dSh,dEh,dIha,dIhc,dRh,dSv,dEv,dIv)
    
    list(output)
  })
}

# The initial values
start<-c(Sh=100000, Eh = 0 ,Iha=1,Ihc=1,Rh=0,Sv=7000, Ev=0, Iv = 1)

# The parameters 
parms <- c(betah = 25.00, betav = 1000, muh = 0.00099, 
           dh = 0.0001, muv = 0.1429, b = 250, qh = 0.01,
           qv = 0.1, x = 0.0238, z = 0.0054, 
           p = 0.8, r = 0.9, k = 0.1, f = 0.238, epsilon = 0.0555, m=0.4) ## vary m from 0.04 to 0.4

logamp   = log(0.6) #amplitude of seasonal function
#logphi   = log(2.9) #phase angle of seasonal fuction

# vector of timesteps
times <- seq(0, 730, 1/30)  

run<-ode(times=times, y=start, func=seir_LF,parms=parms)    
#start=at what stage do you want the model to start?

# Viewing results of numerical solution for all the compartments (optional)
#View(run)


#par(mfrow=c(2,2))

# Plot for all humans population (optional)
plot(times,run[,2], col="green",lwd = 3, type="l",xlab = "Time (days)",ylab = "Population")
lines(times,run[,3], col="blue",lwd = 3)
lines(times,run[,4], col="black",lwd = 3)
lines(times,run[,5], col="red",lwd = 3)
lines(times,run[,6], col="purple",lwd = 3)
legend("topright",legend=c("S_h","E_h","I_ha","I_hc","R_h"),lwd = 3,cex = 1 ,col=c("green","blue", "black","red","purple"), lty=c(1,1,1,1,1))


# Plot for humans population (for infected states)
plot(times,run[,3], col="blue",lwd = 3, type="l",xlab = "Time (days)",ylab = "Population")
lines(times,run[,4], col="black",lwd = 3)
lines(times,run[,5], col="red",lwd = 3)
lines(times,run[,6], col="purple",lwd = 3)
legend("topright",legend=c("E_h","I_ha","I_hc","R_h"),lwd = 3,cex = 1 ,col=c("blue", "black","red","purple"), lty=c(1,1,1,1))

# Plot for mosquitoes (vectors) population (optional)
plot(times,run[,7], col="green",lwd = 3, type="l",xlab = "Time (days)",ylab = "Number of Mosquito", main="Mosquito Distributions")
lines(times,run[,8], col="blue",lwd = 3)
lines(times,run[,9], col="red",lwd = 3)
legend("topright",legend=c("S_v","E_v","I_v"),lwd = 3,cex = 1, col=c("green","blue", "red"), lty=c(1,1,1))
