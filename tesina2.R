rm(list=ls())
library(readr)
library(akima)
library(dplyr)
library(geoR)
library(ggplot2)
library(maps)
library(viridis)
library(spatstat)
library(spdep)
library(scatterplot3d)
library(sgeostat)

setwd("C:/Users/Patrizio/Desktop/TUTTO/Ud'A/CLEBA/statistical learning/HOMEWORK SPATIAL ANALYSIS")
Pollution <- read_csv("Pollution.csv", col_types = cols(mean.obs.value = col_number()))
head(Pollution)
############################################ punto a #############################

dat <- filter(Pollution, Year == 2014, pollutant == "O3")
Y <- 100*dat[,6] 
Y <- as.numeric(unlist(Y))
s <- as.matrix(dat[,3:4])

sum(is.na(Y) == T)
sum(is.na(dat) == T)
sum(is.na(s) == T)

############################################ punto b #############################

# par(mfrow=c(1,1),pty='s',cex.lab=0.7,cex.main=1.0)
hist(Y,freq=FALSE, breaks=20)
lines(density(Y))

############################################ punto c #############################

data <- dat %>% 
   select(c(3,4,6)) %>%
   as.data.frame()
sum(is.na(data))
plot(as.geodata(data),scatter3d=TRUE,highlight=TRUE) # 4 grafici


#variogdata_cloud_smooth<-variog(as.geodata(data),
                         #coords = s,
                         #data = data$mean.obs.value, 
                         #option="smooth") # option smooth o cloud?
#plot(variogdata_cloud_smooth,main='grafico del variogramma',pch='.',cex=0.8)

variogdata_cloud_cloud <- variog(as.geodata(data),
                              coords = s,
                                 data = data$mean.obs.value, 
                                 option="cloud")
plot(variogdata_cloud_cloud,main='nuvola del variogramma',pch='.',cex=0.8)

pollutant.bin<-variog(as.geodata(data),option="bin")
points( pollutant.bin$u , pollutant.bin$v ,pch=19 , col=2) # calcola la media in ogni classe
abline(v= pollutant.bin$bins.lim)
# rappr del variogramma con dati suddivisi per classi
# ogni classe è definita dalla linea abline
# il grafico mostra la continuità spaziale del processo. 
# per valori più bassi li troviamo a distanze più piccole
# ad aumentare della distanza aumenta la diversità tra le coppie dei dati

############################################ punto d #############################

x1 <- s[,1]
x2 <- s[,2]
aq.fit.1 <- lm(Y~x1+x2)
aq.fit.2 <- lm(Y~x1+x2+I(x1^2)+I(x2^2)+x1*x2)
aq.fit.3 <- lm(Y~x1+x2+I(x1^2)+I(x2^2)+x1*x2+I(x1^3)+I(x2^3)+I(x1^2)*x2+I(x2^2)*x1)


# z~x+y+I(x^2)+I(y^2)+x*y+I(x^3)+I(y^3)+I(x^2)*y+I(y^2)*x

AIC(aq.fit.1, aq.fit.2, aq.fit.3)

BIC(aq.fit.1 ,aq.fit.2, aq.fit.3)

##### griglia ######

x1grid<-seq(min(x1),max(x1),length=20)
x2grid<-seq(min(x2),max(x2),length=20)
aq.grid<- expand.grid(x1=x1grid,x2=x2grid) # qui faccio l'incrocio tra tutte le possibili coppie
# di coordinate
aq.grid$yr <- range(aq.grid$x2)
aq.grid$ys <- aq.grid$yr[2] - aq.grid$yr[1]
aq.grid$xr <- range(aq.grid$x1)
aq.grid$xs <- aq.grid$xr[2] - aq.grid$xr[1]
aq.grid$xy <- data.frame(cbind
                         (
                            c
                            (matrix(aq.grid$x1, 
                                    length(aq.grid$x1), length(aq.grid$x2))),
                            c
                            (matrix(aq.grid$x2, 
                                     length(aq.grid$x1), length(aq.grid$x2), 
                                     byrow=T))
                            )
                         )
colnames(aq.grid$xy) <- c("x1", "x2")
locations<-aq.grid$xy

####### primo grado ########

aq.surf<-predict(aq.fit.1,newdata=aq.grid)
---------------------------------------------------------------
graphics.off() # per resettare grafici
par(mar=c(1, 1, 1, 1)) # reset
-----------------------------------------------------------------
persp(x1grid,x2grid,matrix(aq.surf,20,20),xlab="x1",ylab="x2",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Trend lineare",box=TRUE)

dist <- seq(0,10,length=30)
res=aq.fit.1$res
vg <- variog(data = res,coord = s, uvec = dist)
plot(vg,main='Variogramma campionario su res trend',pch='*')

##### secondo grado #####
aq.surf2<-predict(aq.fit.2,newdata=aq.grid)
persp(x1grid,x2grid,matrix(aq.surf2,20,20),xlab="x1",ylab="x2",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Trend quadrato",box=TRUE)
res2=aq.fit.2$res
vg2 <- variog(data = res2,coord = s[,1:2], uvec = dist)
par(mfrow=c(1,1))
plot(vg2,main='Variogramma campionario su res trend',pch='*')

############# terzo grado ####################
aq.surf3<-predict(aq.fit.3,newdata=aq.grid)
persp(x1grid,x2grid,matrix(aq.surf3,20,20),xlab="x1",ylab="x2",zlab=NULL,theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Trend cubico",box=TRUE)
res3=aq.fit.3$res
Wres <- data.frame(coord=s[,1:2], res_reg = res3)
WW = as.geodata(Wres)
vg3 <- variog(WW, uvec = dist) 

---------------------------------------------------------------
graphics.off() # per resettare grafici
par(mfrow=c(1,1)) # reset
-----------------------------------------------------------------
   
plot(vg3,main='Variogramma campionario su res trend',pch='*', )


############# altro grafico #############################
# prendo la quadratica
length(res2)
length(s)
res.interp<-interp(s[,1],s[,2],aq.fit.2$res) 
persp(res.interp$x,res.interp$y,res.interp$z,xlab="x",ylab="y",zlab='residuals',theta=60,
      phi=30,expand=0.9, col = "green",ltheta = 120, shade = 0.75, ticktype = "simple",
      main="Residui",box=TRUE)



############################################ punto e #############################

#parametri fissati ad occhio
tau2=0.075 # nugget
sig2=0.15 # partial sill
phi=3.5 # range
####DA RIVEDERE
vg_fit <- (sig2 + tau2) - sig2*exp(-dist/phi)
lines(vg_fit )
      
############################################ punto f #############################

X0 <- seq(-135,-65,0.1)
Y0 <- seq(25,50,0.1)
s0 <- as.matrix(expand.grid(X0,Y0))

############################################ punto g #############################
---------------------------------------------------------------
   graphics.off() # per resettare grafici
par(mar=c(1,1,1, 1)) # reset
-----------------------------------------------------------------
map("usa")
points(s0,pch=19,cex=0.1)
map("usa",add=TRUE,col="white")
map("state","north carolina")
points(s0,pch=19,cex=0.1)
inusa <- map.where("usa",s0[,1],s0[,2])
s0 <- s0[!is.na(inusa),]
map("usa")
points(s0)
map("state","north carolina")
points(s0,pch=19,cex=0.1)


############################################ punto h #############################
c=as.geodata(data)
mle.wolf<-likfit(c,ini.cov.pars=c(0.1,2), cov.model="exponential",trend= "2nd" ,nugget = 0.1, fix.nugget=FALSE)
# il cubico dà problemi

mle.wolf
# il sigmasq è la varianza, è ottenuto facendo derivate parziali rispetto a ciascuna variabile
# il tausquare è il valore al quale la retta immaginaria interseca le ordinate

# (è una funzione multivariata, non ha una variabile sola ma 6 var)
# il practical range è il valore oltre il quale la correlazione si abbatte
# 5345 è il likelihood, si sceglie quella più grande in un confronto
# ma bisogna sempre usare il BIC per pesare anche la complessità del modello

############################################ punto i #############################



krige.par<-krige.control(type.krige='ok',cov.pars=mle.wolf$cov.pars,
                         cov.model="exponential",
                         trend.d= "2nd" ,
                         trend.l= "2nd")

kriging<-krige.conv(c,locations= s0, krige=krige.par)  
# con s0 gli indico la location dove andare a prevedere 

############################################ punto l #############################
# plot osservato
df <- data.frame(long=s[,1],lat=s[,2],Y=Y)
ggplot(df, aes(long, lat)) +
   borders("state") +
   geom_point(aes(colour = Y)) +
   scale_colour_gradientn(colours = viridis(10)) +
   xlab("")+ylab("")+labs(title="Ozone (ppm), 2014")+
   coord_fixed()
# Plot Kriging predictions
df1 <- data.frame(long=s0[,1],lat=s0[,2],Y=kriging$pred)
ggplot(df1, aes(long, lat)) +
   borders("state") +
   geom_raster(aes(fill = Y)) +
   scale_fill_gradientn(colours = viridis(10))+
   xlab("")+ylab("")+labs(title="Predicted ozone (ppm), 2014")+
   coord_fixed()

############################################ punto m #############################

# Plot dell'errore quadratico medio
se<-sqrt(kriging$krige.var)

plot(aq.grid,
     type="n",xlab="x",ylab="y", 
     xlim=c(aq.grid$xr[1], aq.grid$xr[2]),
     ylim=c(aq.grid$yr[1], aq.grid$yr[2]))	

image(aq.grid$x1,
      aq.grid$x2,
      matrix(kriging$pred,
             length(aq.grid$x1),
             length(aq.grid$x2)), 
      add=T)

contour(aq.grid$x1,
        aq.grid$x2,
        matrix(se,length(aq.grid$x1),length(aq.grid$x2)),
        add=T,
        nlevels=5)

points(s[,1],s[,2],pch=20)
title("Errore quadratico medio/n di previsione")


# 
