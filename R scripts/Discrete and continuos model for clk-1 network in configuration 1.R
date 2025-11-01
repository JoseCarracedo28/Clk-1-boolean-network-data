#cLK1 Networkt configuration 1 script 
#load or install the required libraries 
library(BoolNet)
library(BoolNetPerturb)
#Clk-1 network configuration 1 
redclk1 <- tempfile(pattern = "testNet")
sink(redclk1)
cat("targets, factors\n")
cat("clk1, 0\n")
cat("ETC, clk1\n")
cat("pink1, !ETC\n")
cat("skn1, ros | !ETC\n")
cat("taf4, !clk1\n")
cat("met, !hif1\n")
cat("hlh11, !atfs1 & !mtor\n")
cat("unc51, ampk & !mtor\n")
cat("hlh30, !mtor\n")
cat("creb, taf4 & !crtc1\n")
cat("UP, !clk1 & (ETC2 | !hsp60)\n")
cat("ATG, unc51 & hlh30\n")
cat("lipl4, hlh30\n")
cat("MTG, pink1 & skn1 & unc51\n")
cat("fzo1, creb\n")
cat("clpp1, UP\n") 
cat("LDs, ATG\n")
cat("ETC2, fzo1 & MTG\n")
cat("atfs1, clpp1\n") 
cat("ros, ETC2 | (!ampk & !clk1)\n")
cat("hsp60, atfs1\n")
cat("hif1, ros\n") 
cat("atgl1, !hlh11\n")
cat("lip, lipl4 & LDs & atgl1\n")
cat("betaox, lip\n") 
cat("ATP, (ETC2 & betaox) | ETC\n")
cat("crtc1, !ampk\n")
cat("mtor, !ampk & !unc51\n")
cat("ampk, !ATP | (ros & !hif1)\n")
sink()

Nred <- loadNetwork(redclk1)
print(Nred)
####clk-1 mutant attractors
atracctoresN <- getAttractors(Nred)
plotAttractors(atracctoresN, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")
#Generation wildtype, clk-1;aak-2 double mutant, aak-2 conditions 
clk1nor <- fixGenes(Nred, "clk1", 1)
doblemut <- fixGenes(Nred, c("ampk","clk1"), c(0,0))
aak2nut <- fixGenes(Nred, c("ampk","clk1"), c(0,1))

#######wildtype conditions attractors
atracctoresclk11 <- getAttractors(clk1nor)
plotAttractors(atracctoresclk11, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")

####Double mutant clk-1;aak-2 attractors
atracctoresdoble <- getAttractors(doblemut)
plotAttractors(atracctoresdoble, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")


##################################################
#############set equations in desolve#############
##################################################
library(deSolve)
net.ode <- booleanToODE(Nred, keep.input = TRUE)
out <- ode(func = net.ode$func, 
           parms = net.ode$parameters, 
           y = net.ode$state, 
           times = seq(0, 20, 0.1))
estados_de_interes <- list(c(0,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0), 
                           c(0,1,1,1,1,0,1,1,1,1,0,0,1,0,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0), 
                           c(0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0),
                           c(0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1))

estados_de_interes
for(i in 1:4){print(estados_de_interes[[i]])}

names(estados_de_interes) <- c("Attractor1 10 states", "Attractor2 11 states", "Atractor3 16 states", "initial state")
estados_de_interes
library(ggplot2)
library(reshape2)
################################################################################
##################Continous model in clk-1 mutant conditions####################
################################################################################
for(i in 1:4){state <- validateState(estados_de_interes[[i]], Nred$genes)
net.ode <- booleanToODE(Nred, keep.input = TRUE)
out <- ode(func = net.ode$func, 
           parms = net.ode$parameters, 
           y = state, 
           times = seq(0, 50, 0.01))
outggp <- melt(as.data.frame(as.matrix(out)), id='time')
pdf(paste("/home/Estates configuration 1", names(estados_de_interes)[i] ,"pdf"))
print(ggplot(outggp, aes(time, value, col=variable)) + 
geom_line() +
ggtitle(paste(names(estados_de_interes)[i])) +
theme_bw())
dev.off()
}

################################################################################
########################Continous model in wildtype conditions##################
################################################################################
for(i in 1:4){state <- validateState(estados_de_interes[[i]], clk1nor$genes)
net.odeck <- booleanToODE(clk1nor, keep.input = TRUE)
out2 <- ode(func = net.odeck$func, 
           parms = net.odeck$parameters, 
           y = state, 
           times = seq(0, 10, 0.01))
outggp <- melt(as.data.frame(as.matrix(out2)), id='time')
pdf(paste("Documents/states_clk1_active", names(estados_de_interes)[i] ,"pdf"))
print(ggplot(outggp, aes(time, value, col=variable)) + 
        geom_line() +
        ggtitle(paste(names(estados_de_interes)[i])) +
        theme_bw())
dev.off()
}

################################################################################
##################continuos model for double mutatn conditions##################
################################################################################
for(i in 1:5){state <- validateState(estados_de_interes[[i]], doblemut$genes)
net.odedm <- booleanToODE(doblemut, keep.input = TRUE)
out3 <- ode(func = net.odedm$func, 
            parms = net.odedm$parameters, 
            y = state, 
            times = seq(0, 10, 0.01))
outggp <- melt(as.data.frame(as.matrix(out3)), id='time')
pdf(paste("/Documents/double_mutant", names(estados_de_interes)[i] ,"pdf"))
print(ggplot(outggp, aes(time, value, col=variable)) + 
        geom_line() +
        ggtitle(paste(names(estados_de_interes)[i])) +
        theme_bw())
dev.off()
}




