###Script to run both configurations attractros in an asyncronous update model
library(BoolNet)
####Network configuration 1
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

toSBML(Nred, file = "clk1 network configuration 1.sbml", generateDNFs = FALSE, saveFixed = TRUE)
####clk-1 mutant attractors
atracctoresN <- getAttractors(Nred, type = "asynchronous")
plotAttractors(atracctoresN, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")

#Generation wildtype and clk-1;aak-2 couble mutant conditions 
clk1nor <- fixGenes(Nred, "clk1", 1)
doblemut <- fixGenes(Nred, c("ampk","clk1"), c(0,0))



#######wildtype conditions attractors
atracctoresclk11 <- getAttractors(clk1nor, type = "asynchronous")
plotAttractors(atracctoresclk11, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")

####Double mutant clk-1;aak-2 attractors
atracctoresdoble <- getAttractors(doblemut, type = "asynchronous")
plotAttractors(atracctoresdoble, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")



###########Network configuration 2###########
redclk12 <- tempfile(pattern = "testNet")
sink(redclk12)
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
cat("ros, !ETC & (!ampk | hif1)\n")
cat("hsp60, atfs1\n")
cat("hif1, ros & !ampk\n") 
cat("atgl1, !hlh11\n")
cat("lip, lipl4 & LDs & atgl1\n")
cat("betaox, lip\n") 
cat("ATP, (ETC2 & betaox) | ETC\n")
cat("crtc1, !ampk\n")
cat("mtor, !ampk & !unc51\n")
cat("ampk,  (ros & !hif1)\n") 
sink()

Nred2 <- loadNetwork(redclk12)
print(Nred2)

####clk-1 mutant attractors
atracctoresN2 <- getAttractors(Nred2, type = "asynchronous")
plotAttractors(atracctoresN2, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")

#Generation wildtype and clk-1;aak-2 couble mutant conditions 
clk1nor2 <- fixGenes(Nred2, "clk1", 1)
doblemut2 <- fixGenes(Nred2, c("ampk","clk1"), c(0,0))

#######wildtype conditions attractors
atracctoresclk112 <- getAttractors(clk1nor2, type = "asynchronous")
plotAttractors(atracctoresclk112, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")

####Double mutant clk-1;aak-2 attractors
atracctoresdoble2 <- getAttractors(doblemut2, type = "asynchronous")
plotAttractors(atracctoresdoble2, onColor ="#20A387FF", offColor = "#39568CFF", borderColor = "black")


