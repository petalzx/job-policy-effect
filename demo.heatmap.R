rm(list=ls())
squared.regret.rule<- matrix(c(0.77,0.77,0.94, 0.90,0.94,0.94,0.80,
                               0.42,0.59,0.79,0.75,0.92,0.18,0.89,
                               0.83 ,0.47,0.25,0.00
                               ),3, 6)

library(reshape)                                                # Load reshape package

library(ggplot2)   

our.rule.plot <- melt(squared.regret.rule)
colnames(our.rule.plot) <-c("pre-program income","pre-program education","decision rule")

ggp.1 <- ggplot(our.rule.plot, aes(`pre-program education`,`pre-program income`)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = `decision rule`))
ggp.1 + scale_fill_gradient(low = "white", high = "#d12168",lim = c(0,1))+
 labs(x = "pre-program education bracket", 
       y = "pre-program income bracket", 
       title = "Our proposed approach")+
annotate("text", x=6, y=1, label= "0.47")+
annotate("text", x=3, y=2, label= "0.42")+
annotate("text", x=3, y=3, label= "0.59")  

