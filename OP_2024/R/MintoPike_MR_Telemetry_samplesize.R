# dsftools detection probability
# dsftools stratified multinomial??
# recapr n2rr

library(dsftools)   # for detection probability
library(tidyverse)  # for data operations and plotting
library(patchwork)  # for plotting
library(recapr)     # verifying Robson Regier stuff


write_output <- TRUE  ### to save output to files


## --------------------- detection probability ------------------------ ##

# detection probability: print to console
detection_probability(n_raw=30,
                      prop_usedby = c(0.2, 0.1, 0.05),
                      observe_at_least = 1:3,
                      prop_ofareas = c(0.9, 1))

# detection probability: make a table to write to .csv
table1 <- detection_probability(n_raw=30,
                      prop_usedby = c(0.2, 0.1, 0.05),
                      observe_at_least = 1:3,
                      prop_ofareas = 0.9) %>%
  rename(`Proportion of Population` = prop_usedby) %>%
  mutate(observe_at_least = paste(observe_at_least, "fish")) %>%
  rename(`Detection Threshold` = observe_at_least) %>%
  rename(`Probability of Detecting a Single Area` = p_singlearea) %>%
  rename(`Probability of Detecting 90% of Areas` = p_multipleareas) %>%
  mutate(`Probability of Detecting 100% of Areas` =
           detection_probability(n_raw=30,
                                 prop_usedby = c(0.2, 0.1, 0.05),
                                 observe_at_least = 1:3,
                                 prop_ofareas = 0.9)$p_multipleareas)
table1
if(write_output) {
  write.csv(table1, file="OP_2024/R_output/detection_prob.csv")
}



# detection probability: make a figure!
ggprob <- detection_probability(n_raw=30,
                      prop_usedby = seq(0.01, 0.25, by=0.001),
                      observe_at_least = 1:3,
                      prop_ofareas = c(0.8,0.9, 1)) %>%
  mutate(Threshold = paste(observe_at_least,"fish")) %>%
  mutate(prop_ofareas = factor(paste(100*prop_ofareas, "% of Areas"),
                               levels = paste(seq(80, 100, by=10),"% of Areas")))

plot1 <- ggprob %>%
  ggplot(aes(x = prop_usedby, y=p_singlearea,
             col = Threshold)) +
  geom_line() +
  theme_bw() +
  labs(x="Proportion of Population",
       y="Detection Probability",
       title="Probability of Detecting a Single Area") +
  theme(text = element_text(family = "serif"))

plot2 <- ggprob %>%
  ggplot(aes(x = prop_usedby, y=p_multipleareas,
             col = Threshold)) +
  facet_wrap(~prop_ofareas) +
  geom_line() +
  theme_bw() +
  labs(x="Proportion of Population",
       y="Detection Probability",
       title="Probability of Detecting Multiple Areas") +
  theme(text = element_text(family = "serif")) +
  theme(axis.text.x=element_text(angle=60, hjust=1))

plot1 / plot2

if(write_output) {
  ggsave(plot1 / plot2,
         filename="OP_2024/R_output/detection_prob.png",
         height=6, width=7, units="in")
}





## ------------------- sample size Robson Regier --------------------- ##
n2RR(N=11992, n1=700)
n2RR(N=11992, n1=700, conf=0.9, acc=0.25)

# reproducing Table 4
n1_values <- c(seq(700, 1200, by=100), 1237, 1300, 1400)
N_values <- c(11992, 14675, 17358)

# initializing table
tbl4 <- matrix(nrow=10, ncol=3)
rownames(tbl4) <- c("n1=n2",n1_values)
colnames(tbl4) <- N_values

# filling in top row
tbl4[1,] <- sapply(N_values, \(x) n2RR(N=x, n1=700, conf=0.9, acc=0.25)[[1]][1,3])

# filling in the rest
for(j in 1:3) {
  tbl4[2:10, j] <- sapply(n1_values, \(x) n2RR(N=N_values[j], n1=x, conf=0.9, acc=0.25)[[1]][1,2])
}
tbl4

# make a plot of this??
tbl4df <- as.data.frame(tbl4[-1,])
tbl4df$n1 <- n1_values
nequal <- unlist(unname(tbl4[1,]))

# ggplot version - abandoned
tbl4df %>%
  pivot_longer(names_to = "N", values_to = "n2", cols=1:3)
tbl4df%>%
  ggplot(aes(x = n1, y= n2, col=N)) +
  theme_bw()+
  geom_line() +
  theme(text = element_text(family = "serif"))

tbl4df

# base plot version - used
N_legend <- paste(floor(N_values/1000), N_values %% 1000, sep=",")

if(write_output) {
  png(filename = "OP_2024/R_output/sampsize_RR.png", width=5, height=5, units="in", res=300)
}
par(family="serif")
plot(NA, xlim=range(696,n1_values), ylim=range(tbl4), xlab="n1", ylab="n2")
for(j in 1:3) lines(tbl4df[,4],tbl4df[,j], col=j+1)
points(x=nequal, y=nequal, col=1+(1:3), pch=16)
legend("topright", col=1+(1:3), lwd=1, pch=16, legend=N_legend, title="Assumed N")
if(write_output) {
  dev.off()
}
