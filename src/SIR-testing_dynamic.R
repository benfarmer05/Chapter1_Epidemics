

rm(list=ls())

library(deSolve)
library(tidyverse)
library(ggpmisc)
library(gridExtra)
library(dplyr)
options(scipen = 999)

SIR = function(t,y,p){
  {
    S.LS = y[1]
    I.LS = y[2]
    R.LS = y[3]
    
    S.MS = y[4]
    I.MS = y[5]
    R.MS = y[6]
    
    S.HS = y[7]
    I.HS = y[8]
    R.HS = y[9]
    
    P = y[10]
  }
  with(as.list(p),{
    
    
    
    # # # version with static transmission modifier and a capped disease pool
    # transmission_modifier = (1 / (1 + exp(-l * (C))) + offset)
    # 
    # # Limit the product of susceptibles and disease pool
    # S.LS_P <- ifelse(S.LS * P > N.LS, N.LS, S.LS)
    # S.MS_P <- ifelse(S.MS * P > N.MS, N.MS / P, S.MS)
    # S.HS_P <- ifelse(S.HS * P > N.HS, N.HS / P, S.HS)
    # 
    # # Limit the product of susceptibles and disease pool
    # pool = 1
    # N.group = 2
    # ifelse(pool > N.group, N.group, pool)
    
    # # version with static transmission modifier, and no logic affecting disease pool
    # transmission_modifier = (1 / (1 + exp(-l * (C))) + offset)
    # 
    # dS.LS.dt = -b.LS*S.LS*(P) / (N.LS) * transmission_modifier
    # dI.LS.dt = b.LS*S.LS*(P) / (N.LS) * transmission_modifier - g.LS*I.LS
    # dR.LS.dt = g.LS*I.LS
    # 
    # dS.MS.dt = -b.MS*S.MS*(P) / (N.MS) * transmission_modifier
    # dI.MS.dt = b.MS*S.MS*(P) * transmission_modifier / (N.MS) - g.MS*I.MS
    # dR.MS.dt = g.MS*I.MS
    # 
    # dS.HS.dt = -b.HS*S.HS*(P) / (N.HS) * transmission_modifier
    # dI.HS.dt = b.HS*S.HS*(P) * transmission_modifier / (N.HS) - g.HS*I.HS
    # dR.HS.dt = g.HS*I.HS

    # # version with transmission modifier dependent on cover.group
    # transmission_modifier.LS = (1 / (1 + exp(-l * (C.LS))) + offset)
    # transmission_modifier.MS = (1 / (1 + exp(-l * (C.MS))) + offset)
    # transmission_modifier.HS = (1 / (1 + exp(-l * (C.HS))) + offset)
    # 
    # dS.LS.dt = -b.LS*S.LS*(P) / (N.LS) * transmission_modifier.LS
    # dI.LS.dt = b.LS*S.LS*(P) / (N.LS) * transmission_modifier.LS - g.LS*I.LS
    # dR.LS.dt = g.LS*I.LS
    # 
    # dS.MS.dt = -b.MS*S.MS*(P) / (N.MS) * transmission_modifier.MS
    # dI.MS.dt = b.MS*S.MS*(P) / (N.MS) * transmission_modifier.MS - g.MS*I.MS
    # dR.MS.dt = g.MS*I.MS
    # 
    # dS.HS.dt = -b.HS*S.HS*(P) / (N.HS) * transmission_modifier.HS
    # dI.HS.dt = b.HS*S.HS*(P) / (N.HS) * transmission_modifier.HS - g.HS*I.HS
    # dR.HS.dt = g.HS*I.HS
    
    # # version with transmission modifier dependent on cover.group, but coded better
    # prop.group = c(N.LS / (N.LS + N.MS + N.HS), N.MS / (N.LS + N.MS + N.HS), N.HS / (N.LS + N.MS + N.HS))
    # transmission_modifier = (1 / (1 + exp(-l * (cover.tester/prop.group))) + offset)
    # 
    # dS.LS.dt = -b.LS*S.LS*(P) / (N.LS) * transmission_modifier[1]
    # dI.LS.dt = b.LS*S.LS*(P) / (N.LS) * transmission_modifier[1] - g.LS*I.LS
    # dR.LS.dt = g.LS*I.LS
    # 
    # dS.MS.dt = -b.MS*S.MS*(P) / (N.MS) * transmission_modifier[2]
    # dI.MS.dt = b.MS*S.MS*(P) / (N.MS) * transmission_modifier[2] - g.MS*I.MS
    # dR.MS.dt = g.MS*I.MS
    # 
    # dS.HS.dt = -b.HS*S.HS*(P) / (N.HS) * transmission_modifier[3]
    # dI.HS.dt = b.HS*S.HS*(P) / (N.HS) * transmission_modifier[3] - g.HS*I.HS
    # dR.HS.dt = g.HS*I.HS
    
    # maybe working version
    transmission_modifier.LS = (1 / (1 + exp(-l * (C.LS))) + offset)
    transmission_modifier.MS = (1 / (1 + exp(-l * (C.MS))) + offset)
    transmission_modifier.HS = (1 / (1 + exp(-l * (C.HS))) + offset)
    # transmission_modifier.LS = (1 / (1 + exp(-15 * (C.LS))) + offset)
    # transmission_modifier.MS = (1 / (1 + exp(-300 * (C.MS))) + offset)
    # transmission_modifier.HS = (1 / (1 + exp(-800 * (C.HS))) + offset)

    dS.LS.dt = -b.LS*S.LS*(P) / N.LS * transmission_modifier.LS
    dI.LS.dt = b.LS*S.LS*(P) / N.LS * transmission_modifier.LS - g.LS*I.LS
    dR.LS.dt = g.LS*I.LS

    dS.MS.dt = -b.MS*S.MS*(P) / N.MS * transmission_modifier.MS
    dI.MS.dt = b.MS*S.MS*(P) / N.MS * transmission_modifier.MS - g.MS*I.MS
    dR.MS.dt = g.MS*I.MS

    dS.HS.dt = -b.HS*S.HS*(P) / N.HS * transmission_modifier.HS
    dI.HS.dt = b.HS*S.HS*(P) / N.HS * transmission_modifier.HS - g.HS*I.HS
    dR.HS.dt = g.HS*I.HS

    # # maybe working version 2
    # # transmission_modifier.LS = (1 / (1 + exp(-l * (C.LS))) + offset)
    # # transmission_modifier.MS = (1 / (1 + exp(-l * (C.MS))) + offset)
    # # transmission_modifier.HS = (1 / (1 + exp(-l * (C.HS))) + offset)
    # transmission_modifier.LS = (1 / (1 + exp(-0.0001 * (C.LS))) + offset)
    # transmission_modifier.MS = (1 / (1 + exp(-0.0001 * (C.MS))) + offset)
    # transmission_modifier.HS = (1 / (1 + exp(-0.0001 * (C.HS))) + offset)
    # 
    # dS.LS.dt = -b.LS*S.LS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.LS
    # dI.LS.dt = b.LS*S.LS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.LS - g.LS*I.LS
    # dR.LS.dt = g.LS*I.LS
    # 
    # dS.MS.dt = -b.MS*S.MS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.MS
    # dI.MS.dt = b.MS*S.MS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.MS - g.MS*I.MS
    # dR.MS.dt = g.MS*I.MS
    # 
    # dS.HS.dt = -b.HS*S.HS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.HS
    # dI.HS.dt = b.HS*S.HS*(P) / (N.LS + N.MS + N.HS) * transmission_modifier.HS - g.HS*I.HS
    # dR.HS.dt = g.HS*I.HS
    
    
    # Update dP using the sum of infected individuals
    dP.dt = (I.LS + I.MS + I.HS) - P
    
    
    
    return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, dS.HS.dt, dI.HS.dt, dR.HS.dt, dP.dt)))
  })
}



# #test
# prop.group = c(prop.LS, prop.MS, prop.HS)
# transmission_modifier = (1 / (1 + exp(-lambda * (cover.tester))) + offset)
# transmission_modifier = (1 / (1 + exp(-lambda * (cover.tester/prop.group))) + offset)
# #test

area = 200 * 3 #meters squared. also scaled in a way that translates to CPCe-based cover (1:3 ratio)

lambda = as.numeric(8.0)
offset = 1 - 1 / (1 + exp(-lambda * 1.0))

# # #Nearshore
# N.tester = 172.7661
# # N.tester = 31.64699
# cover.tester = N.tester/area
# # cover.tester = 0.25
# # cover.tester = 0.05
# # prop.LS = 0.1511057435
# # prop.MS = 0.6442768576
# # prop.HS = 0.2046174336
# prop.LS = 0.00000000001
# prop.MS = 0.00000000001
# prop.HS = 0.2046174336
# # prop.LS = 0.01
# # prop.MS = 0.80
# # prop.HS = 0.20
# # prop.LS = 0.2511057435
# # prop.MS = 0.6442768576
# # prop.HS = 0.1046174336 # 0.01 #0.2046174336
# # prop.LS = 0.35
# # prop.MS = 0.60
# # prop.HS = 0.05
# # prop.LS = 0.65
# # prop.MS = 0.30
# # prop.HS = 0.05

# # Midchannel
# N.tester = 31.64699
# # N.tester = 172.7661
# cover.tester = N.tester/area
# # cover.tester = 0.05
# # cover.tester = 0.25
# prop.LS = 0.6784525163
# prop.MS = 0.2340849793
# prop.HS = 0.08746250433
# # prop.LS = 0.7928867609
# # prop.MS = 0.1602192545
# # prop.HS = 0.04689422333
# # prop.LS = 0.000000000001
# # prop.MS = 0.000000000001
# # prop.HS = 0.04689422333

#Offshore
# N.tester = 16.75861
N.tester = 31.64699
cover.tester = N.tester/area
prop.LS = 0.7928867609
prop.MS = 0.1602192545
prop.HS = 0.04689422333
# prop.LS = 0.95
# prop.MS = 0.01
# prop.HS = 0.04
# prop.LS = 0.10
# prop.MS = 0.85
# prop.HS = 0.05
# prop.LS = 0.1511057435
# prop.MS = 0.6442768576
# prop.HS = 0.2046174336
# prop.LS = 0.6784525163
# prop.MS = 0.2340849793
# prop.HS = 0.08746250433

N.LS = prop.LS*N.tester
N.MS = prop.MS*N.tester
N.HS = prop.HS*N.tester
cover.LS = N.LS / area
cover.MS = N.MS / area
cover.HS = N.HS / area

#polyp_SA information:
#Nearshore: 0.010189426 [/ 172.7661 = 0.00005897815602 or 6e-5]
#Midchannel: 0.0004273632 [/ 31.64699 = 0.000013504071 or 1e-5]
#Offshore: 0.0009611146 [/ 16.75861 = 0.00005735049625 or 6e-5]
polyp_SA = 6e-5 * N.tester #*prop.HS
# polyp_SA = 1e-4 #1e-2

#LS
beta.tester.LS = 0.06#/N.LS #0.003
gamma.tester.LS = 0.13 #0.02
# beta.tester.LS = 0.00804#/N.LS
# gamma.tester.LS = 0.08294
I.tester.LS = 0
S.tester.LS = N.LS - I.tester.LS
R.tester.LS = 0

#MS
beta.tester.MS = 0.04#/N.MS#/prop.MS #0.01
gamma.tester.MS = 0.09 #0.04
# beta.tester.MS = 0.09113#/N.MS
# gamma.tester.MS = 0.17553
I.tester.MS = 0
S.tester.MS = N.MS - I.tester.MS
R.tester.MS = 0

#HS
beta.tester.HS = 3.67#/N.HS#/prop.HS #0.30
gamma.tester.HS = 2.44 #0.2
# beta.tester.HS = 1.10503#/N.HSå
# gamma.tester.HS = 1.47514
I.tester.HS = polyp_SA
S.tester.HS = N.HS - I.tester.HS
R.tester.HS = 0

P.tester = I.tester.HS

days.tester = 500
time.tester = 1:days.tester
output.tester = data.frame(ode(c(S.LS = S.tester.LS, I.LS = I.tester.LS, R.LS = R.tester.LS,
                                 S.MS = S.tester.MS, I.MS = I.tester.MS, R.MS = R.tester.MS,
                                 S.HS = S.tester.HS, I.HS = I.tester.HS, R.HS = R.tester.HS,
                                 P = P.tester),
                               time.tester, SIR, c(b.LS = beta.tester.LS, g.LS = gamma.tester.LS,
                                                   b.MS = beta.tester.MS, g.MS = gamma.tester.MS,
                                                   b.HS = beta.tester.HS, g.HS = gamma.tester.HS,
                                                   prop.LS = prop.LS, prop.MS = prop.MS, prop.HS = prop.HS,
                                                   N = N.tester,
                                                   N.LS = N.LS, N.MS = N.MS, N.HS = N.HS,
                                                   C = cover.tester,
                                                   C.LS = cover.LS, C.MS = cover.MS, C.HS = cover.HS,
                                                   l = lambda
                               )))
output.tester = output.tester %>% select(-ncol(.))


betas.tester.LS.adj = beta.tester.LS * (1 / (1 + exp(-lambda * (cover.LS))) + offset)
betas.tester.MS.adj = beta.tester.MS * (1 / (1 + exp(-lambda * (cover.MS))) + offset)
betas.tester.HS.adj = beta.tester.HS * (1 / (1 + exp(-lambda * (cover.HS))) + offset)
betas.tester = c(beta.tester.LS, beta.tester.MS, beta.tester.HS)
betas.tester.adj = c(betas.tester.LS.adj, betas.tester.MS.adj, betas.tester.HS.adj)

gammas.tester = c(gamma.tester.LS, gamma.tester.MS, gamma.tester.HS)
R0s.tester = c((betas.tester.adj[1])/gammas.tester[1], (betas.tester.adj[2])/gammas.tester[2], (betas.tester.adj[3])/gammas.tester[3])
R0s.tester = round(R0s.tester, 2)
covers.tester = c(cover.LS, cover.MS, cover.HS)

# beta.tester.LS.adj = beta.tester.LS * (1 / (1 + exp(-lambda * (cover.LS))) + offset)
# beta.tester.MS.adj = beta.tester.LS * (1 / (1 + exp(-lambda * (cover.MS))) + offset)
# beta.tester.HS.adj = beta.tester.LS * (1 / (1 + exp(-lambda * (cover.HS))) + offset)
# 
# betas.tester = c(beta.tester.LS, beta.tester.MS, beta.tester.HS)
# betas.tester.adj = c(beta.tester.LS.adj, beta.tester.MS.adj, beta.tester.HS.adj)
# gammas.tester = c(gamma.tester.LS, gamma.tester.MS, gamma.tester.HS)
# R0s.tester = c((betas.tester.adj[1])/gammas.tester[1], (betas.tester.adj[2])/gammas.tester[2], (betas.tester.adj[3])/gammas.tester[3])
# R0s.tester = round(R0s.tester, 2)

suscat.names = c('LS', 'MS', 'HS')
tab.tester = tibble(suscat.names, betas.tester, round(betas.tester.adj,2), gammas.tester, R0s.tester, round(covers.tester, 2))
names(tab.tester) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')
# covers.tester = c(cover.LS, cover.MS, cover.HS)
# suscat.names = c('LS', 'MS', 'HS')
# tab.tester = tibble(suscat.names, betas.tester, betas.tester.adj, gammas.tester, R0s.tester, round(covers.tester, 2))
# names(tab.tester) = c('Category', 'beta', 'Adj. beta', 'gamma', 'R0', 'Cover')

output.tester = pivot_longer(output.tester, cols = -1, names_pattern = "(.*)(..)$", names_to = c("Compartment", "Category")) %>%
  mutate(Compartment = ifelse(Compartment == "", "value", Compartment)) %>%
  mutate(Compartment = ifelse(Compartment == 'S.', 'Susceptible',
                              ifelse(Compartment == 'I.', 'Infected',
                                     ifelse(Compartment == 'R.', 'Dead', Compartment))))

colnames(output.tester)[1] = 'days'
colnames(output.tester)[4] = 'prop'

p.fit.tester = ggplot(data = output.tester, aes(days, prop, colour = Compartment, linetype = Category)) +
  xlab("Day of observation period") +
  ylab("Surface area of tissue") +
  geom_line() +
  scale_color_brewer(name = 'Disease compartment', palette = 'Set2') +
  annotate(geom = "table", x = max(output.tester$days), y = max(output.tester$prop)*0.7, label = list(tab.tester),
           vjust = 1, hjust = 1) +
  theme_classic()

p.S.fit.tester = ggplot(data = output.tester %>% filter(Compartment == "Susceptible"), aes(days, prop, linetype = Category)) +
  xlab("Day of observation period") +
  ylab("Surface area of tissue") +
  geom_line() +
  theme_classic()

p.I.fit.tester = ggplot(data = output.tester %>% filter(Compartment == "Infected"), aes(days, prop, linetype = Category)) +
  xlab("Day of observation period") +
  ylab("Surface area of tissue") +
  geom_line() +
  theme_classic()

p.D.fit.tester = ggplot(data = output.tester %>% filter(Compartment == "Dead"), aes(days, prop, linetype = Category)) +
  xlab("Day of observation period") +
  ylab("Surface area of tissue") +
  geom_line() +
  theme_classic()

# p.fit.tester
# p.S.fit.tester
p.I.fit.tester
# p.D.fit.tester
# p.fit.tester;p.S.fit.tester;p.I.fit.tester;p.D.fit.tester
