###Project: Household modeling
###Purpose: Simulate model from https://academic.oup.com/aje/article/186/12/1380/3865598
###Author: Josh Petrie
###Date: 11/29/2022

# ================== Set Up ==================
# Read in proxy community hazard data, one row for each date of study follow up including
# counts of reported influenza cases standardized to the peak (i.e proxy in peak week ==1)
com <- read.csv("hh_model_exampleComData.csv")
# plot(com$DAY,com$COM, type="l")

# Set number of households
n_hh <- 500

# Set minimum household size
min_hh_size <- 3

# Set maximum household size
max_hh_size <- 10

# randomly generate list of sizes of each household (lower prob skews smaller households; higher, larger)
hh_sizes <- rbinom(n_hh, size=max_hh_size-min_hh_size, prob=0.2) + min_hh_size
# hist(hh_sizes)

# get total population size as all household members
n_pop <- sum(hh_sizes)

# make vector of household ids - one element for each individual
# initialize hh_ids, give each individual in first household an id of 1
hh_ids <- rep(1,hh_sizes[1])
# add ids for remaining households 2:n_hh
for(i in 2:n_hh){
  tmp <- rep(i,hh_sizes[i])
  hh_ids <- append(hh_ids,tmp)
}
rm(tmp)

# make dataframe to store data about individuals with initial states
agents <- data.frame( "id" = seq( 1,n_pop ), # individual id
                      "hhid" = hh_ids, # household id
                      "status" = "S", # starting status, S == susceptible
                      "onset_date" = rep( NA,n_pop ), # onset date holder
                      "time_I" = 0, # time since infection holder
                      "source_I" = rep( NA,n_pop ), # community or hh infection holder
                      "exposed" = 0, # ever exposed indicator holder
                      "who_infected" = rep( NA,n_pop ), # who infected me holder
                      "first_hh_case" = rep( NA,n_pop ), # who was index case holder
                      "order" = 0, #indicator of order within chain of infection (index=1, infected by index=2, etc.)
                      stringsAsFactors = FALSE)

# make dataframe of parameter values taken from: 
# https://academic.oup.com/aje/article/186/12/1380/3865598
pars <- data.frame( "community_hazard" = 0.000675, #community scaling parameter
                    "weibull_alpha" = 3.30, #serial interval shape parameter
                    "weibull_gamma" = 2.87, #serial interval shape parameter
                    "household_hazard" = 0.12 ) #household scaling parameter

# how many days to run the model
days <- nrow( com )

# ================== Model function ==================
model <- function( agents, pars, days ){
  
  #Move through time
  for(d in 1:days){
    
    #community infections
    #for all individuals still susceptible...
    for( i in agents$id[ agents$status == "S" ] ){
      #sample random number
      cProb <- runif(1,0,1)
      #if sampled number is less than probability of infection from community...
      if( cProb < ( 1 - exp( -( com$COM[ d ] * pars$community_hazard ) ) ) ){
        #...then individual is infected. Update their data
        agents$status[ agents$id == i ] <- "I"
        agents$onset_date[ agents$id == i ] <- d
        agents$source_I[ agents$id == i ] <- "community"
        agents$who_infected[ agents$id == i ] <- 0
        agents$order[ agents$id == i ] <- 1
        
      }
      
    }
    
    #household infections
    #for all currently infected individuals...
    for( i in agents$id[ agents$status == "I" & agents$time_I > 0 ] ){
      
      #get the household id and time since infection for the current infected person
      current_house <- agents$hhid[agents$id==i]
      current_time_I <- agents$time_I[agents$id==i]
      
      #for each susceptible contact of the infected person
      for(j in agents$id[ agents$hhid == current_house & agents$status == "S" ]){
        
        #indicate they have been exposed
        agents$exposed[ agents$id == j ] <- 1
        #if it hasn't been set already, make the first household case the current infected person
        if( is.na(agents$first_hh_case[ agents$id == j ]) ){
          agents$first_hh_case[ agents$id == j ] <- i
        }
        
        #sample random number
        hProb <- runif(1,0,1)
        #if sampled number is less than probability of infection from household...
        if( hProb < (1-exp(-(exp(-(current_time_I/pars$weibull_alpha)^pars$weibull_gamma)-
                  exp(-((current_time_I+1)/pars$weibull_alpha)^pars$weibull_gamma))*
                pars$household_hazard))){
          #...then individual is infected. Update their data
          agents$status[ agents$id == j ] <- "I"
          agents$onset_date[ agents$id == j ] <- d
          agents$source_I[ agents$id == j ] <- "household"
          agents$who_infected[ agents$id == j ] <- i
          agents$order[ agents$id == j ] <- agents$order[ agents$id == i ] + 1
          
        }
        
      }
      
    }
    
    #advance time since infection counter for all infected individuals
    agents$time_I[ agents$status == "I" ] = agents$time_I[ agents$status == "I" ] + 1
    
    #if time since infection is greater than 14 days then change status to recovered
    agents$status[ agents$status == "I" & 
                     agents$time_I > 14 ] <- "R"
  
  }
  
  return(agents)
  
}

# ================== Run the Model ==================

### Run it once
# out <- model( agents, pars, days )
# table(out$status)
# table(out$who_infected)
# table(out$source_I)
# hist(out$onset_date[!is.na(out$onset_date)])

### Multiple simulations
ptm <- proc.time()
runs <- 100 #100 runs takes about 2 minutes
for(r in 1:runs){
  if(r %in% c(seq(1:(runs/5))*5)){
    message(r)
    message((proc.time() - ptm)[3])
  }
  
  if(r==1){
    Out1 <- model( agents, pars, days )
    Out1$iter <- 1
  } else{
    tmp <- model( agents, pars, days )
    tmp$iter <- r
    Out1 <- rbind(Out1,tmp)
  }
}

rm(tmp)
proc.time() - ptm

# ================== Some basic summaries ==================
library(dplyr) 

###summarizing average number of infections by community or household source
Out1_summary <- Out1 %>%
  group_by(iter) %>%
  summarise(
    infections = sum(!is.na(onset_date)),
    cmnty_acquired = sum(source_I == "community", na.rm = TRUE),
    hh_exposed = sum(exposed),
    hh_acquired = sum(source_I == "household", na.rm = TRUE)
  ) %>% 
  mutate(
    infection_risk = infections/n_pop,
    #
    # I have questions about how best to actually recover SIR here.
    # It seems that the calculation below still assumes all infected by the
    # index case resulting in overestimation.
    #
    secondary_infection_risk = hh_acquired / hh_exposed
  ) %>%
  ungroup()

quantile(Out1_summary$infections, probs = c(.5,.05,.95))
quantile(Out1_summary$cmnty_acquired, probs = c(.5,.05,.95))
quantile(Out1_summary$hh_acquired, probs = c(.5,.05,.95))

###epi curves
library(ggplot2)

iter_week_summary <- Out1 %>%
  filter(!is.na(onset_date)) %>%
  mutate(week = floor(onset_date/7)) %>%
  group_by(iter, week) %>%
  tally() %>%
  ungroup()

all_iter_weeks <- data.frame( iter = rep(seq(1,100,1),each=15),
                         week = rep(seq(1,15,1), 100) )

iter_week_summary <- left_join(all_iter_weeks, iter_week_summary, by=c("iter","week")) %>%
  mutate(n = ifelse(is.na(n),0,n),
         iter = as.character(iter))

all_week_summary <- iter_week_summary %>%
  group_by(week) %>%
  summarize(n = median(n)) %>%
  mutate(iter="MEDIAN") %>%
  ungroup()

iter_week_summary <- bind_rows(iter_week_summary, all_week_summary)

ggplot(data=iter_week_summary, aes(x=week, y=n, group=iter, color=iter)) +
  geom_line(aes(size=iter)) +
  scale_size_manual(values=c(rep(1,100),2)) +
  scale_color_manual(values=c(rep("gray",100),"black")) +
  theme_classic() + 
  theme(legend.position="none") 
  
###plotting transmission networks
library(igraph)

#plot only households with at least one infection from the first run
nodes <- data.frame("id" = Out1$id[Out1$iter==1 & (Out1$status != "S" | Out1$exposed == 1)],
                    "order" = factor(Out1$order[Out1$iter==1 & (Out1$status != "S" | Out1$exposed == 1)]))
edges <- data.frame(
  "to" = ifelse(
    is.na(Out1$who_infected[Out1$iter==1]),
    Out1$first_hh_case[Out1$iter==1],
    Out1$who_infected[Out1$iter==1]
),
"from" = Out1$id[Out1$iter==1]
) %>%
  filter(!is.na(to) & to != 0 & to != from)

net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T) 

colrs <- c("white","black","yellow","orange","red") #,"#FDE725FF"
V(net)$color <- colrs[as.numeric(V(net)$order)+1]

plot(net, edge.arrow.size=.05, vertex.size=2, vertex.label=NA)
plot(net, edge.arrow.size=.01, vertex.size=7, vertex.label.cex=.6)

table(Out1$order,Out1$iter)
