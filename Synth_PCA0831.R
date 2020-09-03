library(Synth)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

# Read data
usa.shen <- read.csv("szcpr-0710.csv",header = T,fileEncoding = "GBK")
usa.shen$per.confirm <- usa.shen$confirmed/usa.shen$popestimate2019*100000 # Calculate the number of confirmed cases per 100,000 people
usa.shen$unit.num <- rep(1:69,each = 29)
usa.shen$day <- as.numeric(rep(1:29,times = 69))
# PCA
datebefore <- c("03/01/2020","03/02/2020","03/03/2020","03/04/2020",
                "01/19/2020","01/20/2020","01/21/2020","01/22/2020")
confirmbef <- usa.shen[usa.shen$fulldate %in% datebefore,c("unit.num","per.confirm","day")]
confirmbef <- spread(confirmbef, key = "day", value = "per.confirm")
confirmbef.cov <- princomp(confirmbef[,2:5],cor=FALSE)
screeplot(confirmbef.cov,type='lines')
summary(confirmbef.cov)
usa.shen$comp1 <- rep(confirmbef.cov$scores[,1],each=29)
usa.shen$comp2 <- rep(confirmbef.cov$scores[,2],each=29)
# Synth
index.usa = 1:69
usa.shen$ctname<- as.character(usa.shen$ctname)
dataprep.out.usashen<-
  dataprep(
    foo = usa.shen,
    predictors = c("density","Latitude","comp1","comp2"),
    predictors.op = "mean",
    dependent = "per.confirm",
    unit.variable = "unit.num",
    time.variable = "day",
    special.predictors = NULL,
    treatment.identifier = 69,
    controls.identifier = index.usa[-69],
    time.predictors.prior = c(1:4),
    time.optimize.ssr =  c(1:4),
    unit.names.variable = "ctname",
    time.plot =  c(1:16)
  )
synth.out.usashen <- synth(dataprep.out.usashen)
path.plot(synth.res = synth.out.usashen, dataprep.res = dataprep.out.usashen)
synth.tables.usashen <- synth.tab(
  dataprep.res = dataprep.out.usashen,
  synth.res = synth.out.usashen)
# Output weight
synthw <- synth.tables.usashen$tab.w
synthw[order(synthw$w.weights, decreasing = T),]
# write.csv(synthw[order(synthw$w.weights, decreasing = T),],"weight0711.csv")
synth.tables.usashen$tab.pred

line.usashen <- data.frame(per.confirm = c(as.numeric(dataprep.out.usashen$Y1plot), as.numeric(dataprep.out.usashen$Y0plot %*% synth.out.usashen$solution.w)), 
                           date = rep((1:16),2),
                           index = rep(c("Shenzhen", "Synthetic"), each = 16))
line.usashen$date <- as.Date(line.usashen$date,origin = "2020-01-18")
ggplot(data = line.usashen, mapping = aes(x = date,y = per.confirm, group = index)) +
  geom_line(aes(color = index), size = 1)+scale_color_manual(values=c("#00BFC4", "#F8766D"))+
  geom_vline(aes(xintercept = as.Date(5,origin = "2020-01-18")), linetype = "dashed")+
  scale_x_date(date_labels ="%m/%d",date_breaks = "1 days")+
   scale_y_continuous(breaks = 1:13)+
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(family = "Times-Roman"),
        plot.title = element_text(hjust = 0.5))+
   annotate(geom = "text", label = "Policy Implementation", x = as.Date(5,origin = "2020-01-18"), 
            y = 7.6, hjust = 1.05, vjust = 0) +
  labs(x = "Date", y = "COVID-19 cases (per 100,000)", 
       title =NULL,
       color = "")

# Placebo test
placeboline.usashen <- {}
for (i in index.usa) {
  dataprep.out<-
    dataprep(
      foo = usa.shen,
      predictors = c("density","Latitude","comp1","comp2"),
      predictors.op = "mean",
      dependent = "per.confirm",
      unit.variable = "unit.num",
      time.variable = "day",
      special.predictors = NULL,
      treatment.identifier = i,
      controls.identifier = index.usa[-i],
      time.predictors.prior = c(1:4),
      time.optimize.ssr =  c(1:4),
      unit.names.variable = "ctname",
      time.plot =  c(1:16)
    )
  synth.out <- synth(dataprep.out)
  placeboline.usashen <- c(placeboline.usashen,list(dataprep.out$Y1plot - dataprep.out$Y0plot %*% synth.out$solution.w))
}
placebo.usashen <- as.data.frame(placeboline.usashen)
placebo.usashen<-gather(placebo.usashen,area,per.confirm,X1:X69)
placebo.usashen$date = rep(1:16,69)
placebo.usashen$date = as.Date(placebo.usashen$date,origin = "2020-01-18")
placebo.usashen$group = c(rep("Others",68*16),rep("Shenzhen",16))
ggplot(data = placebo.usashen,mapping=aes(x = date,y = per.confirm, group = area))+geom_line(aes(color = group,alpha = group),size = 0.8)+
  scale_colour_manual(values = c('grey',"#00BFC4"))+scale_alpha_manual(values = c(0.5,1))+
   geom_vline(aes(xintercept = as.Date(5,origin = "2020-01-18")),linetype = "dashed")+
  annotate(geom = "text", label = "Policy Implementation", x = as.Date(5,origin = "2020-01-18"), 
           y = 7.6, hjust = 1.05, vjust = 0) +
  scale_x_date(date_labels ="%m/%d",date_breaks = "1 days")+
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(x = "Date", y = "Gap in COVID-19 cases (per 100,000)", 
       title = NULL,alpha = NULL,colour = NULL
       )
trend.usa <- spread(line.usashen,index,per.confirm)
trend.usa$Shenzhen <- trend.usa$Shenzhen*130.266
trend.usa$Synthetic <- floor(trend.usa$Synthetic*130.266)
trend.usa
# ggsave("shenzhen_0712.pdf",device=cairo_pdf,width=8,height=6)
