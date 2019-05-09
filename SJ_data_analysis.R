library(data.table)
library(dplyr)
library(stringr)
library(xCell)
library(survival)
library(ggplot2)

##### Step 1: download data #####
# create a folded called Data in your current working directory

# run the following command line

# cd Data

# wget -r  ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/ALL/mRNA-seq/Phase2/L3/expression/StJude/
  
# wget -r  ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/ALL/mRNA-seq/Phase2/L3/expression/BCCA/
  
# wget -r  ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/ALL/clinical/Phase2/
  

##### Step 2: load data #####
fn = list.files("Data/StJude_Phase2_SJ",pattern = "TARGET",full.names = T)

expr0 = lapply(fn, function(x){
  e1 = fread(x,data.table=F)
  e1 = e1%>%select(gene, fpkm_normalized)
  colnames(e1)[2]=basename(x)
  return(e1)
})

expr0 = Reduce(function(x,y){inner_join(x,y,by="gene")},expr0)
expr0 = data.frame(expr0[,-1],row.names = expr0[,1],check.names = F)

clin0 = read.csv("Data/Clinical_Phase2/TARGET_ALL_Validation_ClinicalData_20170525.csv")


##### Step 3: xCell analysis #####
dir.create("Result")
expr1 = xCellAnalysis(expr0,save.raw = T,file.name="Result/Xcell")
expr1 = fread("Result/Xcell_RAW.txt",data.table = F)
expr1 = data.frame(expr1[,-1],row.names = expr1[,1],check.names = F)
expr1 = cbind.data.frame("file_name"=colnames(expr1),t(expr1))
expr1 = cbind.data.frame(TARGET.USI=gsub("(-03A|-09A).*","",expr1$file_name),expr1)
expr1 = cbind.data.frame(supID=str_match(expr1$file_name,"(09A-01R|03A-01R|09B-01R)")[,2],expr1)

df1 = inner_join(clin0,expr1,by = "TARGET.USI")
df1 = df1%>%filter(WBC.at.Diagnosis<200)
df1$SurvObj <- Surv(time = as.integer(df1$Event.Free.Survival.Time.in.Days), 
                    event = (!df1$First.Event %in% c("None","Censored","")))

write.csv(df1,"Result/xCell_and_meta_data.csv",row.names = F)

##### Step 4: plot for macrophages #####
col = "Macrophages"
df1$bin = (df1[,col]>680)

KM<- survfit(SurvObj ~ bin, df1,conf.type = "log-log")

p = survdiff(df1$SurvObj ~ df1$bin)
p <- 1 - pchisq(p$chisq, length(p$n) - 1)

plot(KM,col = c("Black","Red"),lwd=2,main = col,mark=3,
     ylim = c(0.8,1),
     mark.time=df1$Event.Free.Survival.Time.in.Days[df1$First.Event %in% c("None","Censored","")])
legend("bottomright", legend=c("low", "high"),
       col=c("black","red"), lty = 1,cex=2)
text(1000,0.85, paste("p value =",round(p,3)),cex=2)

##### Step 5: plot for B cells #####
col = "B-cells" 
df1$bin = (df1[,col]>median(df1[,col]))

KM<- survfit(SurvObj ~ bin, df1,conf.type = "log-log")

p = survdiff(df1$SurvObj ~ df1$bin)
p <- 1 - pchisq(p$chisq, length(p$n) - 1)

plot(KM,col = c("Black","Red"),lwd=2,main = col,mark=3,
     ylim = c(0.8,1),
     mark.time=df1$Event.Free.Survival.Time.in.Days[df1$First.Event %in% c("None","Censored","")])
legend("bottomright", legend=c("low", "high"),
       col=c("black","red"), lty = 1,cex=2)
text(1000,0.85, paste("p value =",round(p,3)),cex=2)

##### Step 6: plot for NK cells #####
col = "NK cells"
df1$bin = (df1[,col]>median(df1[,col]))

KM<- survfit(SurvObj ~ bin, df1,conf.type = "log-log")

p = survdiff(df1$SurvObj ~ df1$bin)
p <- 1 - pchisq(p$chisq, length(p$n) - 1)

plot(KM,col = c("Black","Red"),lwd=2,main = col,mark=3,
     ylim = c(0.8,1),
     mark.time=df1$Event.Free.Survival.Time.in.Days[df1$First.Event %in% c("None","Censored","")])
legend("bottomright", legend=c("low", "high"),
       col=c("black","red"), lty = 1,cex=2)
text(1000,0.85, paste("p value =",round(p,3)),cex=2)

