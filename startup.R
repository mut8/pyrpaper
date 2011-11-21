source("~/R/functions.R")

addpeaks<-function(source, sep=",") {

  for (i in 1:length(source))
    {
      print(i)
      tmp<-read.csv(source[i], sep=sep, header=F)
      pe<-data.frame(t(tmp[1:10,2:ncol(tmp)]))
      colnames(pe)<-tmp[1:10, 1]
      rownames(pe)<-paste(pe$class, 1:nrow(pe), sep="")

      si<-tmp[12:85, 2:ncol(tmp)]
      colnames(si)<-paste(rownames(pe), "sim", sep="")
      rownames(si)<-rownames(samples)

      ti<-tmp[102:175, 2:ncol(tmp)]
      colnames(ti)<-paste(rownames(pe), "tic", sep="")
      rownames(ti)<-rownames(samples)

      if (i==1)
        {
          peaks<-pe
          sim<-si
          tic<-ti
        } else
      {
        peaks<-rbind(peaks, pe)
        sim<-data.frame(sim, si)
        tic<-data.frame(tic, ti)
      }
    }
  write.csv(peaks, "raw_data/peaks.csv")
  write.csv(sim, "raw_data/sim.csv")
  write.csv(tic, "raw_data/tic.csv")
  print("files added")
}

samples<-read.csv("raw_data/samples.csv", sep=";", header=TRUE)
addpeaks(c("raw_data/furane.csv","raw_data/aromates.csv", "raw_data/phenole.csv", "raw_data/carbos_neu.csv", "raw_data/cyclopentenone.csv", "raw_data/n_compounds_neu.csv", "raw_data/alkan_en_cleaned.csv", "raw_data/the_rest.csv", "raw_data/lignin2.csv", "raw_data/keto-alkyl_alcohole.csv", "raw_data/fa.csv"), sep=";")

minerals<-read.csv("raw_data/minerals.csv", sep=",", header=TRUE)
alldata<-read.csv("raw_data/micdif_alldata.csv", sep=";", header=TRUE)
alldata$days<-rep(c(rep(14,5), rep(97,5), rep(181,5), rep(475,5)),4)
f2b<-read.csv("raw_data/micdif_f2b.csv", sep=",", header=TRUE)


sim<-read.csv("raw_data/sim.csv", sep=",", header=TRUE)
rownames(sim)<-sim[,1]
sim[,1]<-NULL
peaks<-read.csv("raw_data/peaks.csv", sep=",", header=TRUE)

corr<-read.csv("raw_data/corr_data.csv")
corr$pH.CaCl2<-NULL

                                        #remove excluded
sim<-sim[,peaks$exclude==F]
peaks<-peaks[peaks$exclude==F,T]

days<-samples$day
accresp<-samples$accresp
harvest<-samples$harvest
type<-samples$type
harlev<-as.numeric(levels(as.factor(harvest)))
typlev<-levels(type)

peaknr<-length(colnames(sim))

                                        # calculate cTIC
cTIC<-sim
for (i in 1:peaknr)
cTIC[,i]<-sim[,i] * peaks$tic[i]

                                        #subs by categories
colnames(peaks)

rsim<-100*sim/rowSums(sim)
rcTIC<-100*cTIC/rowSums(cTIC)
class_rsim<-sumif(rsim,peaks$class)
orig_rsim<-sumif(rsim,peaks$origin)
class_cTIC<-sumif(rcTIC,peaks$class)
orig_cTIC<-sumif(rcTIC,peaks$origin)



samples$massloss->massloss
samples$masslossSE->masslossSE

processes<-data.frame(corr[,38], corr[,3:30], corr[,51:56], corr[,61:69])
controls<-data.frame(class_cTIC[harvest==2|harvest==6, T], orig_cTIC[harvest==2|harvest==6, T],  corr[,31:34], corr[,39:44], corr[,35:37], corr[,57:59], corr[,45:50], corr[,60])

contnames<-colnames(controls)
procnames<-colnames(processes)

h2<-samples$harvest==2
h3<-samples$harvest==6
h23<-harvest==2|harvest==6
h234<-harvest==2|harvest==6|harvest==15
h34<-harvest==6|harvest==15
h4<-harvest==15

minerals.stat<-data.frame(matrix(nrow=length(typlev), ncol=2*(ncol(minerals[,2:ncol(minerals)]))))
colnames(minerals.stat)<-paste(c(rep("means", ncol(minerals[,2:ncol(minerals)])), rep("se", ncol(minerals[,2:ncol(minerals)]))), rep(colnames(minerals[,2:ncol(minerals)]),2))
rownames(minerals.stat)<-typlev

for (i in 1:length(typlev))
  for (j in 1:(ncol(minerals)-1))
{
  minerals.stat[i,j]<-mean(minerals[minerals[,1]==typlev[i], j+1])
  minerals.stat[i,j+(ncol(minerals)-1)]<-stderr(minerals[minerals[,1]==typlev[i], j+1])
}


alldata.stat<-data.frame(matrix(nrow=length(typlev)*length(levels(alldata$Harvest)), ncol=2+ncol(alldata[3:ncol(alldata)])*2))
colnames(alldata.stat)<-c("type", "harvest", paste(c(rep("means", ncol(alldata[,3:ncol(alldata)])), rep("se", ncol(alldata[,3:ncol(alldata)]))), rep(colnames(alldata[,3:ncol(alldata)]),2)))
alldata.stat$type<-sort(rep(typlev, length(harlev)))
alldata.stat$harvest<-rep(levels(alldata$Harvest), length(typlev))


for (i in 1:nrow(alldata.stat))
  for (k in 3:(ncol(alldata)))
{
  alldata.stat[i,k] <-
    mean(alldata[alldata$Litter==alldata.stat$type[i] & alldata$Harvest==alldata.stat$harvest[i], k])
  alldata.stat[i,k-1+(ncol(alldata.stat)/2)] <-
    stderr(alldata[alldata$Litter==alldata.stat$type[i] & alldata$Harvest==alldata.stat$harvest[i], k])
}

colscale<-c(grey(0), grey(.2), grey(.4), grey(.6))
colscale.all<-c(rep(colscale[1], 4), rep(colscale[2], 5), rep(colscale[3], 5), rep(colscale[4], 5),
                rep(colscale[1], 4), rep(colscale[2], 4), rep(colscale[3], 5), rep(colscale[4], 5),
                rep(colscale[1], 4), rep(colscale[2], 5), rep(colscale[3], 5), rep(colscale[4], 5),
                rep(colscale[1], 4), rep(colscale[2], 5), rep(colscale[3], 5), rep(colscale[4], 4))
pch<-21:24
pch.all<-c(rep(pch[1],19),rep(pch[2],18),rep(pch[3],19),rep(pch[4],18))

initials.rsim<-rsim
initials.orig_cTIC <- orig_cTIC
initials.class_cTIC<- class_cTIC

for (i in 1:nrow(samples))
{
  initials.rsim[i,T] <-      colMeans(rsim      [harvest==0 & type==type[i], T])
  initials.orig_cTIC[i,T] <- colMeans(orig_cTIC [harvest==0 & type==type[i], T])
  initials.class_cTIC[i,T] <-colMeans(class_cTIC[harvest==0 & type==type[i], T])
}


h3init.rsim<-rsim
h3init.orig_cTIC <- orig_cTIC
h3init.class_cTIC<- class_cTIC
samples$h3init.respacc.litC <- samples$cons_acc_resp_litC

closs.corr.rsim <- rsim*(1- samples$cons_acc_resp_litC)
closs.corr.orig_cTIC <- orig_cTIC*(1- samples$cons_acc_resp_litC)
closs.corr.class_cTIC <- class_cTIC* (1- samples$cons_acc_resp_litC)


for (i in 1:nrow(samples))
{
  h3init.rsim[i,T] <-      colMeans(closs.corr.rsim      [harvest==6 & type==type[i], T])
  h3init.orig_cTIC[i,T] <- colMeans(closs.corr.orig_cTIC [harvest==6 & type==type[i], T])
  h3init.class_cTIC[i,T] <-colMeans(closs.corr.class_cTIC[harvest==6 & type==type[i], T])
  samples$h3init.respacc.litC[i]<-mean(samples$cons_acc_resp_litC[harvest==6 & type==type[i]])
}



alldata$CN_inbal <- alldata$C.N_lit/ alldata$C.N_mic
alldata$CP_inbal <- alldata$C.P_lit/ alldata$C.Pmic
alldata$NP_inbal <- alldata$N.P_lit/ alldata$N.Pmic


alld.h4.cond<-alldata$days==475 & c(rep(T,79),F)

alldata$phen2cell<-alldata$Phenoloxidase/alldata$Cellulase
alldata$per2cell<-alldata$Peroxydase/alldata$Cellulase


hmw.proc<-data.frame(
  orig_cTIC$L[h3]-initials.orig_cTIC$L[h3],
  orig_cTIC$C[h3]-initials.orig_cTIC$C[h3],
  orig_cTIC$L[h3] / (orig_cTIC$L[h3]+orig_cTIC$C[h3]) - initials.orig_cTIC$L[h3]/(initials.orig_cTIC$L[h3]+initials.orig_cTIC$C[h3]),
  (-closs.corr.orig_cTIC$L[h3]+initials.orig_cTIC$L[h3])/initials.orig_cTIC$L[h3],
  (-closs.corr.orig_cTIC$C[h3]+initials.orig_cTIC$C[h3])/initials.orig_cTIC$C[h3],
  (-closs.corr.orig_cTIC$L[h3]+initials.orig_cTIC$L[h3])/(initials.orig_cTIC$L[h3]*samples$cons_acc_resp_litC[h3]),
  (-closs.corr.orig_cTIC$C[h3]+initials.orig_cTIC$C[h3])/(initials.orig_cTIC$C[h3]*samples$cons_acc_resp_litC[h3]),
  ((-closs.corr.orig_cTIC$L[h3]+initials.orig_cTIC$L[h3])/initials.orig_cTIC$L[h3])/((-closs.corr.orig_cTIC$C[h3]+initials.orig_cTIC$C[h3])/initials.orig_cTIC$C[h3]), 
  alldata$phen2cell[alldata$days==181]
  , alldata$per2cell[alldata$days==181])

colnames(hmw.proc) <- c("LTIC", "CTIC", "LCI", "Ldec", "Cdec", "Lresp", "Cresp", "LCdec", "Phen2Cell", "Per2Cell"
      #, "log(L:Cdec)"
      )

notlig.proc<-data.frame(
  class_cTIC$al0[h3]-initials.class_cTIC$al0[h3],
  class_cTIC$al1[h3]-initials.class_cTIC$al1[h3],
  class_cTIC$fa[h3]-initials.class_cTIC$fa[h3],
  class_cTIC$phytol[h3]-initials.class_cTIC$phytol[h3],
  (-closs.corr.class_cTIC$al0[h3]+initials.class_cTIC$al0[h3])/initials.class_cTIC$al0[h3],
  (-closs.corr.class_cTIC$al1[h3]+initials.class_cTIC$al1[h3])/initials.class_cTIC$al1[h3],
  (-closs.corr.class_cTIC$fa[h3]+initials.class_cTIC$fa[h3])/initials.class_cTIC$fa[h3],
  (-closs.corr.class_cTIC$phytol[h3]+initials.class_cTIC$phytol[h3])/initials.class_cTIC$phytol[h3],
  (-closs.corr.class_cTIC$al0 [h3]+initials.class_cTIC$al0[h3]) / (initials.class_cTIC$al0[h3]*samples$cons_acc_resp_litC[h3]),
  (-closs.corr.class_cTIC$al1 [h3]+initials.class_cTIC$al1[h3]) / (initials.class_cTIC$al1[h3]*samples$cons_acc_resp_litC[h3]),
  (-closs.corr.class_cTIC$fa [h3]+initials.class_cTIC$fa[h3]) / (initials.class_cTIC$fa[h3]*samples$cons_acc_resp_litC[h3]),
  (-closs.corr.class_cTIC$phytol [h3]+initials.class_cTIC$phytol[h3]) / (initials.class_cTIC$phytol[h3]*samples$cons_acc_resp_litC[h3])
)

colnames(notlig.proc) <-   c("alkanacc","alkenacc", "faacc", "phytolacc", "alkandeg", "alkendeg", "fadeg", "phytoldeg", "alkanresp", "alkenresp", "faresp", "phytolresp")

hmw.proc2 <-
data.frame(  orig_cTIC$L[h4]-h3init.orig_cTIC$L[h4],
orig_cTIC$C[h4]-h3init.orig_cTIC$C[h4],
orig_cTIC$L[h4] / (orig_cTIC$L[h4]+orig_cTIC$C[h4]) - h3init.orig_cTIC$L[h4]/(h3init.orig_cTIC$L[h4]+h3init.orig_cTIC$C[h4]),
(-closs.corr.orig_cTIC$L[h4]+h3init.orig_cTIC$L[h4])/h3init.orig_cTIC$L[h4],
(-closs.corr.orig_cTIC$C[h4]+h3init.orig_cTIC$C[h4])/h3init.orig_cTIC$C[h4],
(closs.corr.orig_cTIC$L[h4]-h3init.orig_cTIC$L[h4])/h3init.orig_cTIC$L[h4]*
  -1/(samples$cons_acc_resp_litC[h4]-samples$h3init.respacc.litC[h4]),
(closs.corr.orig_cTIC$C[h4]-h3init.orig_cTIC$C[h4])/h3init.orig_cTIC$C[h4]*
  -1/(samples$cons_acc_resp_litC[h4]-samples$h3init.respacc.litC[h4]),
((-closs.corr.orig_cTIC$L[h4]+h3init.orig_cTIC$L[h4])/h3init.orig_cTIC$L[h4])/((-closs.corr.orig_cTIC$C[h4]+h3init.orig_cTIC$C[h4])/h3init.orig_cTIC$C[h4]), 
alldata$phen2cell[alld.h4.cond], alldata$per2cell[alld.h4.cond])

colnames(hmw.proc2) <- c("LTIC", "CTIC", "LCI", "Ldec", "Cdec", "Lresp", "Cresp", "LCdec", "Phen2Cell", "Per2Cell"
      #, "log(L:Cdec)"
      )

notlig.proc2<-data.frame(
  class_cTIC$al0[h4]-h3init.class_cTIC$al0[h4],
  class_cTIC$al1[h4]-h3init.class_cTIC$al1[h4],
  class_cTIC$fa[h4]-h3init.class_cTIC$fa[h4],
  class_cTIC$phytol[h4]-h3init.class_cTIC$phytol[h4],
  (-closs.corr.class_cTIC$al0[h4]+h3init.class_cTIC$al0[h4])/h3init.class_cTIC$al0[h4],
  (-closs.corr.class_cTIC$al1[h4]+h3init.class_cTIC$al1[h4])/h3init.class_cTIC$al1[h4],
  (-closs.corr.class_cTIC$fa[h4]+h3init.class_cTIC$fa[h4])/h3init.class_cTIC$fa[h4],
  (-closs.corr.class_cTIC$phytol[h4]+h3init.class_cTIC$phytol[h4])/h3init.class_cTIC$phytol[h4],
  (-closs.corr.class_cTIC$al0[h4]+initials.class_cTIC$al0[h4])/(initials.class_cTIC$al0[h4]*samples$h3init.respacc.litC[h4]),
  (-closs.corr.class_cTIC$al1[h4]+initials.class_cTIC$al1[h4])/(initials.class_cTIC$al1[h4]*samples$h3init.respacc.litC[h4]),
  (-closs.corr.class_cTIC$fa[h4]+initials.class_cTIC$fa[h4])/(initials.class_cTIC$fa[h4]*samples$h3init.respacc.litC[h4]),
  (-closs.corr.class_cTIC$phytol[h4]+initials.class_cTIC$phytol[h4])/(initials.class_cTIC$phytol[h4]*samples$h3init.respacc.litC[h4])
)
colnames(notlig.proc2) <-   c("alkanacc","alkenacc", "faacc", "phytolacc", "alkandeg", "alkendeg", "fadeg", "phytoldeg", "alkanresp", "alkenresp", "faresp", "phytolresp")
