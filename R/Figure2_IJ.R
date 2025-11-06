library(PAFGRS)

h2s=c(.7,.7,.4,.5)
prevs=c(.01,.05,.125,.5)
out <- lapply(1:4, function(i) {
  h2=h2s[i]
  prev=prevs[i]
  dt = rbind(data.table("siblings",1:25,prev=prev , h2=h2, rel_ir=0.5,rel_rr=0.5),
             data.table("half-siblings",1:25,prev=prev , h2=h2, rel_ir=0.25,rel_rr=0.25),
             data.table("cousins",1:25,prev=prev , h2=h2, rel_ir=0.125,rel_rr=0.125),
             data.table("offspring",1:25,prev=prev , h2=h2, rel_ir=0.5,rel_rr=0.25),
             data.table("parents",1:2,prev=prev , h2=h2, rel_ir=0.5,rel_rr=0),
             data.table("MZ",1,prev=prev , h2=h2, rel_ir=1,rel_rr=1),
             data.table("g.parents",1:4,prev=prev , h2=h2, rel_ir=.25,rel_rr=0),
             data.table("g.g.parents",1:8,prev=prev , h2=h2, rel_ir=.125,rel_rr=0))
  dt=cbind(dt,t(apply(dt[,-1],1,function(x) {fgrs_accuracy(nrel = x[1],prev = x[2],h2 = x[3],rel_ir = x[4],rel_rr = x[5],method = "linear2",estimate = "theory")})))
  
  
  dt_gen = data.table("mixed",prev=prev , h2=h2, ngen=1:5,child_pr_gen=rep(1:2,each=5))
  dt_gen= cbind(dt_gen,t(apply(dt_gen[,-1],1,function(x) {fgrs_accuracy(ped=generate_pedigree(ngen = x[3],child_pr_gen = x[4]),prev = x[1],h2 = x[2],,method = "linear2",estimate = "theory")})))
  #dt_gen$n_rel = apply(dt_gen[,-1],1,function(x) { nrow(generate_pedigree(ngen = x[3],child_pr_gen = x[4]))-1})
  
  dt=rbind(dt,dt_gen[,.(paste("child_pr_gen=",child_pr_gen),n_rel=accuracy.n_rel,prev,h2,NA,NA,acc=accuracy.r,sib_eq=sibling_equivalents.prev)],use.names=F)
  dt}) 
dt=Reduce('rbind',out)

names(dt)[7:8] <- c("acc","sib_eq")
library(ggplot2)

### plot the sib_eq for class2: 
g2=ggplot(data=dt[h2==.7&prev==.01],aes(x=factor(V2),y=sib_eq,color=V1,group=V1))+geom_line()+geom_point()+ scale_y_continuous(limits = c(0,10.1))+ #scale_x_continuous(breaks=1:10)+
  theme_bw()+labs(x=expression(N[rel]),y="Sibling equivalents",color="Relative(s)")

g4=ggplot(data=dt[h2==.5&prev==.5],aes(x=factor(V2),y=sib_eq,color=V1,group=V1))+geom_line()+geom_point()+ scale_y_continuous(limits = c(0,10.1))+ #scale_x_continuous(breaks=1:10)+
  theme_bw()+labs(x=expression(N[rel]),y="Sibling equivalents",color="Relative(s)")

g2 
g4 

