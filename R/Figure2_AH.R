library(data.table)
library(PAFGRS)
library(ggplot2)
library(gridExtra)


dt=data.table(h2=c(.7,.7,.4,.5,.5),K=c(.01,.05,.125,.5,.5),method=c(rep(1,4),2),Nsib=rep(c(1:10),each=5) )
est_r_fgrs_emp = apply(dt,1,function(x) fgrs_accuracy(nrel = x[4],prev = x[2],h2 = x[1],rel_rr =.5,rel_ir=.5,N=100000,method = c("pa","quantitative")[x[3]],estimate = "simulate" ))
dt= cbind(dt,t(est_r_fgrs_emp))


dt$r_fgrs_theory= apply(dt,1,function(x) fgrs_accuracy(nrel = x[4],prev = x[2],h2 = x[1],rel_rr =.5,rel_ir=.5,N=100000,method = c("pa","quantitative")[x[3]],estimate = "theory" )[1])
dt$r_fgrs_theory2 =dt$r_fgrs_theory

dt[,class:=paste0("K= ",K,", h2= ", h2)]
dt[method==2,class:=paste0("quantitative",", h2= ", h2)]

g1 <- ggplot(data=dt, aes(x=factor(Nsib),y=r_fgrs_theory,color=class,group=class))+
  geom_line()+
  geom_point(aes(x=factor(Nsib),y = accuracy.estimate,color=class,group=class))+
  geom_errorbar(aes(ymin=accuracy.CI.1,ymax=accuracy.CI.2,width=.25,color=class,group=class))+
  labs(x=expression(N["full siblings"]),y=expression(cor(g,hat(g)[fgrs])),color="")+
  geom_hline(yintercept = sqrt(.5), linetype="dashed")+theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+coord_cartesian(ylim=0:1)+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+labs(x=expression(N["full siblings"]),y=expression(r[paste(G,",",hat(G)[PA-FGRS])]),color="prevalence")+expand_limits(x=0)+
  scale_x_discrete(breaks=seq(0,10,2))


g2 <- ggplot(data=dt, aes(x=factor(Nsib),y=r_fgrs_theory^2*h2,color=class,group=class))+geom_line()+geom_point(aes(x=factor(Nsib),accuracy.estimate^2*h2))+geom_errorbar(aes(ymin=accuracy.CI.1^2*h2,ymax=accuracy.CI.2^2*h2,width=.25))+labs(x=expression(N["full siblings"]),y=expression(R[paste(l,",",hat(G)[PA-FGRS])]^2),color="prevalence")+theme_bw()+
  geom_hline(yintercept = .35,
             linetype="3131",color=scales::hue_pal()(4)[1])+
  geom_hline(yintercept = .35,
             linetype="35",color=scales::hue_pal()(4)[2])+
  geom_hline(yintercept = .2,
             linetype="35",color=scales::hue_pal()(4)[3])+
  geom_hline(yintercept = .25,
             linetype="3131",color=scales::hue_pal()(4)[4])+
  geom_hline(yintercept = .25,
             linetype="35",color="black")+coord_cartesian(ylim=c(0,.75))+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
  scale_x_discrete(breaks=seq(0,10,2))

### Repeat for offspring 

dt=data.table(h2=c(.7,.7,.4,.5,.5),K=c(.01,.05,.125,.5,.5),method=c(rep(1,4),2),Nsib=rep(c(1:10),each=5) )
est_r_fgrs_emp = apply(dt,1,function(x) fgrs_accuracy(nrel = x[4],prev = x[2],h2 = x[1],rel_rr =.25,rel_ir=.5,N=100000,method = c("pa","quantitative")[x[3]],estimate = "simulate" ))
dt= cbind(dt,t(est_r_fgrs_emp))


dt$r_fgrs_theory= apply(dt,1,function(x) fgrs_accuracy(nrel = x[4],prev = x[2],h2 = x[1],rel_rr =.25,rel_ir=.5,N=100000,method = c("pa","quantitative")[x[3]],estimate = "theory" )[1])
dt$r_fgrs_theory2 =dt$r_fgrs_theory

dt[,class:=paste0("K= ",K,", h2= ", h2)]
dt[method==2,class:=paste0("quantitative",", h2= ", h2)]

g3 <- ggplot(data=dt, aes(x=factor(Nsib),y=r_fgrs_theory,color=class,group=class))+
  geom_line()+
  geom_point(aes(x=factor(Nsib),y = accuracy.estimate,color=class,group=class))+
  geom_errorbar(aes(ymin=accuracy.CI.1,ymax=accuracy.CI.2,width=.25,color=class,group=class))+
  geom_hline(yintercept = 1, linetype="dashed")+theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+coord_cartesian(ylim=0:1)+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+
  labs(x=expression(N["half-siblings offspring"]),y=expression(r[paste(G,",",hat(G)[PA-FGRS])]),color="prevalence")+expand_limits(x=0)+
  scale_x_discrete(breaks=seq(0,10,2))


g4 <- ggplot(data=dt, aes(x=factor(Nsib),y=r_fgrs_theory^2*h2,color=class,group=class))+geom_line()+geom_point(aes(x=factor(Nsib),accuracy.estimate^2*h2))+geom_errorbar(aes(ymin=accuracy.CI.1^2*h2,ymax=accuracy.CI.2^2*h2,width=.25))+labs(x=expression(N["half-siblings offspring"]),y=expression(R[paste(l,",",hat(G)[PA-FGRS])]^2),color="prevalence")+theme_bw()+
  geom_hline(yintercept = .7,
             linetype="3131",color=scales::hue_pal()(4)[1])+
  geom_hline(yintercept = .7,
             linetype="35",color=scales::hue_pal()(4)[2])+
  geom_hline(yintercept = .4,
             linetype="35",color=scales::hue_pal()(4)[3])+
  geom_hline(yintercept = .5,
             linetype="3131",color=scales::hue_pal()(4)[4])+
  geom_hline(yintercept = .5,
             linetype="35",color="black")+coord_cartesian(ylim=c(0,.75))+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
  scale_x_discrete(breaks=seq(0,10,2))


### Generate a pedigree 

dt=data.table(h2=c(.7,.7,.4,.5),K=c(.01,.05,.125,.5),Nchild=rep(c(1:2),each=4),Ngen=rep(c(1:5),each=8))
est_r_fgrs = apply(dt,1,function(x) PAFGRS:::r_fgrs_gen_2nd(ngen = x[4],child_pr_gen = x[3],prev = x[2],h2 = x[1]))

r_fgrs_gen_quant= function(ngen=2, child_pr_gen=2, prev=0.2 , h2=.5){
  ped=PAFGRS:::generate_pedigree(ngen,child_pr_gen)
  r_fgrs_pedigree_quant(h2=h2,ped=ped)
}
r_fgrs_pedigree_quant=
  function(h2=.5, ped=NULL,rel_matrix=NULL){
    if(is.null(rel_matrix))
      rel <- kinship(pedigree(id=ped$id,dadid = ped$dadid,momid = ped$momid,sex=ped$sex))*2
    else rel <- rel_matrix
    #P <- h2*rel[-1,-1]*(dnorm(qnorm(1-prev)))^2
    P <- h2*rel[-1,-1]
    diag(P) <- 1
    G <- h2*rel[-1,1]
    c(r=sqrt(t(G)%*%solve(P)%*%G/h2),n_rel=sum(G>0))}
est_r_fgrs_q = apply(dt[h2==0.5],1,function(x) r_fgrs_gen_quant(ngen = x[4],child_pr_gen = x[3],h2 = 0.5))

dt_q=dt[h2==0.5]
dt_q[,K:=0]
dt= rbind(cbind(dt,t(est_r_fgrs)),cbind(dt_q,t(est_r_fgrs_q)))

dt[,class:=paste0("K= ",K,", h2= ", h2)]
dt[K==0,class:="quantitative"]


g5.1 <- ggplot(data=dt[Nchild==1], aes(x=factor(Ngen+1),y=r,color=class,group=class,label=n_rel))+
  geom_line()+geom_point()+
  labs(x=expression(paste("Generations (",N[rel],")")),y=expression(r[paste(G,",",hat(G)[PA-FGRS])]),color="")+theme_bw()+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+expand_limits(y=c(0,.5))+
  scale_x_discrete(labels=paste0(dt[Nchild==1&K==0,Ngen+1]," (",dt[Nchild==1&K==0,n_rel],")"))

g5.2 <- ggplot(data=dt[Nchild==1], aes(x=factor(Ngen+1),y=r^2*h2,color=class,group=class,label=n_rel))+
  geom_line()+geom_point()+
  labs(x=expression(paste("Generations (",N[rel],")")),y=expression(R[paste(l,",",hat(G)[PA-FGRS])]^2),color="")+theme_bw()+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+expand_limits(y=c(0,.17))+
  scale_x_discrete(labels=paste0(dt[Nchild==1&K==0,Ngen+1]," (",dt[Nchild==1&K==0,n_rel],")"))


g6.1 <- ggplot(data=dt[Nchild==2], aes(x=factor(Ngen+1),y=r,color=class,group=class,label=n_rel))+
  geom_line()+geom_point()+
  labs(x=expression(paste("Generations (",N[rel],")")),y=expression(r[paste(G,",",hat(G)[PA-FGRS])]),color="")+theme_bw()+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+expand_limits(y=c(0,.5))+
  scale_x_discrete(labels=paste0(dt[Nchild==1&K==0,Ngen+1]," (",dt[Nchild==2&K==0,n_rel],")"))

g6.2 <- ggplot(data=dt[Nchild==2], aes(x=factor(Ngen+1),y=r^2*h2,color=class,group=class,label=n_rel))+
  geom_line()+geom_point()+
  labs(x=expression(paste("Generations (",N[rel],")")),y=expression(R[paste(l,",",hat(G)[PA-FGRS])]^2),color="")+theme_bw()+scale_color_manual(values = c(scales::hue_pal()(4),"black"))+theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+expand_limits(y=c(0,.17))+
  scale_x_discrete(labels=paste0(dt[Nchild==1&K==0,Ngen+1]," (",dt[Nchild==2&K==0,n_rel],")"))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(g1+ theme(legend.position="bottom"))

grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                         g2 + theme(legend.position="none"),
                         g3 + theme(legend.position="none"),
                         g4 + theme(legend.position="none"),
                         g5.1 + theme(legend.position="none"),
                         g5.2 + theme(legend.position="none"),
                         g6.1 + theme(legend.position="none"),
                         g6.2 + theme(legend.position="none"),ncol = 2),
             mylegend,
             nrow=2,heights=c(10,1))

ggsave(arrangeGrob(g1 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g2 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g3 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g4 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g5.1 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g5.2 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g6.1 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),
                   g6.2 + theme(legend.position="none",aspect.ratio = 1,axis.text=element_text(size=7)),ncol = 2),filename = "fig2_AH.pdf",height = 10,width=6)
