plot_ed5<-function(df,p){
  dodge_width<-.8
  ggplot(df,aes(x=type,y=value,fill=Distance))+
    geom_bar(stat="identity",position=position_dodge(width=dodge_width),width=.7,linewidth=.5)+
    geom_errorbar(aes(ymin=value-Dse,ymax=value+Dse),position=position_dodge(width=dodge_width),width=.3,linewidth=.5)+
    facet_wrap(rho~.,nrow=1)+coord_cartesian(ylim=c(.5,.75))+
    theme_few()+ylab("C-index")+xlab("")+
    scale_fill_npg(labels=c("DISCO","DM"))+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=10),
      axis.title=element_text(size=12),legend.text=element_text(size=10),
      legend.title=element_text(size=12),strip.text=element_text(size=14),
      text=element_text(family="Helvetica"))+
    labs(title=paste0(p," variables, Noise ~ N(0,1)"))
}
