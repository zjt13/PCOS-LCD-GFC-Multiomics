#!/usr/bin/env Rscript
# Figure 1G,H — one observation per animal; two-sided Welch t-tests (UNADJUSTED P values)
args_all <- commandArgs(trailingOnly=FALSE)
script_flag <- grep('^--file=',args_all,value=TRUE)
if(length(script_flag)!=1L)stop('Use: Rscript scripts/02_HE_quantification.R')
repo <- normalizePath(file.path(dirname(sub('^--file=','',script_flag)),'..'),mustWork=TRUE)
infile <- file.path(repo,'data','HE_raw_counts.csv')
out <- file.path(repo,'.work','Fig1','histology');dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(pkg in c('ggplot2','dplyr','tidyr','ggpubr','readr'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install missing package: ',pkg)
library(ggplot2);library(dplyr);library(tidyr);library(ggpubr);library(readr)
groups <- c('Control','PCOS','LCD','GFC','Combine')
colors <- c(Control='#3155E7',PCOS='#B18459',LCD='#9C38A6',GFC='#FF922B',Combine='#F97B7D')
he <- read_csv(infile,show_col_types=FALSE,col_types=cols(Group=col_character(),Mouse=col_character(),Corpus_Lutea=col_double(),Ovarian_Cysts=col_double()))
if(anyNA(he) || any(!he$Group %in% groups) || any(he$Corpus_Lutea<0) || any(he$Ovarian_Cysts<0) || anyDuplicated(he[c('Group','Mouse')]))stop('Invalid HE animal-level input')
he$Group <- factor(he$Group,levels=groups)
long <- he %>% pivot_longer(cols=c('Corpus_Lutea','Ovarian_Cysts'),names_to='Outcome',values_to='Count')
sm <- long %>% group_by(Outcome,Group) %>% summarise(n=n(),mean=mean(Count),sd=sd(Count),sem=sd/sqrt(n),.groups='drop')
write_csv(sm,file.path(out,'HE_summary_mean_SEM.csv'))
comparisons <- c('Control','LCD','GFC','Combine')
stats <- bind_rows(lapply(unique(long$Outcome),function(o){bind_rows(lapply(comparisons,function(g){
 x <- long$Count[long$Outcome==o & long$Group==g]
 y <- long$Count[long$Outcome==o & long$Group=='PCOS']
 tst <- t.test(x,y,var.equal=FALSE,alternative='two.sided')
 data.frame(Outcome=o,group1='PCOS',group2=g,n_PCoS=length(y),n_group=length(x),difference=mean(x)-mean(y),t=unname(tst$statistic),df=unname(tst$parameter),p_unadjusted=tst$p.value)
}))}))
stats$sig <- ifelse(stats$p_unadjusted<.0001,'****',ifelse(stats$p_unadjusted<.001,'***',ifelse(stats$p_unadjusted<.01,'**',ifelse(stats$p_unadjusted<.05,'*','ns'))))
write_csv(stats,file.path(out,'HE_Welch_vs_PCOS_unadjusted.csv'))
make_plot <- function(outcome,ylabel,y_start,y_step){
 sub <- long %>% filter(Outcome==outcome)
 sub_sm <- sm %>% filter(Outcome==outcome)
 ss <- stats %>% filter(Outcome==outcome) %>% mutate(y.position=y_start+(match(group2,comparisons)-1)*y_step)
 ggplot()+geom_col(data=sub_sm,aes(Group,mean,fill=Group),width=.62,alpha=.92)+
 geom_errorbar(data=sub_sm,aes(Group,ymin=pmax(0,mean-sem),ymax=mean+sem,color=Group),width=.16,linewidth=.6,show.legend=FALSE)+
 geom_jitter(data=sub,aes(Group,Count,color=Group),width=.08,height=0,size=2.4,show.legend=FALSE)+
 scale_fill_manual(values=colors,guide='none')+scale_color_manual(values=colors,guide='none')+
 stat_pvalue_manual(ss,label='sig',xmin='group1',xmax='group2',y.position='y.position',tip.length=.01,bracket.size=.5,size=4,hide.ns=FALSE)+
 coord_cartesian(ylim=c(0,max(ss$y.position)+.8),clip='off')+
 labs(x=NULL,y=ylabel)+theme_classic(base_size=13)+theme(axis.text.x=element_text(angle=45,hjust=1,color='black'),plot.margin=margin(10,10,8,8))
}
p_l <- make_plot('Corpus_Lutea','Number of Corpora Lutea',7.2,.72)
p_c <- make_plot('Ovarian_Cysts','Number of Ovarian Cysts',5.3,.62)
p_both <- ggarrange(p_l,p_c,ncol=1,nrow=2,align='v')
saveplot <- function(p,n,w,h){figdir<-file.path(repo,'figures');dir.create(figdir,recursive=TRUE,showWarnings=FALSE);ggsave(file.path(figdir,paste0(n,'.pdf')),p,width=w,height=h,bg='white');ggsave(file.path(figdir,paste0(n,'.png')),p,width=w,height=h,dpi=300,bg='white')}
saveplot(p_l,'Fig1G',4,4)
saveplot(p_c,'Fig1H',4,4)
# Combined diagnostic plot not exported as a separate manuscript panel.
message('Histology plots saved to ',out)
