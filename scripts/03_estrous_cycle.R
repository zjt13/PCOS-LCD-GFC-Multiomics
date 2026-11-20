#!/usr/bin/env Rscript
# Figure 1C,D: reads externally supplied OBSERVED cycle-stage data.
args_all <- commandArgs(trailingOnly=FALSE)
script_flag <- grep('^--file=',args_all,value=TRUE)
if(length(script_flag)!=1L)stop('Use: Rscript scripts/03_estrous_cycle.R [--demo]')
repo <- normalizePath(file.path(dirname(sub('^--file=','',script_flag)),'..'),mustWork=TRUE)
demo <- '--demo' %in% commandArgs(trailingOnly=TRUE)
input <- file.path(repo,'data','estrous_cycle.csv')
out <- file.path(repo,'.work','Fig1',if(demo)'estrous_DEMONSTRATION_ONLY' else 'estrous')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(pkg in c('ggplot2','readr','dplyr'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install missing package: ',pkg)
library(ggplot2);library(readr);library(dplyr)
groups <- c('Control','PCOS','LCD','GFC','Combine')
dat <- read_csv(input,show_col_types=FALSE,col_types=cols(Group=col_character(),Mouse=col_character(),Day=col_integer(),Stage=col_character()))
if(anyNA(dat) || any(!dat$Group %in% groups) || any(!dat$Stage %in% c('P','E','M','D')) || any(dat$Day<1) || anyDuplicated(dat[c('Group','Mouse','Day')]))stop('Invalid estrous observation data')
if(!demo && nrow(dat)<=40)warning('Only ',nrow(dat),' records. Verify that these are genuine observations, not the demonstration dataset.')
dat <- dat %>% mutate(Group=factor(Group,levels=groups),Stage=factor(Stage,levels=c('M','E','P','D')))
# Representative selection: first lexically sorted Mouse ID within each group; disclose in README.
ids <- dat %>% distinct(Group,Mouse) %>% arrange(Group,Mouse) %>% group_by(Group) %>% slice(1) %>% ungroup()
repdat <- inner_join(dat,ids,by=c('Group','Mouse')) %>% mutate(StageY=case_when(Stage=='D'~1,Stage=='M'~2,Stage=='E'~3,Stage=='P'~4))
stage <- dat %>% count(Group,Stage,name='n_observations',.drop=FALSE) %>% group_by(Group) %>% mutate(Percent=100*n_observations/sum(n_observations)) %>% ungroup()
write_csv(stage,file.path(out,'estrous_stage_percentage.csv'))
write_csv(ids,file.path(out,'representative_mouse_ids.csv'))
p_line <- ggplot(repdat,aes(Day,StageY,group=interaction(Group,Mouse)))+geom_line(linewidth=.65)+geom_point(size=2)+facet_wrap(~Group,nrow=1)+
 scale_y_continuous(breaks=1:4,labels=c('D','M','E','P'),limits=c(.75,4.25))+scale_x_continuous(breaks=sort(unique(dat$Day)))+
 labs(x='Time (days)',y=NULL)+theme_classic(base_size=12)+theme(strip.text=element_text(face='bold'),panel.spacing.x=grid::unit(1,'lines'))
stage_colors <- c(M='#AF25C9',E='#28A745',P='#DD3333',D='#315ED1')
p_bar <- ggplot(stage,aes(Group,Percent,fill=Stage))+geom_col(width=.7,color='white',linewidth=.25)+scale_fill_manual(values=stage_colors,drop=FALSE)+
 scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),expand=c(0,0))+
 labs(x=NULL,y='Estrous-stage percentage (%)',fill=NULL)+theme_classic(base_size=13)+theme(axis.text.x=element_text(angle=30,hjust=1))
if(demo){p_line<-p_line+labs(subtitle='DEMONSTRATION DATA — NOT EXPERIMENTAL OBSERVATIONS')+theme(plot.subtitle=element_text(color='red',face='bold',size=10));p_bar<-p_bar+labs(subtitle='DEMONSTRATION DATA — NOT EXPERIMENTAL OBSERVATIONS')+theme(plot.subtitle=element_text(color='red',face='bold',size=10))}
saveplot <- function(p,n,w,h){figdir<-if(demo)file.path(out,'figures') else file.path(repo,'figures');dir.create(figdir,recursive=TRUE,showWarnings=FALSE);ggsave(file.path(figdir,paste0(n,'.pdf')),p,width=w,height=h,bg='white');ggsave(file.path(figdir,paste0(n,'.png')),p,width=w,height=h,dpi=300,bg='white')}
saveplot(p_line,'Fig1C',10.5,3)
saveplot(p_bar,'Fig1D',5.4,4)
message(if(demo)'DEMO plots saved (NOT VALID FOR MANUSCRIPT): ' else 'Observed-data plots saved: ',out)
