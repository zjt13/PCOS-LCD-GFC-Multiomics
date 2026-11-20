#!/usr/bin/env Rscript
# Figure 1E — body-weight trajectories; mean +/- SEM
# Run from ANY directory: Rscript scripts/01_bodyweight.R
args_all <- commandArgs(trailingOnly = FALSE)
script_flag <- grep('^--file=', args_all, value = TRUE)
if (length(script_flag) != 1L) stop('Use: Rscript scripts/01_bodyweight.R')
repo <- normalizePath(file.path(dirname(sub('^--file=', '', script_flag)), '..'), mustWork = TRUE)
input <- file.path(repo, 'data', 'bodyweight_long.csv')
out <- file.path(repo, '.work', 'Fig1', 'bodyweight')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (pkg in c('ggplot2', 'dplyr', 'readr')) if (!requireNamespace(pkg, quietly = TRUE)) stop('Install missing R package: ',pkg)
library(ggplot2);library(dplyr);library(readr)
groups <- c('Control','PCOS','LCD','GFC','Combine')
group_cols <- c(Control='#4D4D4D',PCOS='#D55E00',LCD='#0072B2',GFC='#009E73',Combine='#7B2CBF')
dat <- read_csv(input, show_col_types=FALSE, col_types=readr::cols(Group=col_character(),Mouse=col_character(),Week=col_integer(),Weight_g=col_double()))
stopifnot(all(c('Group','Mouse','Week','Weight_g') %in% names(dat)), nrow(dat)>0L)
if (anyNA(dat) || any(!dat$Group %in% groups) || any(!dat$Week %in% 1:4) || any(dat$Weight_g<=0)) stop('Invalid body weight data')
if (anyDuplicated(dat[c('Group','Mouse','Week')])) stop('Duplicate animal-week records')
dat <- dat %>% mutate(Group=factor(Group, levels=groups), Subject=interaction(Group, Mouse, drop=TRUE))
summary <- dat %>% group_by(Group, Week) %>% summarise(n=n(), mean=mean(Weight_g),sd=sd(Weight_g),sem=sd/sqrt(n), .groups='drop')
write_csv(summary,file.path(out,'bodyweight_summary_mean_SEM.csv'))
p_mean <- ggplot(summary, aes(Week,mean,color=Group,group=Group))+
 geom_ribbon(aes(ymin=mean-sem,ymax=mean+sem,fill=Group),alpha=.10,color=NA,show.legend=FALSE)+
 geom_line(linewidth=1.1)+geom_point(size=2.6)+
 scale_color_manual(values=group_cols,drop=FALSE)+scale_fill_manual(values=group_cols,drop=FALSE)+
 scale_x_continuous(breaks=1:4, labels=paste('Week',1:4))+
 labs(x='Time',y='Body weight (g)',color=NULL)+theme_classic(base_size=13)+theme(legend.position='right')
p_individual <- ggplot(dat,aes(Week,Weight_g,group=Subject))+
 geom_line(color='grey60',linewidth=.5,alpha=.75)+geom_point(color='grey40',size=1.5)+
 geom_line(data=summary,aes(x=Week,y=mean,group=1),inherit.aes=FALSE,linewidth=1.0)+
 geom_point(data=summary,aes(x=Week,y=mean),inherit.aes=FALSE,shape=21,fill='white',size=2.4)+
 facet_wrap(~Group,nrow=1,scales='free_y')+
 scale_x_continuous(breaks=1:4)+labs(x='Week',y='Body weight (g)')+theme_classic(base_size=12)
base <- dat %>% group_by(Group,Mouse) %>% mutate(Baseline=Weight_g[Week==1][1],Delta_g=Weight_g-Baseline) %>% ungroup()
if (anyNA(base$Baseline)) stop('Each mouse must have a Week 1 measurement')
delta <- base %>% group_by(Group,Week) %>% summarise(n=n(),mean=mean(Delta_g),sd=sd(Delta_g),sem=sd/sqrt(n),.groups='drop')
write_csv(delta,file.path(out,'bodyweight_delta_summary.csv'))
p_delta <- ggplot(delta,aes(Week,mean,color=Group,group=Group))+
 geom_hline(yintercept=0,linetype='dashed',color='grey65')+
 geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=.07)+geom_line(linewidth=.9)+geom_point(size=2.3)+
 scale_color_manual(values=group_cols)+scale_x_continuous(breaks=1:4,labels=paste('Week',1:4))+
 labs(x='Time',y='Change in body weight (g)',color=NULL)+theme_classic(base_size=12)
saveplot <- function(p,n,w,h){figdir<-file.path(repo,'figures');dir.create(figdir,recursive=TRUE,showWarnings=FALSE);ggsave(file.path(figdir,paste0(n,'.pdf')),p,width=w,height=h,bg='white');ggsave(file.path(figdir,paste0(n,'.png')),p,width=w,height=h,dpi=300,bg='white')}
saveplot(p_mean,'Fig1E',4.9,3.3)
# Additional diagnostics are computed but not exported as manuscript panels.

message('Fig1E written to ',file.path(repo,'figures'))
