if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
}
if (!require("Rsamtools")) { 
  BiocManager::install("Rsamtools")
  library('Rsamtools')
}
if (!require("data.table")) { 
  install.packages("data.table")
  library('data.table')
}
if (!require("tidyr")) { 
  install.packages("tidyr")
  library('tidyr')
}
if (!require("cowplot")) { 
  install.packages("cowplot")
  library('cowplot')
}
if (!require("ggplot2")) { 
  install.packages("ggplot2")
  library('ggplot2')
}



RPM_summary<-read.xlsx('/Users/gerbix/Documents/vikas/scratch/pavian_in/RPM_summary.xlsx', sheetIndex = 1)

signifcant_list<-c('1747','480','11216','694009',)
CA<-'1747'
MC<-'480'
HPIV3<-'11216'
CV<-'694009'
RVA<-'147711'
RVC<-'463676'

colnames(RPM_summary)<-c('taxid','name','WA6-UW3', 'WA7-UW4', 'WA4-UW2', 'SC5683', 'WA3-UW1', 'SC5698', 'WA9-UW6', 'WA8-UW')


tax_pull<-function(df, tax_class){ 
  keep_list<-which(df$taxid == tax_class)
  temp_df<-df[keep_list,]
  t_df<-t(temp_df)
  final_temp<-t_df[3:nrow(t_df),]
  final_temp<-as.data.frame(final_temp)
  colnames(final_temp)[1]<-'RPM'
  final_temp$RPM<-as.character(final_temp$RPM)
  return(final_temp)
  }

CA_df<-(tax_pull(RPM_summary, CA))
MC_df<-tax_pull(RPM_summary, MC)
HPIV3_df<-tax_pull(RPM_summary, HPIV3)
CV_df<-tax_pull(RPM_summary, CV)
RVA_df<-tax_pull(RPM_summary, RVA)
RVC_df<-tax_pull(RPM_summary, RVC)


colors<-c('#003f5c','#444e86','#955196','#dd5182','#ff6e54','#ffa600')

CA_plot<-ggplot(CA_df, aes(x = rownames(CA_df), y = as.numeric(CA_df$RPM))) + 
  geom_bar(stat = "identity",fill=colors[1]) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(expand = c(0, 0),limits = c(1,1000))+ 
  xlab('Sample') + 
  ylab('RPM') + 
  labs(title="Cutibacterium acnes") + 
  theme(plot.title = element_text(size=8))
CA_plot

MC_plot<-ggplot(MC_df, aes(x = rownames(MC_df), y = as.numeric(MC_df$RPM))) + 
  geom_bar(stat = "identity",fill=colors[2]) + 
  labs(title="Moraxella catarrhalis")  + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  #scale_y_continuous(expand = c(0, 0),limits = c(0,6000))+ 
  xlab('Sample') + 
  ylab('RPM') +
  scale_y_log10(expand = c(0,0),limits = c(1,10000))+
  theme(plot.title = element_text(size=8))
MC_plot

HPIV3_plot<-ggplot(HPIV3_df, aes(x = rownames(HPIV3_df), y = as.numeric(HPIV3_df$RPM))) + 
  geom_bar(stat = "identity",fill=colors[3]) + 
  labs(title="HPIV3")  + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(expand = c(0, 0),limits = c(1,10000))+ 
  xlab('Sample') + 
  ylab('RPM') +
  theme(plot.title = element_text(size=8))

HPIV3_plot


CV_plot<-ggplot(CV_df, aes(x = rownames(CV_df), y = as.numeric(CV_df$RPM))) + 
  geom_bar(stat = "identity",fill=colors[4]) + 
  labs(title="SARS-CoV-2 ")  + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(expand = c(0, 0),limits = c(1,10000))+ 
  xlab('Sample') + 
  ylab('RPM') +
  theme(plot.title = element_text(size=8))

CV_plot

RVA_plot<-ggplot(RVA_df, aes(x = rownames(RVA_df), y = as.numeric(RVA_df$RPM))) + 
  geom_bar(stat = "identity",fill=colors[5]) + 
  labs(title="Rhinovirus A")  + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(expand = c(0, 0),limits = c(1,10000))+ 
  xlab('Sample') + 
  ylab('RPM') +
  theme(plot.title = element_text(size=8))

RVA_plot

RVC_plot<-ggplot(RVC_df, aes(x = rownames(RVC_df), y = as.numeric(RVC_df$RPM))) + 
  geom_bar(stat = "identity",fill=colors[6]) + 
  labs(title="Rhinovirus C")  + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(expand = c(0, 0),limits = c(1,10000))+ 
  xlab('Sample') + 
  ylab('RPM') +
  theme(plot.title = element_text(size=8))

RVC_plot

panel<-plot_grid(CV_plot, MC_plot, HPIV3_plot, RVA_plot, RVC_plot,CA_plot,  labels = c('A', 'B', 'C', 'D', 'E', 'F'), label_size = 12)

ggsave(plot = panel, filename = 'plot.pdf',height = 10, width = 10)

save_plot("p4.pdf", p4, ncol = 2, base_asp = 1.1)





