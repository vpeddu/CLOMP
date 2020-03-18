if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
}
if (!require("cowplot")) { 
  install.packages("cowplot")
  library('cowplot')
}
if(!require("vridis")){ 
  install.packages("viridis")
  library("viridis")
  
  }


#PANEL A 
rpm_ct<-read.xlsx('/Users/gerbix/Documents/vikas/scratch/draft_6/RPM_CT_values.xlsx', sheetIndex = 1)

rpm_ct_zero_removed<-rpm_ct[which(rpm_ct$SARS.CoV.2_RPM > 0),]


colors<-c('#e03210',
  '#e87252',
  '#e5a592',
  '#d4d4d4',
  '#a8b3d9',
  '#7594dc',
  '#1277de') 

rpm_ct_zero_removed$RdRp.gene.CT<-as.numeric(as.character(rpm_ct_zero_removed$RdRp.gene.CT))
rpm_ct_zero_removed$Sample<-as.character(rpm_ct_zero_removed$Sample)

#rpm_ct_zero_removed<-rpm_ct_zero_removed[-which(rpm_ct_zero_removed$Sample == 'WA3-UW1'),]

#best_fit = lm((RdRp.gene.CT) ~ log(SARS.CoV.2_RPM - min(rpm_ct_zero_removed$RdRp.gene.CT + 1), base = 10), data=rpm_ct_zero_removed)

rl <- lm(SARS.CoV.2_RPM ~ RdRp.gene.CT , data = rpm_ct_zero_removed)
x<-summary(rl)
x$r.squared



RPM_CT_plot<-ggplot(rpm_ct_zero_removed, aes(x = rpm_ct_zero_removed$RdRp.gene.CT, y = rpm_ct_zero_removed$SARS.CoV.2_RPM)) + 
  geom_point(aes(color = rpm_ct_zero_removed$Sample, shape = rpm_ct_zero_removed$Gene), size = 4) + 
  theme_classic() + 
  ggtitle('RPM of assigned SARS-CoV-2 reads') + 
  xlab('RdRp gene CT') + 
  ylab('RPM') + 
  scale_color_viridis(discrete=TRUE) +
  #scale_color_manual(values = colors) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  theme(plot.title = element_text(size = 8)) + 
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  #geom_smooth(method = "lm", formula= (y ~ log(x - min(rpm_ct_zero_removed$RdRp.gene.CT) + 1, base = 2)),color = 'black', se = FALSE) +
  scale_y_log10(limits = c(1,8000)) + 
  #ylim(c(0,max(rpm_ct_zero_removed$SARS.CoV.2_RPM)))
  xlim(c(10,40)) 
RPM_CT_plot

ggsave(plot = RPM_CT_plot, file = 'Figure_2_draft_8.pdf', height = 5, width = 5)


# PANEL B 
# 
# reads_assembly<-read.xlsx('/Users/gerbix/Documents/vikas/scratch/pavian_in/RPM_CT_values.xlsx', sheetIndex = 2)
# reads_assembly<-reads_assembly[complete.cases(reads_assembly),]
# 
# reads_assembly$Read_counts<-as.numeric(as.character(reads_assembly$Read_counts))
# reads_assembly$percent_genome_cover<-as.numeric(as.character(reads_assembly$percent_genome_cover))
# 
# reads_assembly_zero_removed<-reads_assembly[-which(reads_assembly$percent_genome_cover==0),]
# 
# reads_assembly$Read_counts<-as.numeric(as.character(reads_assembly$Read_counts))
# 
# reads_assembly_plot<-ggplot(reads_assembly_zero_removed, aes(x = reads_assembly_zero_removed$Sample, y = reads_assembly_zero_removed$percent_genome_cover)) + 
#   geom_bar(aes(fill = reads_assembly_zero_removed$Read_counts), stat = 'identity') + 
#   theme_classic() + 
#   xlab('Sample') + 
#   ylab('% Assembled') +
#   #scale_y_continuous(expand=c(0,100)) + 
#   #ylim(c(0,0)) + 
#   scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
#   #scale_color_gradient(high = "#7d376a", low = '#f4b3df', aesthetics = 'fill') +
#   #scale_color_manual(values = colors) +
#   theme(legend.position="bottom") +
#   scale_fill_viridis(option="cividis", trans = 'reverse',limits=c(2e7,9e5), breaks = c(2e7, 1e7 ,9e5)) +
#   ggtitle('Percent of SARS-CoV-2 genome assembled from unassigned reads') + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(plot.title = element_text(size = 8)) + 
#   labs(fill = 'Total reads on sample' ) 
#   #scale_fill_brewer(palette = 2)
# #scale_fill_continuous(high = "#9ad4a5", low = '#1db85b', aesthetics = 'fill', trans = 'reverse') 
#   #scale_x_log10()
#   #geom_smooth(method = "lm", formula= (y ~ log(x - min(reads_assembly_zero_removed$RdRp.gene.CT) + 1, base = 2)),color = 'black', se = FALSE) +
#   #scale_y_log10(limits = c(1,8000)) + 
#   #ylim(c(0,max(rpm_ct_zero_removed$SARS.CoV.2_RPM)))
#   #xlim(c(0,max(reads_assembly_zero_removed$RdRp.gene.CT))) 
# reads_assembly_plot
# 
# panel <-plot_grid(reads_assembly_plot,RPM_CT_plot, labels = c('A', 'B'), label_size = 12)

ggsave(plot = panel, filename = 'figure_2_draft_6.pdf', height = 5, width = 8)

