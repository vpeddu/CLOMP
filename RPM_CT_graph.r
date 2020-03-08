if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
}


rpm_ct<-read.xlsx('/Users/gerbix/Documents/vikas/scratch/pavian_in/RPM_CT_values.xlsx', sheetIndex = 1)

rpm_ct_zero_removed<-rpm_ct[which(rpm_ct$SARS.CoV.2_RPM > 0),]

colors<-c('#e03210',
  '#e87252',
  '#e5a592',
  '#d4d4d4',
  '#a8b3d9',
  '#7594dc',
  '#1277de') 

RPM_CT_plot<-ggplot(rpm_ct_zero_removed, aes(x = rpm_ct_zero_removed$SARS.CoV.2_avg_CT, y = rpm_ct_zero_removed$SARS.CoV.2_RPM)) + 
  geom_point(aes(color = rpm_ct_zero_removed$Sample)) + 
  theme_classic() + 
  xlab('CT') + 
  ylab('RPM') + 
  scale_color_manual(values = colors) +
  theme(legend.position="bottom", legend.title = element_blank())
RPM_CT_plot

ggsave(plot = RPM_CT_plot, file = 'RPM_CT_plot.pdf', height = 5, width = 5)
