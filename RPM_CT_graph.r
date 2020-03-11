if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
}


rpm_ct<-read.xlsx('/Users/gerbix/Documents/vikas/scratch/pavian_in/RPM_CT_values.xlsx', sheetIndex = 1)

rpm_ct_zero_removed<-rpm_ct[which(rpm_ct$SARS.CoV.2_RPM > 0),]
# for(i in 1:nrow(rpm_ct_zero_removed)){ 
#   if(rpm_ct_zero_removed$RdRp.gene.CT[i] == '-' | rpm_ct_zero_removed$E.gene.CT[i] == '-'){ 
#     rpm_ct_zero_removed$
#     }
#   }


colors<-c('#e03210',
  '#e87252',
  '#e5a592',
  '#d4d4d4',
  '#a8b3d9',
  '#7594dc',
  '#1277de') 

rpm_ct_zero_removed$RdRp.gene.CT<-as.numeric(as.character(rpm_ct_zero_removed$RdRp.gene.CT))

best_fit = lm((RdRp.gene.CT) ~ log(SARS.CoV.2_RPM - min(rpm_ct_zero_removed$RdRp.gene.CT + 1), base = 10), data=rpm_ct_zero_removed)
x<-summary(best_fit)
x$r.squared



RPM_CT_plot<-ggplot(rpm_ct_zero_removed, aes(x = rpm_ct_zero_removed$RdRp.gene.CT, y = rpm_ct_zero_removed$SARS.CoV.2_RPM)) + 
  geom_point(aes(color = rpm_ct_zero_removed$Sample, shape = rpm_ct_zero_removed$Gene), size = 4) + 
  theme_classic() + 
  xlab('RdRp gene CT') + 
  ylab('log10(RPM)') + 
  scale_color_manual(values = colors) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_smooth(method = "lm", formula= (y ~ log(x - min(rpm_ct_zero_removed$RdRp.gene.CT) + 1, base = 2)),color = 'black', se = FALSE) +
  scale_y_log10(limits = c(1,8000)) + 
  #ylim(c(0,max(rpm_ct_zero_removed$SARS.CoV.2_RPM)))
  xlim(c(0,max(rpm_ct_zero_removed$RdRp.gene.CT))) 
RPM_CT_plot

ggsave(plot = RPM_CT_plot, file = 'RPM_CT_plot.pdf', height = 5, width = 5)
