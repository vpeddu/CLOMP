if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
}

RPM_summary<-read.xlsx('/Users/gerbix/Documents/vikas/scratch/pavian_in/RPM_summary.xlsx', sheetIndex = 1)

#Staphyloccus epidermidis = 1282
# Streptococcus = 1301 
# Cutibacterium acnes = 1747 
# Moraxella catarrhalis = 480 

common_flora_taxid<-c('1282', '1301', '1747')

flora_df<-RPM_summary[RPM_summary$taxid %in% common_flora_taxid,][,c(3:ncol(RPM_summary))]

common_flora_list<-c()
for(i in 1:nrow(flora_df)){ 
  for(j in 1:ncol(flora_df)){ 
    if(flora_df[i,j] > 10){
    common_flora_list<-append(common_flora_list, flora_df[i,j])
      }
    }
  }


flora_median<-median(as.numeric(common_flora_list))


moraxella_df<-RPM_summary[which(RPM_summary$taxid == '480'),]

t.test(common_flora_list, mu = moraxella_df$WA6.UW3[1])

t.test(common_flora_list, mu = mean(moraxella_df$WA9.UW6[1], moraxella_df$WA8.UW[1]))




