# Read Kraken tsvs and give back excel file with RPM calculations 

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
args = commandArgs(trailingOnly=TRUE)

path<-args[1]

setwd(path)

files<-list.files(pattern = '*.tsv')

taxa_detect<-function(df, taxid){ 
  temp_rpm<-df$RPM[which(df$taxid == taxid)]
  if(identical(temp_rpm, numeric(0))){ 
    temp_rpm <- 0
  }
  return(temp_rpm)
}

T1_RPM<-c()
MS2_RPM<-c()
for(i in 1:length(files)){ 
  print(files[i])
  temp_tsv<-read.csv(files[i], sep = "\t", col.names = c('percent_clade_reads', 'number_clade_reads_rooted_at_taxon','number_clade_reads_this_taxon', 'taxa', 'taxid', 'name'), header = FALSE)
  total_reads = temp_tsv$number_clade_reads_rooted_at_taxon[2] + temp_tsv$number_clade_reads_rooted_at_taxon[1]
  temp_tsv$RPM = temp_tsv$number_clade_reads_this_taxon  / (total_reads / 1e6)
  temp_tsv$cumulative_RPM<-temp_tsv$number_clade_reads_rooted_at_taxon / (total_reads / 1e6)
  temp_tsv$taxa<-trimws(temp_tsv$taxa , which = "both", whitespace = "\t")
  temp_tsv$name<-as.character(temp_tsv$name)
  file_name = strsplit(files[i],"_final")[[1]][1]
  if( i == 1 ){ 
    final_tsv<-temp_tsv[,c(5,6,7)] 
    colnames(final_tsv)[3]<-file_name
  }
  else{ 
    final_tsv[,(i+2)]<-0
    new_list<-c()
    for(j in 1:nrow(temp_tsv)){ 
      index<-which(temp_tsv[j,5] == final_tsv[,1])
      #print(index)
      if(identical(index, integer(0))){ 
        final_tsv[(nrow(final_tsv)+ 2), ]<- 0
        final_tsv[nrow(final_tsv),1]<-temp_tsv[j,5]
        final_tsv[nrow(final_tsv),2]<-as.character(temp_tsv[j,6])
        final_tsv[nrow(final_tsv),ncol(final_tsv)]<-temp_tsv[j,7]
        
        next
      }
      else{ 
        final_tsv[index,i+2]<-temp_tsv[j,7]
      }
      colnames(final_tsv)[i+2]<-file_name
    }
  }
  
  T1_RPM<-append(T1_RPM,taxa_detect(temp_tsv, 187217))
  MS2_RPM<-append(MS2_RPM,taxa_detect(temp_tsv, 329852))
}


final_tsv<-final_tsv[complete.cases(final_tsv),]

to_remove<-c()

for(i in 1:nrow(final_tsv)){ 
  if( all(final_tsv[i,3:ncol(final_tsv)] < 10 )){ 
    to_remove<-append(to_remove, i)
    }
  }

zero_removed<-final_tsv[-to_remove,]

wb = createWorkbook()

sheet = createSheet(wb, "RPM < 10")

addDataFrame(zero_removed, sheet=sheet, startColumn=1, row.names=FALSE)

sheet = createSheet(wb, "All RPM values")

addDataFrame(final_tsv, sheet=sheet, startColumn=1, row.names=FALSE)



saveWorkbook(wb, "RPM_summary.xlsx")







