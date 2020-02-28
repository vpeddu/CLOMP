# Read Kraken tsvs and give back excel file with RPM calculations 

if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
  }


setwd('/Users/gerbix/Documents/vikas/scratch/respo')

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
  temp_tsv<-read.csv(files[i], sep = '\t', col.names = c('percent_clade_reads', 'number_clade_reads_rooted_at_taxon','number_clade_reads_this_taxon', 'taxa', 'taxid', 'name'), header = FALSE)
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

# QC hardcoded to fail if RPM of T1 or MS2 < 100 
qc_df<-data.frame(files,T1_RPM,MS2_RPM)
qc_df$files<-as.character(qc_df$files)
qc_df$QC<-NA
for(i in 1:nrow(qc_df)){ 
  qc_df$files[i]<-strsplit(qc_df$files[i],"_final")[[1]][1]
  if(any(qc_df[i,c(2:3)] > 100 )){ 
    qc_df$QC[i]<-'PASS'
  }
  else { 
    qc_df$QC[i]<-'FAIL'
    }
  }

DNA_df<-final_tsv[,c(1,2,which(grepl('DNA',colnames(final_tsv))))]
RNA_df<-final_tsv[,c(1,2,which(grepl('RNA',colnames(final_tsv))))]

# Function to add in custom flags later
flag<-function(flagname,cutoff,df){ 
  df[,(ncol(df + 1))] <- NA
  colnames(df)[ncol(df) + 1]<-paste0(flagname)
  
  for(i in 1:nrow(df)){ 
    if(any(df[i,] > cutoff)){ 
      df[i,(ncol(df)+1)]<-'FLAG'
      }
  }
  return(df)
  }

# Function to flag columns where water control has RPM > cutoff value 
water_flag<-function(cutoff,df){ 
  df[,(ncol(df) + 1)] <- NA
  colnames(df)[ncol(df)]<-'WATER_PASS'
  water_cols<-which(grepl('*H2O*', colnames(df)))
  for(i in 1:nrow(df)){ 
    temp_names<-c()
    temp_list<-c()
    for(j in 3:(ncol(df) - 1)){
      #print(df[i,water_cols[1]])
    if(j == water_cols[1]){ 
      next
    }
    else { 
      if(df[i,water_cols[1]] < 1){ 
        next
      }
    #print(j)  
    temp_ratio<-df[i,j] / df[i,water_cols[1]]
    if(temp_ratio > cutoff){ 
      temp_list<-append(colnames(df)[j], temp_list)
      }
    #temp_list<-append(temp_ratio, temp_list) 
    #temp_names<-append(colname
      }
    }
    #print(temp_list)
    #if(any(temp_list < cutoff )){ 
      df[i,(ncol(df))]<-paste0((temp_list),collapse=', ' )
    #}
  }
  return(df)
}


DNA_df_flag<-water_flag(10,DNA_df)
RNA_df_flag<-water_flag(10,RNA_df)



wb = createWorkbook()

sheet = createSheet(wb, "QC")

addDataFrame(qc_df, sheet=sheet, startColumn=1, row.names=FALSE)

sheet = createSheet(wb, "DNA RPM values")

addDataFrame(DNA_df_flag, sheet=sheet, startColumn=1, row.names=FALSE)

sheet = createSheet(wb, "RNA RPM values")

addDataFrame(RNA_df_flag, sheet=sheet, startColumn=1, row.names=FALSE)

# Highlight failed QC rows
# fo <- Fill(foregroundColor="yellow")
# cs <- CellStyle(wb, fill=fo)
# sheets <- getSheets(wb)    
# sheet <- sheets[["QC"]]  
# 
# rows <- getRows(sheet, rowIndex=2:(nrow(qc_df)+1))
# cells <- getCells(rows, 4)
# values <- lapply(cells, getCellValue)
# 
# highlight <- "test"
# for (i in names(values)) {
#   x <- as.character(values[i])
#   print(x)
#   if (x == 'FAIL') {
#     highlight <- c(highlight, i)
#   }    
# }
# highlight <- highlight[-1]
# 
# lapply(names(cells[highlight]),
#        function(ii)setCellStyle(cells[[ii]],cs))

# Highlight water flagged rows in DNA sheet
# 
# fo <- Fill(foregroundColor="teal")
# cs <- CellStyle(wb, fill=fo)
# sheets <- getSheets(wb)    
# sheet <- sheets[["DNA RPM values"]]  
# 
# rows <- getRows(sheet, rowIndex=1:(nrow(DNA_df_flag)))
# cells <- getCells(rows, ncol(DNA_df_flag))
# values <- lapply(cells, getCellValue)
# 
# highlight <- "test"
# for (i in names(values)) {
#   x <- as.character(values[i])
#   if (x == 'FLAG') {
#   print(x)
#     highlight <- c(highlight, i)
#   }    
# }
# highlight <- highlight[-1]
# 
# lapply(names(cells[highlight]),
#        function(ii)setCellStyle(cells[[ii]],cs))


saveWorkbook(wb, "QC_data.xlsx")







