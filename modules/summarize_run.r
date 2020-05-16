        #Read Kraken tsvs and give back excel file with RPM calculations 
        library('xlsx')
        library('data.table')
        library('tidyr')



        files<-list.files(path = ".", pattern = '*.tsv')
        print(files)
        taxa_detect<-function(df, taxid){ 
        temp_rpm<-df$RPM[which(df$taxid == taxid)]
        if(identical(temp_rpm, numeric(0))){ 
          temp_rpm <- 0
        }
        return(temp_rpm)
        }



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
      }


      final_tsv<-final_tsv[complete.cases(final_tsv),]
      cname<-c()
      for(i in 3:length(colnames(final_tsv))){ 
        cname<-append(cname, strsplit(colnames(final_tsv)[i], '_')[[1]][1])
        }
      colnames(final_tsv)[3:length(colnames(final_tsv))]<-cname

      to_remove<-c()

      for(i in 1:nrow(final_tsv)){ 
        if( all(final_tsv[i,3:ncol(final_tsv)] < 10 )){ 
          to_remove<-append(to_remove, i)
        }
      }

      if( length(to_remove) > 0 ){
      zero_removed<-final_tsv[-to_remove,]
      } else{ 
          zero_removed<-final_tsv
      }
  write.csv(final_tsv, 'RPM_summary.csv')