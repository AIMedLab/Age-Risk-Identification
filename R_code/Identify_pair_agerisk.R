library(dplyr)
library(reshape2)
library(stringr)
library(hash)
library(rcompanion)
library(speedglm)
library(lmtest)
library(compare)


#group by age
drug_outcome_demo = drug_outcome_demo[complete.cases(drug_outcome_demo[ , 'AGE']),]  #remove all rows with na in age
drug_outcome_demo = drug_outcome_demo[c('drug_concept_id','outcome_concept_id','ID','AGE')]

drug_outcome_demo$age_code = cut(abs(drug_outcome_demo$AGE),c(0,14,24,64,120))
drug_outcome_demo$age_code = as.character(drug_outcome_demo$age_code)
drug_outcome_demo$age_code[drug_outcome_demo$age_code == '(64,120]'] = 'senior'
drug_outcome_demo$age_code[drug_outcome_demo$age_code == '(24,64]'] = 'adult'
drug_outcome_demo$age_code[drug_outcome_demo$age_code == '(14,24]'] = 'youth'
drug_outcome_demo$age_code[drug_outcome_demo$age_code == '(0,14]'] = 'children'



#step1: overall chi-squared test for each drug

drug_outcome_demo['count'] = 1

adrpair_age = data.frame(drug_outcome_demo %>% group_by(drug_concept_id,outcome_concept_id,age_code) %>% summarize(count = sum(count))) # unique drug-ADR-group count
adrpair_age = adrpair_age[complete.cases(adrpair_age), ]

drug_overall_test = data.frame(unique(adrpair_age$drug_concept_id))  

drug_overall_test['X-squared'] = 0
drug_overall_test['df'] = 0
drug_overall_test['p_value'] = 0
drug_overall_test['p_value_adj'] = 0
drug_overall_test['overall.diff'] = 0
list=data.frame()
k=1

for (drug in unique(adrpair_age$drug_concept_id)){
  list = adrpair_age[adrpair_age$drug_concept_id==drug,]
  if (length(unique(list$age_code)) == 1){         
    drug_overall_test[k,'p_value']='no'
  }
  else{
    if(nrow(list)>1){
      contingency_table = dcast(list, age_code~outcome_concept_id)[-1]
      contingency_table[is.na(contingency_table)]=0
      matrix = matrix (chisq.test(contingency_table))
      drug_overall_test[k,'X-squared'] = matrix[1]
      drug_overall_test[k,'df'] = matrix[2]
      drug_overall_test[k,'p_value'] = matrix[3]
    }
  }
  k=k+1
}

drug_overall_test = drug_overall_test[which(drug_overall_test$p_value!='no'),] # remove drugs exist only in one group
drug_overall_test['p_value_adj']=p.adjust(drug_overall_test$p_value, method="bonferroni") 

drug_overall_test[which(drug_overall_test$p_value_adj <= 0.05),'overall.diff'] = 'TRUE'
drug_overall_test[which(drug_overall_test$p_value_adj > 0.05),'overall.diff'] = 'FALSE'


#step2: overall chi-squared test for drug-ADR pair

adrpair_age_2  = adrpair_age[which(adrpair_age$drug_concept_id %in% drug_overall_test[which(drug_overall_test$overall.diff=='TRUE'),'drug_concept_id']=='TRUE'),] 
adrpair_age_2 = adrpair_age_2 %>% group_by(drug_concept_id,outcome_concept_id)%>% summarise(total_count = sum(count)) # unique sig.drug-ADR count
adrpair_age_2 = data.frame(adrpair_age_2%>% filter(total_count>=50)) # filter pairs with occurrence less than 50

adrpair_age_2['X_squared'] = 0 
adrpair_age_2['df'] = 0  
adrpair_age_2['p_value'] = 0  
adrpair_age_2['p_value_adj'] = 0  
adrpair_age_2['overall.diff'] = 0 #col8
k=1


for (k in 1:nrow(adrpair_age_2)){
  list = adrpair_age[which(adrpair_age$drug_concept_id==adrpair_age_2[k,'drug_concept_id'] & adrpair_age$outcome_concept_id==adrpair_age_2[k,'outcome_concept_id']),]
  if (length(unique(list$age_code)) == 1){        
    adrpair_age_2[k,'p_value']='no'
  }
  else{
    if(nrow(list)>1){
      list_2 = adrpair_age[which(adrpair_age$drug_concept_id==adrpair_age_2[k,'drug_concept_id']),] 
      contingency_table = dcast(list_2, age_code~outcome_concept_id)[-1]
      contingency_table[is.na(contingency_table)]=0
      contingency_table['sum'] = rowSums (contingency_table[which(adrpair_age_2[k,'outcome_concept_id'] != colnames(contingency_table))]) 
      contingency_table = contingency_table[c(adrpair_age_2[k,'outcome_concept_id'],'sum')] 
      matrix = matrix (chisq.test(contingency_table))
      adrpair_age_2[k,'X_squared'] = matrix[1]
      adrpair_age_2[k,'df'] = matrix[2]
      adrpair_age_2[k,'p_value'] = matrix[3]
    }
  }
}

adrpair_age_2 = adrpair_age_2[which(adrpair_age_2$p_value!='no'),] # remove pairs exist only in one group
adrpair_age_2['p_value_adj']=adrpair_age_2$p_value * nrow(drug_outcome_demo)

adrpair_age_2[which(adrpair_age_2$p_value_adj <= 0.05),'overall.diff'] = 'TRUE'
adrpair_age_2[which(adrpair_age_2$p_value_adj > 0.05),'overall.diff'] = 'FALSE'


#step3: chi-squared tests for each drug-ADR pair

adrpair_sig = adrpair_age_2[which(adrpair_age_2$overall.diff=='TRUE'),] # all sig.pairs
adrpair_age_3 = adrpair_age[which(adrpair_age$drug_concept_id %in% drug_overall_test[which(drug_overall_test$overall.diff=='TRUE'),'drug_concept_id']=='TRUE'),] # all sig.drugs
adrpair_age_3 = data.frame(adrpair_age_3 %>% filter(count>=5))

pairwise_2 = data.frame(matrix(NA, nrow = 1, ncol = 12)) # pairs exist in two age groups
Colnames(pairwise_2) = c('drug_concept_id','outcome_concept_id','adult : children','adult : senior','adult : youth','children : senior','children : youth','senior : youth',
			'children','youth','adult','senior')
pairwise_3=pairwise_2 # pairs exist in three age groups
pairwise_4=pairwise_2 # pairs exist in four age groups

pairwise2_result = function(k, pairwise, contingency_table){
	  pairwise_2[nrow(pairwise_2)+1,] = -1
          pairwise_2[nrow(pairwise_2),'drug_concept_id'] = adrpair_sig[[k,'drug_concept_id']]
          pairwise_2[nrow(pairwise_2),'outcome_concept_id'] = adrpair_sig[[k,'outcome_concept_id']]
          
          col=pairwise[[1]]
          pairwise_2[nrow(pairwise_2),c(col)]=as.numeric(pairwise[1,3])
	  contingency_table_new = as.data.frame(t(contingency_table))
          names(contingency_table_new) = contingency_table_new[1,]
          for (group in names(contingency_table_new)){
            pairwise_2[nrow(pairwise_2),c(group)] = as.numeric(contingency_table_new[4,c(group)])
          }
        }


pairwise3_result = function(k, pairwise, contingency_table){
          pairwise_new = as.data.frame(t(pairwise))
          names(pairwise_new) = pairwise_new[1,]
          pairwise_new = pairwise_new[3,]
          
          pairwise_3[nrow(pairwise_3) + 1,]=-1
          pairwise_3[nrow(pairwise_3),'drug_concept_id'] = adrpair_sig[[k,'drug_concept_id']]
          pairwise_3[nrow(pairwise_3),'outcome_concept_id'] = adrpair_sig[[k,'outcome_concept_id']]          
          for (col in names(pairwise_new)){
            pairwise_3[nrow(pairwise_3),c(col)] = as.numeric(pairwise_new[1,c(col)])
          }
          
          contingency_table_new = as.data.frame(t(contingency_table))
          names(contingency_table_new) = contingency_table_new[1,]
          for (group in names(contingency_table_new)){
            pairwise_3[nrow(pairwise_3),c(group)] = as.numeric(contingency_table_new[4,c(group)])
          }
        }


pairwise41_result = function(k, pairwise, contingency_table){
	  pairwise_new = as.data.frame(t(pairwise))
          names(pairwise_new) = pairwise_new[1,]
          pairwise_new = pairwise_new[3,]
          
          pairwise_4[nrow(pairwise_4) + 1,]=-1
          pairwise_4[nrow(pairwise_4),'drug_concept_id'] = adrpair_sig[[k,'drug_concept_id']]
          pairwise_4[nrow(pairwise_4),'outcome_concept_id'] = adrpair_sig[[k,'outcome_concept_id']]   
          
          for (col in names(pairwise_new)){
            pairwise_4[nrow(pairwise_4),c(col)] = as.numeric(pairwise_new[1,c(col)])
          }
          contingency_table_new = as.data.frame(t(contingency_table))
          names(contingency_table_new) = contingency_table_new[1,]
          for (group in names(contingency_table_new)){
            pairwise_4[nrow(pairwise_4),c(group)] = as.numeric(contingency_table_new[4,c(group)])
          }
        }

pairwise42_result = function(k, pairwise, contingency_table){          
          pairwise_new = as.data.frame(t(pairwise))
          names(pairwise_new) = pairwise_new[1,]
          pairwise_new = pairwise_new[3,]
          
          sig_col = names(pairwise_new)[which(pairwise$p.adj.Chisq<=0.5)]
          sig_col=paste(sig_col[1],sig_col[2])
          if(length(unique(strsplit(sig_col, " ")[[1]]))==5){
            pairwise_4[nrow(pairwise_4) + 1,] = -1
            pairwise_4[nrow(pairwise_4),'drug_concept_id'] = adrpair_sig[[k,'drug_concept_id']]
            pairwise_4[nrow(pairwise_4),'outcome_concept_id'] = adrpair_sig[[k,'outcome_concept_id']]
            for (col in names(pairwise_new)){
              pairwise_4[nrow(pairwise_4),c(col)] = as.numeric(pairwise_new[1,c(col)])
            }

            contingency_table_new = as.data.frame(t(contingency_table))
            names(contingency_table_new) = contingency_table_new[1,]
            for (group in names(contingency_table_new)){
              pairwise_4[nrow(pairwise_4),c(group)] = as.numeric(contingency_table_new[4,c(group)])
            }            
          }
        }



for (k in 1:nrow(adrpair_sig)){
  list = adrpair_age_3[which(adrpair_age_3$drug_concept_id==adrpair_sig[[k,'drug_concept_id']] & adrpair_age_3$outcome_concept_id==adrpair_sig[[k,'outcome_concept_id']]),] # pair count
  all_list = adrpair_age_3[which(adrpair_age_3$drug_concept_id==adrpair_sig[[k,'drug_concept_id']]),] # drug count

  if(nrow(list)>1){
    adult_count = sum( all_list[which( all_list$age_code=='adult'),'count'])
    senior_count = sum( all_list[which( all_list$age_code=='senior'),'count'])
    children_count = sum( all_list[which( all_list$age_code=='children'),'count'])
    youth_count = sum( all_list[which( all_list$age_code=='youth'),'count'])
    contingency_table = dcast(list, age_code~outcome_concept_id)
    contingency_table=contingency_table[order(contingency_table$age_code),]
    contingency_table$other =0
    contingency_table$ratio=0
      
    for(i in 1:nrow(contingency_table)){
      if (contingency_table[i,1]=='adult'){contingency_table[i,3]=adult_count-contingency_table[i,2]
      contingency_table[i,4]=contingency_table[i,2]/contingency_table[i,3]}
      if (contingency_table[i,1]=='senior'){contingency_table[i,3]=senior_count-contingency_table[i,2]
      contingency_table[i,4]=contingency_table[i,2]/contingency_table[i,3]}
      if (contingency_table[i,1]=='youth'){contingency_table[i,3]=youth_count-contingency_table[i,2]
      contingency_table[i,4]=contingency_table[i,2]/contingency_table[i,3]}
      if (contingency_table[i,1]=='children'){contingency_table[i,3]=children_count-contingency_table[i,2]
      contingency_table[i,4]=contingency_table[i,2]/contingency_table[i,3]}
      }
      
    for(i in 0:nrow(contingency_table)){
      if(i==0){
        read_text = paste("\nOutcome" ,"target","others")}
      else{
        string = paste(contingency_table[i,1],contingency_table[i,2],contingency_table[i,3])
        read_text = paste(read_text,string,sep='\n')
        }
      }
    Matriz = as.matrix(read.table(textConnection(read_text),header=TRUE,row.names=1))
    pairwise = pairwiseNominalIndependence(Matriz,fisher = FALSE, gtest = FALSE, chisq = TRUE, method = "fdr")   #pairwise comparisons
      
    # put results in the tables
    if(nrow(list)==2){ # pairs in two groups
      if (pairwise$p.adj.Chisq<=0.05){
        pairwise2_result(k, pairwise, contingency_table)
          }        
      }         
     else if(nrow(list)==3){ # pairs in three groups
       if (length(which(pairwise$p.adj.Chisq<=0.05))!=1){
         pairwise3_result(k, pairwise, contingency_table)
        }
      }
     else if(nrow(list)==4){ # pairs in four groups
       if (length(which(pairwise$p.adj.Chisq<=0.05))==0|length(which(pairwise$p.adj.Chisq<=0.05))==6|length(which(pairwise$p.adj.Chisq<=0.05))==5){
         pairwise41_result(k, pairwise, contingency_table)
        }        
       else if (length(which(pairwise$p.adj.Chisq<=0.05))==4){
          pairwise42_result(k, pairwise, contingency_table)
        }
    }
   }
  }   


# format

pairwise_2$result=0
for (i in 1:nrow(pairwise_2)){
  ratio = pairwise_2[i,c('children','youth','adult','senior')]  
  max = colnames(ratio)[apply(ratio,1,which.max)] # get col name of max
  result = max
  pairwise_2$sig_group[i] = result
}

pairwise_3$result = 0
for (i in 1:nrow(pairwise_3)){
  p_val=pairwise_3[i,c('adult : children','adult : senior','adult : youth','children : senior','children : youth','senior : youth')]
  name = colnames(p_val)[which(p_val<=0.05 & p_val>=0)]
  ratio = pairwise_3[i,c('children','youth','adult','senior')]  
  max = colnames(ratio)[apply(ratio,1,which.max)] # get col name of max
  if(length(name)==2){
    name = gsub(':','',name)
    str = paste(name[1],name[2],sep=' ')
    list = as.list(strsplit(str, '\\s+')[[1]])    
    h=hash()
    for (j in 1:4){h[[list[[j]]]] = str_count(str,list[[j]])}
    two_class = invert(h)[['2']]   
    if(two_class==max){
      result=max
    }else{
      result=paste(invert(h)[['1']][1],invert(h)[['1']][2],sep='_')}
  }
    else if(length(name)==3){
      result = max
  }  #==0, no such cases  
  pairwise_3$sig_group[i] = result
}



pairwise_4$result=0   
for (i in 1:nrow(pairwise_4)){
  p_val = pairwise_4[i,c('adult : children','adult : senior','adult : youth','children : senior','children : youth','senior : youth')]
  name = colnames(p_val)[which(p_val>0.05)]  # equal cases
  ratio = pairwise_4[i,c('children','youth','adult','senior')]  
  max = colnames(ratio)[apply(ratio,1,which.max)] # get col name of max
  if(length(name)==1){
    name = gsub(':','',name)
    list = as.list(strsplit(name, '\\s+')[[1]])    
    if(max==list[[1]]|| max==list[[2]]){
      result=paste(list[[1]],list[[2]],sep='_')
    }else{result=max}
    }
    else if(length(name)==2){
      name = gsub(':','',name)
      list = as.list(strsplit(name, '\\s+'))
      if(max==list[[1]][1]||max==list[[1]][2]){
        result = paste(list[[1]][1],list[[1]][2],sep='_')
      }else{
        result = paste(list[[2]][1],list[[2]][2],sep='_')
      }
    }else{result = max}  
  pairwise_4$sig_group[i]=result
}

pairwise_regression = rbind(pairwise_2, pairwise_3)
pairwise_regression = rbind(pairwise_regression, pairwise_4)


#step4: logistic regression & lrt test

sig_pairs = pairwise_regression[,c('drug_concept_id','outcome_concept_id','sig_group')]
sig_pairs$lrt_p = -1
interaction_regr = drug_outcome_demo[c('drug_concept_id','outcome_concept_id','age_code')]

for (i in 1:nrow(sig_pairs)){
   drug_name = sig_pairs$drug_concept_id[i]
   condition_name = sig_pairs$outcome_concept_id[i]    
   interaction_regr$drug=0
   interaction_regr$condition=0   
   interaction_regr[which(interaction_regr$drug_concept_id== drug_name),c('drug')]=1
   interaction_regr[which(interaction_regr$outcome_concept_id== condition_name),c('condition')]=1     
   nested = speedglm(interaction_regr$condition ~ interaction_regr$drug+interaction_regr$age_code,family = binomial(link="logit"))  
   complex = speedglm(interaction_regr$condition ~ interaction_regr$drug+interaction_regr$age_code+interaction_regr$drug:interaction_regr$age_code,family = binomial(link="logit"))  
   result = lrtest(nested, complex)
   sig_pairs$lrt_p[i]=result$Pr[[2]]
  } 
}

sig_pairs['adj_p'] = sig_pairs$lrt_p * nrow(drug_outcome_demo)
sig_pairs = sig_pairs[which(sig_pairs$adj_p<=0.001),]









