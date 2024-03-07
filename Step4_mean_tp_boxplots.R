## mean transition probabilities ##

act_id_1 <- sc_5@active.ident
act_id_2 <- sc_6@active.ident

rv_1 <- nsc1_tp_5  #Step 3 result for control
rv_2 <- nsc1_tp_6  #Step 3 result for TBI

exp_1 <- "Control"
exp_2 <- "TBI"

int_cl_name = "A-stage1" #interested cluster name
neighbours = c("RG-like","A-stage1","A-stage2")


tp_boxplot(int_cl_name,neighbours,act_id_1,act_id_2,
           exp_1,exp_2,rv_1,rv_2)



tp_boxplot <- function(int_cl_name,neighbours,act_id_1,act_id_2,
                       exp_1,exp_2,rv_1,rv_2) {
  
  #subseting cluster cells form columns
  cluster1 <- which(act_id_1 == int_cl_name)
  cluster2 <- which(act_id_2 == int_cl_name)
  
  
  #subsetting neighbour cells from rows
  neighbours_1 <- grep(pattern = paste(neighbours,collapse = "|"),x = act_id_1)
  neighbours_2 <- grep(pattern = paste(neighbours,collapse = "|"),x = act_id_2)
  
  
  #transition probabilities 
  tp_1 <- as.matrix(rv_1$tp[,cluster1])
  tp_2 <- as.matrix(rv_2$tp[,cluster2])
  
  
  tp_1_df <- data.frame(cluster = act_id_1,tp_1)%>%group_by(cluster)%>%summarise_all(sum)
  tp_2_df <- data.frame(cluster = act_id_2,tp_2)%>%group_by(cluster)%>%summarise_all(sum)
  
  
  
  tp_df <- rbind(reshape2::melt(cbind(sample = rep(exp_1,times = length(neighbours)),tp_1_df)),
                 reshape2::melt(cbind(sample = rep(exp_2,times = length(neighbours)),tp_2_df)))
  
  
  ggboxplot(tp_df, x = "cluster", y = "value",ylab = "Transition Probability",xlab = "Cluster",bxp.errorbar = T,
            fill = "sample", palette = c("#F8766D","#00BFC4"), title = paste(int_cl_name,"Mean Transition Probability")) +
    stat_compare_means(aes(group = sample),label = "p.signif",method = "wilcox.test",size = 6,label.y.npc = c("top"),vjust = 0.4) +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30,vjust = 0.5),
          axis.title.x = element_blank())
