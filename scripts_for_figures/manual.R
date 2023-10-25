

fig_df_s =  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_ds20.RData")

fig_df=  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens2 2.RData")

  #fig_df <- loadRData(input_file)
  
  fig_df_mut <- fig_df$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
  fig_df_mut_s <- fig_df_s$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
  
  
  
  
  
  
  K_ = max(fig_df_mut$K)
  n_vec = c( fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])
  alpha_val = fig_df_mut$alpha[1]
  
  
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  
  
  julia_assign("K_bound", K_)
  julia_assign("N1", n_vec[1])
  julia_assign("N2", n_vec[2])
  julia_assign("N3", n_vec[3])
  julia_assign("N4", n_vec[4])
  # alpha = e0*K
  julia_assign("alpha", fig_df$line$alpha[1]*K_)
  print( fig_df$line$alpha[1]*K_)
  #a = julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N1, N1, alpha)")
  
  df_prior = tibble(K= 1:K_, 
                    Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N1, K_bound, alpha)"),3),
                    Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N2, K_bound, alpha)"),3), 
                    Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N3, K_bound, alpha)"),3),
                    Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N4, K_bound, alpha)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
   
  #replace
  fig_df_mut$density[1:10] = fig_df_mut_s$density[1:10]
  fig_df_mut$Rh[1:10] = fig_df_mut_s$Rh[1:10]
  fig_df_mut$pkn[1:10] = fig_df_mut_s$pkn[1:10]
  
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`Pkn_1` = paste0("n = ", fig_df_mut$N[1]), `Pkn_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`Pkn_3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`Pkn_4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  p_h<- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
    ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\n =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$alpha[1],fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$=%3.f',fig_df_mut$N[1]),sprintf('$n$=%3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$=%3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$=%3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  p_h
  if (is.null(Sigma_coef)){
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"_2.pdf" ))
    plot(p_h)
    dev.off()
  }else{
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"R_coef_",Sigma_coef,"_2.pdf" ))
    plot(p_h)
    dev.off()
  }
  
  ### Sampe plot without title 
  p_h2 <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
    #ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$alpha[1],fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',fig_df_mut$N[1]),sprintf('$n$ = %3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  p_h2
  
  if (is.null(Sigma_coef)){
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"_2_.pdf" ))
    plot(p_h2)
    dev.off()
  }else{
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"R_coef_",Sigma_coef,"_2_.pdf" ))
    plot(p_h2)
    dev.off()
  }
  
  
  E_k<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
  write.table(E_k, file =paste0(fig_path,"Prior_Posterior_exp_alpha_",fig_df_mut$alpha[1],"_2.csv"), sep = ",", col.names = NA,
              qmethod = "double")
  
  
  
  
  ## Posterior of the weights
  
  pkn.labs<- c("W_1","W_2","W_3","W_4")
  names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ",n_vec[2]),paste0("n = ", n_vec[3]),paste0("n = ", n_vec[4]))
  pkn_names <- as_labeller(
    c(`W_1` = paste0("n = ", n_vec[1]), `W_2` = paste0("n = ", n_vec[2]),`W_3` = paste0("n = ", n_vec[3]),`W_4` = paste0("n = ", n_vec[4])))
  
  weights_fig <-fig_df$weights 
  weights_fig_s <- fig_df_s$weights 
  
  weights_fig$n <-as.factor(weights_fig$W_val) 
  weights_fig_s$n <-as.factor(weights_fig_s$W_val) 
  
  
  
  weights_fig_thin<- weights_fig%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)
  weights_fig_thin_s<- weights_fig_s%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)
  Weights = weights_fig_thin[order(weights_fig_thin$n,weights_fig_thin$it,weights_fig_thin$K),] 
  Weights[1:80000,] = weights_fig_thin_s[order(weights_fig_thin_s$it, weights_fig_thin_s$K),]
  
  pw <- ggplot(Weights, aes(x=K, y=weights, group =K, fill = n )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',n_vec[1]),sprintf('$n$ = %3.f',n_vec[2]),sprintf('$n$ = %3.f',n_vec[3]),sprintf('$n$ = %3.f',n_vec[4])))))+
    ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$alpha[1],n_vec[1],n_vec[2],n_vec[3],n_vec[4])))+
    theme_minimal()+  ylim(0,0.65)+facet_wrap(~W_,labeller = pkn_names)
  pw
  
  
  
  if (is.null(Sigma_coef)){
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"_3.pdf" ))
    plot(pw)
    dev.off()
  }else{
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"R_coef_",Sigma_coef,"_3.pdf" ))
    plot(pw)
    dev.off()
  }
  
  pw2 <- ggplot(Weights, aes(x=K, y=weights, group =K, fill = n )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',n_vec[1]),sprintf('$n$ = %3.f',n_vec[2]),sprintf('$n$ = %3.f',n_vec[3]),sprintf('$n$ = %3.f',n_vec[4])))))+
    #ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$alpha[1],n_vec[1],n_vec[2],n_vec[3],n_vec[4])))+
    theme_minimal()+  ylim(0,0.65)+facet_wrap(~W_,labeller = pkn_names)
  pw2
  
  
  if (is.null(Sigma_coef)){
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"_3_.pdf" ))
    plot(pw2)
    dev.off()
  }else{
    pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"R_coef_",Sigma_coef,"_3_.pdf" ))
    plot(pw2)
    dev.off()
  }
  
  ### MTM 
  if (MTM){
    df_k_20<- df_post(n_vec[1], c_vec=c_vec,ind=1, post = fig_df  )
    df_k_200<- df_post(n_vec[2], c_vec = c_vec,ind=2, post = fig_df  )
    df_k_2000<- df_post(n_vec[3], c_vec=c_vec,ind=3, post = fig_df  )
    df_k_20000<- df_post(n_vec[4],c_vec=c_vec,ind=4, post = fig_df  )
    df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)
    
    
    df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)
    pkn.labs<- n_vec
    names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ",  n_vec[2]),paste0("n = ",  n_vec[3]),paste0("n = ",  n_vec[4]))
    pkn_names <- as_labeller(
      c(`20` = paste0("n = ",  n_vec[1]), `200` = paste0("n = ",  n_vec[2]),`2000` = paste0("n = ",  n_vec[3]),`20000` = paste0("n = ",  n_vec[4])))
    
    
    df_fin_bar = df_fin %>%
      group_by(Process_type,N, pkn) %>%
      summarize(count=n())%>% mutate(pkn_dens =count/sum(count))
    K_ = max(fig_df_mut$K)
    alpha_ = fig_df$line$alpha[1]
    
    
    pm <- ggplot(df_fin_bar, aes(pkn,pkn_dens,color =Process_type, linetype = Process_type))+ geom_bar(aes(linetype=Process_type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
      geom_vline(xintercept=3,  linetype="dashed")+
      ggtitle(TeX(sprintf('PD for the num of clusters for MTM $\\alpha =%.3f$,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',fig_df_mut$alpha[1],c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
      scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
      theme_minimal()+ ylab("Density")+xlab(expression(tilde(K))) +
      scale_linetype_manual(values=c("solid", "dashed","longdash","dotted"),name = TeX(sprintf('$c$')),labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4]))))) +
      scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
      facet_wrap(~N,labeller = pkn_names)
    pm
    
    if (is.null(Sigma_coef)){
      pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"_4.pdf" ))
      plot(pm)
      dev.off()
    }else{
      pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"R_coef_",Sigma_coef,"_4.pdf" ))
      plot(pm)
      dev.off()
    }
    
    
    
    pm2 <- ggplot(df_fin_bar, aes(pkn,pkn_dens,color =Process_type, linetype = Process_type))+ geom_bar(aes(linetype=Process_type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
      geom_vline(xintercept=3,  linetype="dashed")+
      #   ggtitle(TeX(sprintf('PD for the num of clusters for MTM $\\alpha =%.3f$,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',fig_df_mut$alpha[1],c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
      scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
      theme_minimal()+ ylab("Density")+xlab(expression(tilde(K))) +
      scale_linetype_manual(values=c("solid", "dashed","longdash","dotted"),name = TeX(sprintf('$c$')),labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4]))))) +
      scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
      facet_wrap(~N,labeller = pkn_names)
    pm2
    
    if (is.null(Sigma_coef)){
      pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"_4_.pdf" ))
      plot(pm2)
      dev.off()
    }else{
      pdf(paste0(fig_path,"Figure1_alpha_",fig_df_mut$alpha[1],"R_coef_",Sigma_coef,"_4_.pdf" ))
      plot(pm2)
      dev.off()
    }
  }
  
}
