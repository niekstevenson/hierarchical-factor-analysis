library(corrplot)
library(colorspace)

plot_relations <- function(sampled = NULL, loadings = NULL, corrs = NULL, cred = TRUE, do_corr = F){
  if(do_corr || !is.null(corrs)){
    if(is.null(corrs)){
      values <- sampled$samples$theta_var[,,sampled$samples$stage == "sample", drop = F]
    } else{
      values <- corrs
    }
    means <- cov2cor(apply(values, 1:2, mean))
    corrplot(means, addCoef.col ="black", number.cex = .75, tl.col = "black")
  } else{
    # Now we assume loadings
    cols <- diverging_hcl(200, palette = "Red-Green")
    if(is.null(loadings)){
      values <- sampled$samples$theta_lambda[,,sampled$samples$stage == "sample", drop = F]
    } else{
      values <- loadings
    }
    means <- apply(values, 1:2, mean)
    corrplot(means, is.corr = F, col = cols, col.lim = c(-1.2, 1.2), cl.pos = "n", addCoef.col ="black", number.cex = .75, tl.col = "black")
    # You might have to play around with the x and y limits of the legend.
    colorlegend(xlim=c(2 + .5*ncol(values),4 + .5*ncol(values)), ylim=c(2,6), cols, c(seq(-1,1,.25)), align="l", vertical=TRUE, addlabels=TRUE) 
  }
  
  
  # You might have to play around with the x and y limits of the legend.
  
  # Only add this if you want confidence intervals
  if(cred){
    if(do_corr){
      for(i in 1:dim(values)[3]){
        values[,,i] <- cov2cor(values[,,i])
      }
    }
    
    cred <- aperm(apply(values, 1:2, quantile, probs = c(0.025, 0.975)))
    cred <- round(cred, 2)
    print(mean(abs(cred[,,1] - cred[,,2])))
    conf <- paste0("[", format(cred[,,1,drop = F], digits=1), ":", format(cred[,,2,drop = F], digits=1), "]")
    
    
    xs <- row(cred[,,1])
    ys <- (ncol(cred[,,1])+1) - col(cred[,,2]) - 0.15
    text(xs, ys, conf, pos=1, cex=0.55, font = 2)
  }
  
}


# Make some nice GGplot, could honestly be done with a lot less work with the bayesplot package
lambda_recovery_plot <- function(sampled, lambda_input, stage = "sample"){
  lambda_recovered <- sampled$samples$theta_lambda[,,sampled$samples$stage == stage]
  values <- as.numeric(lambda_recovered[,,1])
  for(x in 2:dim(lambda_recovered)[3]){
    lambda_tmp <- lambda_recovered[,,x]
    values <- c(values, as.numeric(lambda_tmp))
  }
  # First get the parameter names
  parNames <- expand.grid(rownames(lambda_recovered), 1:n_factors)
  loadings <- rep(parNames[,1], dim(lambda_recovered)[3])
  factors <- rep(parNames[,2], dim(lambda_recovered)[3])
  iteration <- rep(1:dim(lambda_recovered)[3], each = nrow(lambda_recovered)*n_factors)
  all_data <- data.frame(iteration = iteration, loadings = loadings, factors = factors, value = values)
  # We use the means to color the distributions
  means <- aggregate(value ~ loadings*factors, data = all_data, mean)
  all_data$means <- rep(means$value, dim(lambda_recovered)[3])
  
  true_lambda <- as.numeric(lambda_input)
  true_data <- data.frame(iteration = rep(-1, length(true_lambda)), loadings = parNames[,1], factors = parNames[,2], value = true_lambda, means = means$value)
  all_data <- rbind(all_data, true_data)
  
  #Remove constraint parameters
  all_data <- all_data[all_data$means != 0,]
  
  library(ggplot2)
  
  library(colorspace)
  cols <- diverging_hcl(200, palette = "Red-Green")
  cols <- cols[-c(80:120)]
  # The ggplot itself
  tmp <- ggplot(all_data, aes(value)) +   facet_grid(loadings ~ factors, switch = "y") + 
    geom_density(data = subset(all_data, iteration != -1, ), aes(colour=means, fill=means))+
    geom_vline(data = subset(all_data, iteration == -1), aes(xintercept = value), linewidth = 1.1, lty = "11") + 
    xlim(c(-1.1,1.1)) + 
    theme_bw() +
    scale_fill_gradientn(colours=cols, limits = c(-1, 1)) + 
    scale_colour_gradientn(colours=cols, limits = c(-1, 1)) +
    xlab("loadings") + 
    ylab("parameters") +
    theme(
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5),
      panel.spacing = unit(0, "lines"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
      axis.text.y = element_blank(),
      axis.ticks.y=element_blank(),
      axis.title=element_text(size=16),
      strip.text = element_text(size = 12)) 
  plot(tmp)
  
}

