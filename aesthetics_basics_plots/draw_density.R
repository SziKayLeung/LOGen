
suppressMessages(library(grid))
# density plot
density_plot <- function(dat,x.var,y.var, x_lab, y_lab,title){
  
  print(paste0(title))
  print(paste0("Correlation between", x.var, "and", y.var))
  
  print(cor.test(dat[[x.var]],dat[[y.var]]))
  cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  
  corr.value <- cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  p.value <- cor.test(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")$p.value 
  
  
  corr <- grobTree(textGrob(paste("r = ", round(corr.value, 2)), 
                            x = 0.05, y = 0.80, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 14, fontface = "italic",family="CM Roman")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  print(paste0("corr.value", corr.value))
  print(paste0("p.value", p.value))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1, name = "Density") +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = paste(title,"")) + 
    geom_smooth(method=lm, colour = "black") + 
    mytheme + 
    theme(legend.position = "none")
  
  return(p)
}


corr_plot <- function(dat,x.var,y.var, x_lab, y_lab,title){
  
  print(paste0(title))
  print(paste0("Correlation between", x.var, "and", y.var))
  
  print(cor.test(dat[[x.var]],dat[[y.var]]))
  cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  
  corr.value <- cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  p.value <- cor.test(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")$p.value 
  
  corr <- grobTree(textGrob(paste0("r = ", round(corr.value, 2), "; P = ", round(p.value,2)), 
                            x = 0.05, y = 0.90, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 14, fontface = "italic",family="CM Roman")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  print(paste0("corr.value", corr.value))
  print(paste0("p.value", p.value))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    geom_point(size = 3) +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = paste(title,"")) + 
    geom_smooth(method=lm, colour = "black", se = FALSE) + 
    mytheme + 
    theme(legend.position = "none")
  
  return(p)
}


# Aim: plot the histogram of the feature (x.var), with Dataset as colour identifier 
# Input:
  # dat: df to plot, ensure "Dataset" column exists
  # x.var = str: the column of df for plotting as x-variable
  # xlab = str: x-label axis
# Output: 
  # histogram of feature distribution, colour-coded by dataset

geom_density_plot <- function(dat, x.var, xlab){
  
  p <- ggplot(dat, aes(x = !! rlang::sym(quo_name(enquo(x.var))), fill = Dataset)) + 
    geom_density(alpha = 0.2) + mytheme + labs(x = xlab, y = "Density") 
  
  return(p)
}
