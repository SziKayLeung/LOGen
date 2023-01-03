
suppressMessages(library(viridis)) 
suppressMessages(library(wesanderson)) 
suppressMessages(library(extrafont))
suppressMessages(loadfonts())

# plot aesthetics
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=18,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 18,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"),
                 plot.title = element_text(size=16))


# subplot aesthetics
mytheme_font <- theme(text=element_text(size=18,  family="CM Roman"),
                      axis.title.x = element_text(vjust=-0.5, colour = "black"),
                      axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                      legend.text = element_text(size = 12,family="CM Roman"),
                      axis.text.x= element_text(size=16,family="CM Roman"),
                      axis.text.y= element_text(size=16,family="CM Roman"))

# scale axis into 1000s
ks <- function(x){ format(x/1000, big.mark=",")} 

# change x axis
perc_lab <- function(x){ format(x* 100, big.mark=",")} 

