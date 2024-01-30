
suppressMessages(library(viridis)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(wesanderson)) 
suppressMessages(library(extrafont))
suppressMessages(loadfonts())

# negate
`%ni%` <- Negate(`%in%`)

# plot aesthetics
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=16),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 12),
                 axis.text.x= element_text(size=16),
                 axis.text.y= element_text(size=16),
                 plot.title = element_text(size=16),
                 plot.subtitle = element_text(size=14))

mythemeNoLegend <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=16),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.text = element_text(size = 16),
                 axis.text.x= element_text(size=16),
                 axis.text.y= element_text(size=16),
                 plot.title = element_text(size=16),
                 plot.subtitle = element_text(size=14))



# subplot aesthetics
mytheme_font <- theme(text=element_text(size=16),
                      axis.title.x = element_text(vjust=-0.5, colour = "black"),
                      axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                      legend.text = element_text(size = 12),
                      axis.text.x= element_text(size=16),
                      axis.text.y= element_text(size=16))

# scale axis into 1000s 
# e.g. scale_y_continuous(labels = ks) 
ks <- function(x){ format(x/1000, big.mark=",")} 

# change x axis
perc_lab <- function(x){ format(x* 100, big.mark=",")} 

# ticks for x axis
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# use for plotting (scale_x_discrete) to determine how many breaks to be used
# e.g scale_x_discrete(breaks = every_nth(n = 2)) 
plot_break_every_nth <- function(n) {return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})}
