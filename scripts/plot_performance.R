#!/usr/bin/env Rscript

##
#  Index construction time plot
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  The input is a CSV file with this scheme:
#    d2,method,const-time
##

library('viridis')
library('tibble')
library('tidyverse')
library('hrbrthemes')
library('ggplot2')
library('optparse')
library('egg')

## Default values
in_file <- 'benchmark.csv'
p_exclude_space <- 'cuda'
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-e', '--exclude'), type='character', default=p_exclude_space, help='Space value to exclude for index size and query time plots\n\t\t[default=\'%default\']', metavar='character'))
opt_parser <- OptionParser(option_list=option_list)
## Parsing command-line arguments
opt <- parse_args(opt_parser)

# Fix Roboto Condense Font ##############################
hrbrthemes::import_roboto_condensed()
#b <- read.csv(extrafont:::fonttable_file(), stringsAsFactors = FALSE)
#b[grepl("Light", b$FontName),]$FamilyName <- font_rc_light
#write.csv(b,extrafont:::fonttable_file(), row.names = FALSE)
extrafont::loadfonts()
#########################################################

## Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)

title <- 'Index Construction Time'
g1 <- ggplot(d, aes(y=const_time, x=d2, color=method, shape=method, linetype=space)) +
  geom_line(linewidth=1) +
  geom_point(size=4) +
  labs(title=title, y="Construction time (s)", x=expression("d"["2"]),
       color="Method", shape="Method", linetype="Execution Space") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(breaks=c(128, 256, 512, 1024), trans='log2') +
  scale_shape_manual(values=c(20, 18)) +
#  scale_color_viridis(discrete=TRUE) +
#  coord_cartesian(xlim=c(110, 1100)) +
  coord_fixed(ratio=1) +
  theme_ipsum_rc(base_size=16, plot_title_size=20, axis_title_size=17,
                 axis_title_just = "cc" )

ggsave('construction_time.pdf', device=cairo_pdf)

# Drop `const_time` column
ds <- d[d$space == opt$exclude, -which(names(d) %in% c("const_time"))]

title <- 'Query Time'
g2 <- ggplot(ds, aes(y=query_time, x=d2, color=method, shape=method)) +
  geom_line(linewidth=1) +
  geom_point(size=4) +
  labs(title=title, y="Query time (ns)", x=expression("d"["2"]),
       color="Method", shape="Method") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(breaks=c(128, 256, 512, 1024), trans='log2') +
  scale_shape_manual(values=c(20, 18)) +
#  scale_color_viridis(discrete=TRUE) +
#  coord_cartesian(xlim=c(110, 1100)) +
#  coord_fixed(ratio=2) +
  theme_ipsum_rc(base_size=16, plot_title_size=20, axis_title_size=17,
                 axis_title_just = "cc" )

ggsave('query_time.pdf', device=cairo_pdf)

title <- 'Index Size'
g3 <- ggplot(ds, aes(y=index_size, x=d2, color=method, shape=method)) +
  geom_line(linewidth=1) +
  geom_point(size=4) +
  labs(title=title, y="Index size (MB)", x=expression("d"["2"]),
       color="Method", shape="Method") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(breaks=c(128, 256, 512, 1024), trans='log2') +
  scale_shape_manual(values=c(20, 18)) +
#  scale_color_viridis(discrete=TRUE, direction=-1) +
#  coord_cartesian(xlim=c(110, 1100)) +
  coord_fixed(ratio=1) +
  theme_ipsum_rc(base_size=16, plot_title_size=20, axis_title_size=17,
                 axis_title_just = "cc" )

ggsave('index_size.pdf', device=cairo_pdf)

g2 <- g2 + guides(shape='none', color='none', linetype='none')
g3 <- g3 + guides(linetype='none')
gs <- ggarrange(g2, g3, g1, ncol=3, nrow=1)
ggsave('all.pdf', gs, device=cairo_pdf, width=20, height=6)
