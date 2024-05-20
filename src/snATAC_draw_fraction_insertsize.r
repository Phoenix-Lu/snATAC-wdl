library(tidyverse)

Args <- commandArgs(trailingOnly = TRUE)

Args <- c("filter_all_filtered_sorted.insertlen")
data <- read.table(Args[1], header = F, sep = "\t")
data <- data[2]
# 去除含零行
data <- tibble::as_tibble(data) %>% filter(V2!=0)
# 设置插入片段长度的阈值，过滤掉太长的片段


length_cutoff <- 1200
fragment <- data$V2[data$V2 <= length_cutoff]


######
##Part1：基础语法画图
######
# 利用直方图统计频数分布，设置柱子个数
breaks_num <- 500
res <- hist(fragment, breaks = breaks_num, plot = FALSE)


# # 添加坐标原点
# plot(x = c(0, res$breaks),
#      y = c(0, 0, res$counts) / 10^2,
#      type = "l", col = "red",
#      xlab = "Fragment length(bp)",
#      ylab = "Normalized read density 10^2",
#      main = paste("Sample Fragment sizes\n",
#        basename(Args[1])))

######
##Part2：ggplot2 画图及其拼图
######
## 不同数据分布
DATA <- data.frame(x1 = c(0, res$breaks),y1=c(0, 0, res$counts) / 10^2)
p1 <- ggplot(DATA,aes(x =x1,y = y1 )) +
  geom_line(col="red") +
  xlab("Fragment length(bp)") +
  ylab("Normalized read density 10^2") +
  ggtitle("Sample Fragment sizes") +
  theme_classic() 
  
## 画小图
DATA2 <- data.frame(x1 = c(0, res$breaks),y1=log10(c(0, 0, res$counts) / 10^2)+1)
p2 <- ggplot(DATA2,aes(x =x1,y = y1 ))+
  geom_line(col="red")+
  xlab("Fragment length(bp)")+
  ylab(expression("Normalized read count (lg)"))+
  ggtitle("Sample Fragment sizes")+
  theme_classic()


## 小图插入右上角
library(cowplot)
p3 <- ggdraw() +
  draw_plot(p1, 0, 0, 1, 1) +
  draw_plot(p2, x=0.5,y=0.5, width = 0.5, height = 0.5) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), size = 15) +
  draw_figure_label(basename(Args[1]),
                    position = 'bottom.left',
                    size = 15,
                    colour = 'red',
                    fontface = "bold") 

ggsave(paste(basename(Args[1]), ".png", sep=''), p3, width = 10, height = 8)
