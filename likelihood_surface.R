source("functions.R")
library(ggplot2)
library(directlabels)
library(cowplot)
# Generate combinations of alt allele and depths
alt <- seq(1,10,1) 
depth <- seq(100,1000,50)

fr <- crossing(alt,depth)

fr_wgs <- fr %>% mutate(lik= map2_dbl(.$depth,.$alt, model_likelihood),
                    null_lik = map2_dbl(.$depth,.$alt,null_likelihood),
                    log_expected_errors = map2_dbl(.$alt,.$depth,error_num_expected_wgs),
                    tlod = log10(lik) - log10(null_lik),
                    vaf = alt/depth,
                    log_expected_true = truth_expected_wgs(vaf))



p2 <- ggplot(fr_wgs,aes(x=depth,y=alt)) + geom_contour(aes(z=vaf,colour=stat(level)),binwidth=0.005)  + scale_y_continuous(breaks=seq(1,10,1))#+ geom_contour(aes(z=tlod,colour=stat(level)))
p3 <- direct.label(p2,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("Vaf contours")

p1 <- ggplot(fr_wgs,aes(x=depth,y=alt,z=tlod)) + geom_contour(aes(colour=stat(level)),binwidth=1) + scale_y_continuous(breaks=seq(1,10,1))
p4 <- direct.label(p1,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("TLOD contours")

p5 <- ggplot(fr_wgs,aes(x=depth,y=alt,z=log_expected_errors)) + geom_contour(aes(colour=stat(level))) + scale_y_continuous(breaks=seq(1,10,1))
p6 <- direct.label(p5,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("Log expected false positives")

p7 <- ggplot(fr_wgs,aes(x=depth,y=alt,z=log_expected_true)) + geom_contour(aes(colour=stat(level)),binwidth=1) + scale_y_continuous(breaks=seq(1,10,1))
p8 <- direct.label(p7,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("Log expected true positives")


p_final_wgs <- cowplot::plot_grid(p3,p4,p6,p8,nrow=2)


# WES
fr_wes <- fr %>% mutate(lik= map2_dbl(.$depth,.$alt, model_likelihood),
                        null_lik = map2_dbl(.$depth,.$alt,null_likelihood),
                        log_expected_errors = map2_dbl(.$alt,.$depth,error_num_expected_wes),
                        tlod = log10(lik) - log10(null_lik),
                        vaf = alt/depth,
                        log_expected_true = truth_expected_wes(vaf),
                        tf_ratio = 10**(log_expected_true-log_expected_errors))



p9 <- ggplot(fr_wes,aes(x=depth,y=alt)) + geom_contour(aes(z=vaf,colour=stat(level)),binwidth=0.005)  + scale_y_continuous(breaks=seq(1,10,1))#+ geom_contour(aes(z=tlod,colour=stat(level)))
p10 <- direct.label(p9,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("Vaf contours")

p11 <- ggplot(fr_wes,aes(x=depth,y=alt,z=tlod)) + geom_contour(aes(colour=stat(level)),binwidth=1) + scale_y_continuous(breaks=seq(1,10,1))
p12 <- direct.label(p11,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("TLOD contours")

p13 <- ggplot(fr_wes,aes(x=depth,y=alt,z=log_expected_errors)) + geom_contour(aes(colour=stat(level))) + scale_y_continuous(breaks=seq(1,10,1))
p14 <- direct.label(p13,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("Log expected false positives")

p15 <- ggplot(fr_wes,aes(x=depth,y=alt,z=log_expected_true)) + geom_contour(aes(colour=stat(level)),binwidth=1) + scale_y_continuous(breaks=seq(1,10,1))
p16 <- direct.label(p15,list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust=1,vjust=1,
                           box.color = NA, fill = "transparent", "draw.rects")) + 
  theme_bw() + ggtitle("Log expected true positives")

p17 <- ggplot(fr_wes %>% dplyr::filter(tlod < 8 & tlod > 4),aes(x=tlod,y=alt/depth)) + geom_point(aes(colour=as.factor(tf_ratio)))
p18 <- ggplot(fr_wes %>% dplyr::filter(tlod < 8 & tlod > 4),aes(x=depth,y=alt/depth,z=tf_ratio)) + geom_contour(aes(colour=stat(level)))

p_final_wes <- cowplot::plot_grid(p10,p12,p14,p16,nrow=2)