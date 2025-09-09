ind1 <- c(1, 3, 5, 2, 4, 6)
ind2 <- c(c(1, 4, 7), c(1, 4, 7) + 7, c(1, 4, 7) + 14, c(1, 4, 7) + 21)
output1 <- t(read.table("output/output_S1.txt")[ind1, ind2])
output2 <- t(read.table("output/output_S2.txt")[ind1, ind2])
output3 <- t(read.table("output/output_S3.txt")[ind1, ind2])
output4 <- t(read.table("output/output_S4.txt")[ind1, ind2])

plot.fdr.fun <- function(output){
  l <- nrow(output) /  2
  means <- output[1 : l, ]
  errors <- output[(l + 1) : (2 * l), ]
  rownames(means) <- c("LinDA", "LinDA-BMDD", "LinDA-SAVER")
  colnames(means) <- c("m=50;n=50", "m=200;n=50", "m=500;n=50",
                       "m=50;n=200", "m=200;n=200", "m=500;n=200")
  bar_positions <- barplot(means, beside = TRUE,  
                           col = c("lightgrey", "lightcoral", "lightblue"),
                           xlab = 'Empirical False Discovery Rate',
                           horiz = TRUE,
                           xlim = c(0, max(means + errors) + 0.01),
                           yaxt = "n",
                           cex.axis = 1.2, cex.lab = 1.2)
  
  for (i in 1:nrow(means)) {
    arrows(means[i, ] - 1.96*errors[i, ], 
           bar_positions[i, ], 
           means[i, ] + 1.96*errors[i, ], 
           bar_positions[i, ], 
           angle = 90, 
           code = 3, 
           length = 0.1)
  }
  
  abline(v = 0.05, col = "black", lty = 2, lwd = 1)
  axis(2, at = bar_positions[2, ], labels = FALSE, 
       las = 1, cex.axis = 1.2, tick = FALSE)
  if(setup == 'S1') {
    a <- 0.008
  } else if(setup == 'S2') {
    a <- 0.016 
  } else if(setup == 'S3') {
    a <- 0.055
  } else if(setup == 'S4') {
    a <- 0.032
  }
  text(par("usr")[1]-a, bar_positions[2, ], colnames(means), srt = 60, xpd = TRUE, cex = 1.2)
}


plot.power.fun <- function(output){
  l <- nrow(output) /  2
  means <- output[1 : l, ]
  errors <- output[(l + 1) : (2 * l), ]
  rownames(means) <- c("LinDA", "LinDA-BMDD", "LinDA-SAVER")
  colnames(means) <- c("m=50;n=50", "m=200;n=50", "m=500;n=50",
                       "m=50;n=200", "m=200;n=200", "m=500;n=200")
  bar_positions <- barplot(means, beside = TRUE,  
                           col = c("lightgrey", "lightcoral", "lightblue"),
                           xlab = 'True Positive Rate',
                           horiz = TRUE,
                           xlim = c(0, max(means + errors) + 0.1),
                           yaxt = "n",
                           cex.axis = 1.2, cex.lab = 1.2)
  
  for (i in 1:nrow(means)) {
    arrows(means[i, ] - 1.96*errors[i, ], 
           bar_positions[i, ], 
           means[i, ] + 1.96*errors[i, ], 
           bar_positions[i, ], 
           angle = 90, 
           code = 3, 
           length = 0.1)
  }
  
  legend("bottomright", legend = rownames(means), 
         fill = c("lightgrey", "lightcoral", "lightblue"), 
         cex = 1.2,  # Adjust text size
         bty = "n",  # No box around the legend
         inset = c(-0.02, 0.02)) # Horizontal, vertical movement
}

combine.fun <- function(setup, output) {
  pdf(paste0(setup, '_barplot.pdf'), width = 11, height = 8)
  par(mfrow = c(1, 2))
  par(mar = c(5, 4, 4-2, 2-1) + 0.1)
  plot.fdr.fun(output[c(1 : 3 , 7 : 9), ])
  par(mar = c(5, 4-3, 4-2, 2) + 0.1)
  plot.power.fun(output[c(4 : 6 , 10 : 12), ])
  dev.off()
}

setup <- 'S1'
output <- output1
combine.fun(setup, output)

setup <- 'S2'
output <- output2
combine.fun(setup, output)

setup <- 'S3'
output <- output3
combine.fun(setup, output)

setup <- 'S4'
output <- output4
combine.fun(setup, output)




########################################
plot.power.fun <- function(output, s){
  l <- nrow(output) /  2
  means <- output[1 : l, ]
  errors <- output[(l + 1) : (2 * l), ]
  rownames(means) <- c("LinDA", "LinDA-BMDD", "LinDA-SAVER")
  colnames(means) <- c("m=50;n=50", "m=200;n=50", "m=500;n=50",
                       "m=50;n=200", "m=200;n=200", "m=500;n=200")
  bar_positions <- barplot(means, beside = TRUE,  
                           col = c("lightgrey", "lightcoral", "lightblue"),
                           xlab = 'True Positive Rate',
                           horiz = TRUE,
                           xlim = c(0, max(means + errors) + 0.1),
                           yaxt = "n",
                           cex.axis = 1.2, cex.lab = 1.2)
  
  for (i in 1:nrow(means)) {
    arrows(means[i, ] - 1.96*errors[i, ], 
           bar_positions[i, ], 
           means[i, ] + 1.96*errors[i, ], 
           bar_positions[i, ], 
           angle = 90, 
           code = 3, 
           length = 0.1)
  }
  
  if (s == 4) {
    legend("right", legend = rownames(means), 
           fill = c("lightgrey", "lightcoral", "lightblue"), 
           cex = 1.2,  # Adjust text size
           bty = "n",  # No box around the legend
           inset = c(-0.4, -0.4), # Horizontal, vertical movement
           xpd = TRUE) 
  }

}

combine.fun <- function(setup, output, s) {
  if(s == 4) {
    pdf(paste0(setup, '_barplot_2.pdf'), width = 13, height = 8)
  } else {
    pdf(paste0(setup, '_barplot_2.pdf'), width = 11, height = 8)
  }
  par(mfrow = c(1, 2))
  par(mar = c(5, 4, 4-2, 2-1) + 0.1)
  plot.fdr.fun(output[c(1 : 3 , 7 : 9), ])
  if(s == 4) {
    par(mar = c(5, 4-3, 4-2, 10) + 0.1)
  } else {
    par(mar = c(5, 4-3, 4-2, 2) + 0.1)
  }
  plot.power.fun(output[c(4 : 6 , 10 : 12), ], s)
  dev.off()
}


setup <- 'S2'
output <- output2
combine.fun(setup, output, 2)

setup <- 'S3'
output <- output3
combine.fun(setup, output, 3)

setup <- 'S4'
output <- output4
combine.fun(setup, output, 4)




