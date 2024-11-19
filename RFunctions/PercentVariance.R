percent.variance <- function(stand.dev, return.val = FALSE, plot.var = TRUE, plot.title = "Percent variance"){
	percent.variance <- stand.dev^2/sum(stand.dev^2) * 100
	if(plot.var == TRUE){
		plyr::round_any(max(percent.variance), 5, f = ceiling)
		barplot(percent.variance, 
						xlab = "PC",
						ylab = "Percent variance",
						las = 1,
						ylim = c(0, plyr::round_any(max(percent.variance), 3, f = ceiling)), 
						col = "blue", 
						main = plot.title
		)
	}
	if(return.val == TRUE){
		return(percent.variance)
	}
}


