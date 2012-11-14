PlotResults <-
function(r, ActualValues, LocalCI=NA, GlobalCI=NA, xlab="r", ylab="", ReferenceValue=NA, Legend=FALSE, LegendItems=c("Value", "Local CI", "Global CI"), LegendPosition="topright") {
	if (is.na(LocalCI)[1]) {
		AllValues <- unlist(c(ActualValues, GlobalCI))
	} else {
		AllValues <- unlist(c(ActualValues, LocalCI$Min, LocalCI$Max, GlobalCI))
	}
	LegendColors <- c("red", "black", "black")
	LegendLineTypes <- c(1, 3, 4)
	plot(r, ActualValues, type="l", col="red", xlab=xlab, ylab=ylab, ylim=c(min(AllValues, na.rm = TRUE), max(AllValues, na.rm = TRUE)))
	LegendToShow <- 1
	if (!is.na(LocalCI)[1]) {
		lines(r, LocalCI$Min, lty=3)
		lines(r, LocalCI$Max, lty=3)
		LegendToShow <- c(LegendToShow, 2)
	} 
	if (!is.na(GlobalCI)[1]) {
		lines(r, GlobalCI$Min, lty=4)
		lines(r, GlobalCI$Max, lty=4)
		LegendToShow <- c(LegendToShow, 3)
	}
	abline(h=ReferenceValue, lty=2)
	if (Legend) {
		LegendItems <- LegendItems[1:length(LegendToShow)]
		LegendColors <- LegendColors[LegendToShow]
		LegendLineTypes <- LegendLineTypes[LegendToShow]
		legend(LegendPosition, legend = LegendItems, col = LegendColors, lty = LegendLineTypes, inset=0.02)
	}
}
