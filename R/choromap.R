choromap <- function(sf, values, breaks, n = 9,
                     name = "PuBu", divcol = FALSE,
                     placelegend = "right",
                     title = "", space = 8,
                     showlegend = TRUE, ...){

  class <- classInt::classIntervals(values, style = "fixed", fixedBreaks = breaks)
  pal <- RColorBrewer::brewer.pal(n, name)
  if (divcol) {pal <- rev(pal)}
  col.pal <- colorRampPalette(pal)((length(breaks)-1))
  which.col <- classInt::findColours(class, col.pal)

  op <- par(
    oma = c(0, 0, 0, space), # Room for title and legend
    mfrow = c(1, 1) # Add more plots if you want
  )
  plot(sf, col = which.col, ...)
  title(main = paste(title))

  par(op) # Leave the last plot

  op <- par(usr = c(0, 1, 0, 1), # Reset the coordinates
            xpd = NA) # Allow plotting outside the plot region

  if (showlegend) {
    legend(placelegend, bty = "n",
           fill = rev(attr(which.col, "palette")),
           legend = rev(names(attr(which.col, "table"))))
  }
}
