isomap2 <- function(sf, values, breaks, x.pred, y.pred, n = 9, name = "PuBu", divcol = FALSE, placelegend = "right", titlelegend = "", space = 8, places = NULL, names.places = "", col.places = "black", showlegend = TRUE, ...){

  class <- classInt::classIntervals(values,
                                    style = "fixed",
                                    fixedBreaks = breaks)
  pal <- RColorBrewer::brewer.pal(n, name)
  if (divcol) {pal <- rev(pal)}
  col.pal <- colorRampPalette(pal)((length(breaks)-1))
  which.col <- classInt::findColours(class, col.pal)

  grid <- expand.grid(sort(unique(x.pred)), sort(unique(y.pred)))
  x.grid <- grid[ ,1]
  y.grid <- grid[ ,2]
  cells <- c(length(unique(x.grid)), length(unique(y.grid)))
  values.grid <- rep(NA, prod(cells))

  for (i in 1:(length(x.pred))){
    id <- which(x.grid == x.pred[i] & y.grid == y.pred[i])
    values.grid[id] <- values[i]
  }

  Show <- matrix(values.grid, cells[1], cells[2])

  op <- par(
    oma = c(0, 0, 0, space), # Room for title and legend
    mfrow = c(1, 1) # Add more plots if you want
  )

  plot(sf)
  image(unique(x.grid), unique(y.grid), Show, col = col.pal,
        breaks = breaks, xlab = "", ylab = "", axes = FALSE,
        add = TRUE)
  plot(sf, col = NA, add = TRUE, ...)
  points(places, col = col.places, pch = 15, cex = 1)
  calibrate::textxy(places[, 1], places[, 2], labs = names.places, col = col.places, cex = 1, font = 1.5)

  par(op) # Leave the last plot

  op <- par(
    usr = c(0, 1, 0, 1), # Reset the coordinates
    xpd = NA # Allow plotting outside the plot region
  )

  if (showlegend) {
    legend(placelegend, bty = "n",
           fill = rev(attr(which.col, "palette")),
           legend = rev(names(attr(which.col, "table"))),
           title = titlelegend)
  }
}
