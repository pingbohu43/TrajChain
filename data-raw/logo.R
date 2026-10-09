# Creates man/figures/logo.png, the package logo (run from the package root).
png("man/figures/logo.png", width = 518, height = 600, bg = "transparent", res = 144)
op <- par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(-1, 1), ylim = c(-1.155, 1.155), asp = 1, axes = FALSE, xlab = "", ylab = "")
ang <- seq(90, 450, by = 60) * pi / 180
hx <- 1.1 * cos(ang)
hy <- 1.1 * sin(ang)
polygon(hx, hy, col = "#16365c", border = "#2a78d6", lwd = 7)
# naive age curve (bends down with selection)
xs <- seq(-0.62, 0.62, length.out = 100)
naive <- -0.05 + 0.36 * tanh(3 * xs) - 1.1 * pmax(xs + 0.05, 0)^2
lines(xs, naive, col = "#eb6834", lwd = 4, lty = 2)
# chained shape function: grid points linked by ratios
xg <- seq(-0.62, 0.62, length.out = 7)
yg <- -0.05 + 0.36 * tanh(3 * xg)
segments(xg[-7], yg[-7], xg[-1], yg[-1], col = "white", lwd = 7)
points(xg, yg, pch = 21, bg = "#2a78d6", col = "white", cex = 2.3, lwd = 3)
text(0, -0.62, "TrajChain", col = "white", cex = 1.65, font = 2)
par(op)
dev.off()
