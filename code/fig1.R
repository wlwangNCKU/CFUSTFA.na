library(mvtnorm)
library(mnormt)
library(relliptical)
library(rgl)
library(MASS)
library(ks)
library(misc3d)

################################################################################
#
#   Filename: fig1.R
#   Purpose: Generate 3D contour plots of CFUST distributions under different
#            lambda matrices to illustrate the influence of skewness structure
#   Input data files: function/cfustfn.R
#   Output data files: results/fig1.png
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, mnormt, relliptical, rgl, MASS, ks, misc3d
#
################################################################################

source(paste(SPATH, '/function/cfustfn.R', sep=''))

#True parameters
set.seed(123)
n = 500
mu.ast = c(0, 0, 0)
xi = mu.ast
Sigma = diag(3)
nu = 3
la1 = matrix(c(2, 2, -2, -2, 0, 0), 3, 2)
la2 = matrix(c(-2, 8, 0, 2, -2, -2), 3, 2)
la3 = matrix(c(-2, -2, 4, 0, -3, -2), 3, 2)
la4 = matrix(c(2, -2, -2, 2, 0, 0), 3, 2)
Y1 = rumst(n, xi, Sigma, la1, nu)
Y2 = rumst(n, xi, Sigma, la2, nu)
Y3 = rumst(n, xi, Sigma, la3, nu)
Y4 = rumst(n, xi, Sigma, la4, nu)

le <- rep(5e-9, 5)
colors = c('coral', '#00FF00', 'deeppink', '#00FFFF')
x=seq(-500, 500, length = 200)
y=seq(-500, 500, length = 200)
z=seq(-500, 500, length = 200)

open3d()
par3d(windowRect = c(50, 50, 1050, 1050), zoom = 6)
par3d(cex = 0.6)
mfrow3d(nr = 2, nc = 2, sharedMouse = TRUE)
bg3d("white")

make_lambda_title <- function(index, v1, v2) {
  bquote(
    "(" * .(letters[index]) * ") " * bold(Lambda)[.(index)] * " = {" ~
      "(" * .(v1[1]) * "," * .(v1[2]) * "," * .(v1[3]) * ")"^T * "," ~
      "(" * .(v2[1]) * "," * .(v2[2]) * "," * .(v2[3]) * ")"^T * "}"
  )
}

plot3d(Y1, box = FALSE, size = 0, type = 'n',
       xlim = range(Y1[,1]), ylim = range(Y1[,2]), zlim = range(Y1[,3]),
       xlab = '', ylab = '', zlab = '')
contour3d(function(x, y, z) dCFUST(cbind(x, y, z), mu.ast = xi, Sigma = Sigma, alpha = la1, nu = nu),
          level = le, x = x, y = y, z = z, add = TRUE, fill = TRUE, alpha = 0.3, color = colors[1])
view3d(theta = 75, phi = 20)
axes3d(edges = "bbox")
bgplot3d({
  plot.new()
  title(
    main = make_lambda_title(1, c(2, 2, -2), c(-2, -2, 0)),
    cex.main = 1.8, line = 1, font.main = 2
  )
})

next3d()
plot3d(Y2, box = FALSE, size = 0, type = 'n',
       xlim = range(Y2[,1]), ylim = range(Y2[,2]), zlim = range(Y2[,3]),
       xlab = '', ylab = '', zlab = '')
contour3d(function(x, y, z) dCFUST(cbind(x, y, z), mu.ast = xi, Sigma = Sigma, alpha = la2, nu = nu),
          level = le, x = x, y = y, z = z, add = TRUE, fill = TRUE, alpha = 0.3, color = colors[2])
view3d(theta = 75, phi = 20)
axes3d(edges = "bbox")
bgplot3d({
  plot.new()
  title(
    main = make_lambda_title(2, c(-2, 8, 0), c(2, -2, -2)),
    cex.main = 1.8, line = 1, font.main = 2
  )
})

next3d()
plot3d(Y3, box = FALSE, size = 0, type = 'n',
       xlim = range(Y3[,1]), ylim = range(Y3[,2]), zlim = range(Y3[,3]),
       xlab = '', ylab = '', zlab = '')
contour3d(function(x, y, z) dCFUST(cbind(x, y, z), mu.ast = xi, Sigma = Sigma, alpha = la3, nu = nu),
          level = le, x = x, y = y, z = z, add = TRUE, fill = TRUE, alpha = 0.3, color = colors[3])
view3d(theta = 75, phi = 20)
axes3d(edges = "bbox")
bgplot3d({
  plot.new()
  title(
    main = make_lambda_title(3, c(-2, -2, 4), c(0, -3, -2)),
    cex.main = 1.8, line = 1, font.main = 2
  )
})

next3d()
plot3d(Y4, box = FALSE, size = 0, type = 'n',
       xlim = range(Y4[,1]), ylim = range(Y4[,2]), zlim = range(Y4[,3]),
       xlab = '', ylab = '', zlab = '')
contour3d(function(x, y, z) dCFUST(cbind(x, y, z), mu.ast = xi, Sigma = Sigma, alpha = la4, nu = nu),
          level = le, x = x, y = y, z = z, add = TRUE, fill = TRUE, alpha = 0.3, color = colors[4])
view3d(theta = 75, phi = 20)
axes3d(edges = "bbox")
bgplot3d({
  plot.new()
  title(
    main = make_lambda_title(4, c(2, -2, -2), c(2, 0, 0)),
    cex.main = 1.8, line = 1, font.main = 2
  )
})

rgl.snapshot(paste(SPATH, '/results/fig1.png', sep=''), fmt = 'png')
