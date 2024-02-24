##
## imagematrix class definition
##
## $Header: /database/repository/rimage/R/Attic/imagematrix.R,v 1.1.2.5 2004/03/17 06:35:18 tomo Exp $
##
## Copyright (c) 2003 Nikon Systems Inc.
## For complete license terms see file LICENSE



imagematrix <- function(mat, type = NULL, ncol = dim(mat)[1], nrow = dim(mat)[2],
                        noclipping = FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
  if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
  if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
  if (type != "rgb" && type != "grey") stop("Type is incorrect.")
  if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")

  imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)

  if (length(imgdim) == 3 && type == "grey") {
    # force to convert grey image
    mat <- rgb2grey(mat)
  }

  if (noclipping == FALSE && ((min(mat, na.rm = TRUE) < 0) || (1 < max(mat, na.rm = TRUE)))) {
    warning("Pixel values were automatically clipped because of range over.")
    mat <- clipping(mat)
  }

  mat <- array(mat, dim = imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix", class(mat))
  mat
}




plot.imagematrix <- function(x, significance_level = 0.05, colors = colorRampPalette(c("#FE9929", "#253494"))(100), ...) {
  if (is.null(attr(x, "type"))) stop("Type should be specified.")

  colmat <- matrix(0, nrow = nrow(x), ncol = ncol(x))

  # 
  colmat[x <= significance_level] <- x[x <= significance_level]
  colmat[x > significance_level] <- x[x > significance_level]

  
  color_palette <- colorRampPalette(colors)(100)
 
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1))
  par(mar = c(0, 0, 0, 0))
  par(oma = c(0, 0, 0, 0))
  par(omi = c(0, 0, 0, 0))
  # Imagen principal
  image(x = seq(0, ncol(colmat) - 1), y = seq(0, nrow(colmat) - 1),
        z = t(colmat[nrow(colmat):1, , drop = FALSE]), col = color_palette,
        xlab = "", ylab = "", axes = FALSE, asp = 1, ...)

  

  legend_image <- as.raster(matrix(color_palette, ncol = 1))
  plot(c(0, 15), c(0, 1), type = 'n', axes = F, xlab = '', ylab = '', main = '')
  rasterImage(legend_image,  1, 1, 0, 0)
  text(x = 3, y = seq(0, 1, l = 5), labels = sprintf("%.1f", seq(0, 1, l = 5)))

}








imageType <- function(x) {
  attr(x, "type")
}

rgb2grey <- function(img, coefs=c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
  if (length(dim(img))<3) stop("image matrix isn't rgb image.")
  imagematrix(coefs[1] * img[,,1] + coefs[2] * img[,,2] + coefs[3] * img[,,3],
              type="grey")
}
clipping <- function(img, low=0, high=1) {
  img[img < low] <- low
  img[img > high] <- high
  img
}
normalize <- function(img) {
  (img - min(img))/(max(img) - min(img))
}

# the end of file

### Added by Alejandro C. Frery
### 24 April 2014

imagematrixPNG <- function(x, name){
  dimensions <- dim(x)
  zero4 <- rep(0,4)
  png(file=name, width=dimensions[2], height=dimensions[1])
  par(mar=zero4, oma=zero4, omi=zero4)
  plot(x)
  dev.off()
}


imagematrixEPS <- function(x, name){
  dimensions <- dim(x)
  zero4 <- rep(0,4)
  postscript(file=name, width=dimensions[1], height=dimensions[2], paper = "special")
  par(mar=zero4, oma=zero4, omi=zero4)
  plot(x)
  dev.off()
}

# Finalmente! FunC'C#o que equaliza uma imagem de uma banda
equalize <- function(imagem) {
  imagemeq <- ecdf(imagem)(imagem)
  dim(imagemeq) <- dim(imagem)
  return(imagemeq)
  
}

### Added by Alejandro C. Frery
### 31 January 2018
# FunC'C#o que equaliza trC*s bandas independentemente
equalize_indep <- function(imagem) {
	req <- ecdf(imagem[,,1])(imagem[,,1])
	geq <- ecdf(imagem[,,2])(imagem[,,2])
	beq <- ecdf(imagem[,,3])(imagem[,,3])

  return(
  	imagematrix(
  		array(c(req, geq, beq), dim=dim(imagem))
    )
  )
}

### Added by Alejandro C. Frery
### 5 March 2018
# FunC'C#o que lineariza trC*s bandas independentemente
normalize_indep <- function(imagem) {
  rlin <- normalize(imagem[,,1])
  glin <- normalize(imagem[,,2])
  blin <- normalize(imagem[,,3])
  
  return(
    imagematrix(
      array(c(rlin, glin, blin), dim=dim(imagem))
    )
  )
}

### Histogram Matching
### Purpose: forces an output (destination) image to have the same histogram
### as the input (reference) image
### Works on images with one band
### Returns an imagematrix object
### Maceio, 23 August 2014
### Alejandro C. Frery 

HistogramMatching <- function(reference, destination) {
  
  return(imagematrix(
    HistToEcdf(
      hist(reference, breaks='FD', plot=FALSE), inverse=TRUE)(equalize(destination)),
    ncol=dim(reference)[1], nrow=dim(reference)[2], type='grey')
  )
  
}

# Copyright 2013 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: mstokely@google.com (Murray Stokely)


