##               COPYRIGHT NOTICE
##
##  Copyright (C) 2014 Edwin L. (Ted) Woollett
##  http://www.csulb.edu/~woollett
##
##  This program is free software; you can redistribute
##  it and/or modify it under the terms of the
##  GNU General Public License as published by
##  the Free Software Foundation; either version 2 
##  of the License, or (at your option) any later version. 
##
##  This program is distributed in the hope that it
##  will be useful, but WITHOUT ANY WARRANTY;
##  without even the implied warranty of MERCHANTABILITY
##  or FITNESS FOR A PARTICULAR PURPOSE. See the 
##  GNU General Public License for more details at
##  http://www.gnu.org/copyleft/gpl.html



##   cp1code.R 
##   Code used in:
##     Computational Physics with R and Maxima,
##     Chapter 1, Numerical Differentiation, 
##        Quadrature, and Roots
##
##     Use source("c:/k1/cp1code.R") (for example)
##        to load code definitions into R  

## sec. 1.2  R Interfaces, R Language

mytest = function() {
     xv = seq(0,1,by = 0.25)
     yv = sin(xv)
     cat(" xv = ",xv," yv = ",yv)
     zv = exp(yv)
     cat(" zv = ",zv)
     xv = (1:4)/4
     cat(" xv = ", xv)
     yv = sin(xv)
     cat( " yv = ", yv ," \n\n")}
     
## sec. 1.4.2  Simple Numerical Derivative Methods

tryh = function() {
     x = 1
     exact = cos(x)     
     cat(" enter h <= 0 to stop \n")
     repeat {
       h = as.numeric(readline(" input h: "))
       if (h <= 0) break       
       fprime = 0.5*(sin(x+h) - sin(x-h))/h
       diff = exact - fprime
       cat(" h = ",h," error = ",diff," \n\n")}}
       

D1c = function (func, x, h = 1e-8) {
        (func(x + h) - func(x - h))/(2*h) }
        
D1f = function(func,x,h = 1e-8) {
        (func(x + h) - func(x))/h }
        
D1b = function(func,x,h = 1e-8) {
        (func(x) - func(x - h))/h }
        
## sec. 1.5.5   Trapezoidal Rule: Uniform grid

trap = function(xv,yv){
   n = length(xv)
   (xv[2]-xv[1])*((yv[1]+yv[n])/2 + sum(yv[2:(n-1)]))}
   
trap2 = function(func, a, b, N) { # N is number of panels, N+1 is length of xv
   h = (b - a)/N
   xv = seq(a, b, by = h)
   yv = func(xv)
   h*((yv[1]+yv[N+1])/2 + sum(yv[2:N]))}
   
## sec. 1.5.7  Trapezoidal Rule: non-uniform grid

trapz = function(xv,yv) {
   idx = 2:length(xv)
   as.numeric( (xv[idx] - xv[idx - 1]) %*% (yv[idx - 1] + yv[idx])/2)}
   
## sec. 1.5.9   Simpson's 1/3 Rule

simp = function(xv,yv) {
   N = length(xv) - 1
   if (N %% 2 != 0) {
        return(" length(xv) should be an odd integer ")}        
   h = xv[2] - xv[1]
   if (N == 2) s = yv[1]+4*yv[2]+yv[3]  else {   
       s = yv[1] + yv[N+1] + 4*sum(yv[seq(2,N,by=2)]) +
                2*sum(yv[seq(3,N-1,by=2)])}
   s*h/3}
   
simp2 = function(func, a, b, N) {
    if (N %% 2 != 0) {
        return(" N should be an even integer ")}        
    h = (b - a)/N
    xv = seq(a, b, by=h)
    yv = func(xv)
    if (N == 2) s = yv[1]+4*yv[2]+yv[3]  else {   
          s = yv[1] + yv[N+1] + 4*sum(yv[seq(2,N,by=2)]) +
                2*sum(yv[seq(3,N-1,by=2)])}
    s*h/3}
    
##  sec. 1.6.4  Newton-Raphson root


newton = function(f,x0,eps = 1e-5,h = 1e-4,small = 1e-14) {
    xn = x0
    repeat {
       if ( abs(f(xn)) < eps ) return(xn)
       s = (f(xn+h) - f(xn-h))/(2*h)
       if (abs(s) <= small) {
          cat(" derivative at",xn,"is zero \n")
          break}
       xn = xn - f(xn)/s}}
       
## sec. 1.6.5    Secant search for  root

secant = function(f,x0,x1,eps = 1e-5, small = 1e-14) {
    repeat {
      if ( abs(f(x1)) < eps ) return(x1)      
      s = f(x1) - f(x0)
      if (abs(s) <= small) {
           cat(" derivative near",(x0+x1)/2,"is zero \n")
           break}
      x2 = x1 - f(x1)*(x1 - x0)/s
      x0 = x1
      x1 = x2}}
      
##  sec. 1.6.10  divide and conquer root search

rtsearch = function(func,x,dx,xacc) {
    fold = func(x)
    repeat {
      if (abs(dx) <= xacc) break   
      x = x + dx   
      if (fold*func(x) < 0) {
          x = x - dx
          dx = dx/2}}
    x}
    
    