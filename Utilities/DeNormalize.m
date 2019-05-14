function y = DeNormalize(x,xmin,xmax,shift)
y = x*(xmax -xmin)+shift;