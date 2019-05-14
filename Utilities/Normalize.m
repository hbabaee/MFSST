function y = Normalize(X,Xmin,Xmax,shift)
[n , D] = size(X);
x = X - repmat(shift,n,1);
xmin = Xmin - shift;
xmax = Xmax - shift;

y = x;
for n=1:D
    y(:,n) = (x(:,n))/(xmax(n) - xmin(n));
end
