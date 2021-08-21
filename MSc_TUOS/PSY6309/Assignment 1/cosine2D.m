% Part I - Functions
% 1.2 Creating a two-dimensional sine wave
%{
Q1:
>>[X]=meshgrid(-5:.2:5,-5:.2:5);
It generates a 2D meshgrid using vector x = [-5:.2:5] and y = [-5:.2:5].
Then the meshgrid created by x is assigned to X.

Q2:
>>imagesc(cosine2D(X,2));
It generates an image consistes of vertical bars whose color changes
periodically from blue to yellow then back to blue again. This is output
fits my anticipation since cosine is a periodic function and the result
also shows this property.
%}
function z = cosine2D(k, x)
    z = cos(k*x);
end
