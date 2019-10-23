function [J, grad] = costFunctionReg(theta, X, y, lambda)
%COSTFUNCTIONREG Compute cost and gradient for logistic regression with regularization
%   J = COSTFUNCTIONREG(theta, X, y, lambda) computes the cost of using
%   theta as the parameter for regularized logistic regression and the
%   gradient of the cost w.r.t. to the parameters. 

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta.
%               You should set J to the cost.
%               Compute the partial derivatives and set grad to the partial
%               derivatives of the cost w.r.t. each parameter in theta
h_theta = sigmoid(X*theta);
reg_theta = theta(2:end);
sum_reg_theta_sqr = reg_theta'*reg_theta;
J = (-y'*log(h_theta)-(1-y)'*log(1-h_theta))/m+lambda/2/m*sum_reg_theta_sqr;

grad_theta_temp = ((h_theta-y)'*X)'/m;
grad_theta1 = grad_theta_temp(1);
grad_theta_2toend = grad_theta_temp(2:end)+lambda/m*reg_theta;
grad = [grad_theta1; grad_theta_2toend];

% =============================================================

end
