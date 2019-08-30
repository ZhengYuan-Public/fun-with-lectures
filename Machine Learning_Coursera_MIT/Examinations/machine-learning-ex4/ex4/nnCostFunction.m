function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%% Implements the neural network cost function for a two layer neural network which performs classification
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)),hidden_layer_size, (input_layer_size + 1));%size = 25x401
Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end),num_labels, (hidden_layer_size + 1));%size = 10x26
% Setup some useful variables
m = size(X, 1);
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));
%% ====================== YOUR CODE HERE ======================
% --forward propagation for prediction#begin--
A1 = X;                                                                 %Original Input, size = 5000x400
A1_intercept = [ones(size(X,1),1),X];              %Add bias x_0 = 1, size = 5000x401
Z2 = A1_intercept * Theta1';                           %size = 5000x25
A2 = sigmoid(Z2);                                               %size = 5000x25
A2_intercept = [ones(size(A2,1),1), A2];        %Add bias a_0^(1), size = 5000x26
Z3 = A2_intercept * Theta2';                           %size = 5000x10
A3 = sigmoid(Z3);                                               %size = 5000x10
[~, Index] = max(A3');                                         
%[Y,I] = max(X) returns the indices of the maximum values in vector I.
%If the values along the first non-singleton dimension contain more
%than one maximal element, the index of the first one is returned.
prediction = Index(:);
% --forward propagation for prediction#end--
%% Cost function without regularization term
temp = A3';
h_x_vector = temp(:);             %A(:) show all elements(cloumn1, cloumn2, ...) in A in a column, size = (10*5000)x1
y_vector = [ ];
for i = 1:m
    y_vector = [y_vector 1:num_labels==y(i)];
    %Just learn this!
    %y = [1:num_labels == y(i)];
    %i =1, y(1) = 10, return a 10x1 matrix = [0 0 0 0 0 0 0 0 0 1]
    %i =2, y(2) = 10, return a 10x1 matrix = [0 0 0 0 0 0 0 0 0 1]
    %.........................
    %i =m=5000, y(5000) = 9, return a 10x1 matrix = [0 0 0 0 0 0 0 0 1 0]
    %size = 1x50000
end
%---cost function---
%J(¦È) = (-y*log(h_¦È(x))-(1-y)*(1-log(h_¦È(x)))/m;
J = (- y_vector*log(h_x_vector) - (1 - y_vector)*log(1 - h_x_vector)) / m;
%% --Regulariation term
Theta1_ni = Theta1(:, 2:end);   % Theta1 with no intercept terms. size = 25x400
Theta2_ni = Theta2(:, 2:end);   % Theta2 with no intercept terms. size = 10x25
regularization_term = lambda/(2*m) *(Theta1_ni(:)'*Theta1_ni(:) + Theta2_ni(:)'*Theta2_ni(:));
J = J + regularization_term;
%% Part 2: Implement the backpropagation algorithm to compute the gradients
Y = reshape(y_vector(:), size(A3)); % size = 5000 x 10
Delta3 = A3 - Y; % size = 5000 x 10
Theta2_grad = Theta2_grad + Delta3'*A2_intercept; % 10 x 26
Delta2_intercept = Theta2'*Delta3' .* sigmoidGradient([ones(1, size(Z2', 2)); Z2']);	% 26 x 5000
Delta2 = Delta2_intercept(2:end, :);                % 25 x 5000
Theta1_grad = Theta1_grad + Delta2*A1_intercept;   % 25 x 401
Theta1_grad = Theta1_grad ./ m;
Theta2_grad = Theta2_grad ./ m;
%% Part 3: Implement regularization with the cost function and gradients.
Theta1_regularization = lambda/m*[zeros(size(Theta1, 1), 1) Theta1(:, 2:end)];
Theta2_regularization = lambda/m*[zeros(size(Theta2, 1), 1) Theta2(:, 2:end)];
Theta1_grad = Theta1_grad + Theta1_regularization;
Theta2_grad = Theta2_grad + Theta2_regularization;
% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];
end
