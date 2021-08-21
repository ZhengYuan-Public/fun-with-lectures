% Part I - Functions
% 1.1 Rectify the input matrix of arbitrary size.
function output = rectify(input_matrix)
output = zeros(size(input_matrix)); % Preallocation.
for i=1:numel(input_matrix)
        if input_matrix(i) > 0
            output(i) = input_matrix(i);
        end
end