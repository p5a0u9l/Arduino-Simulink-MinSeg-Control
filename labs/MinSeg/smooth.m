% The purpose of smooth is to find the moving average of data. 
% Developer: HRLin 
function [output] = smooth(data,sampling)
n_data = length(data);
output = zeros(n_data,1);
for i = 1:1:n_data
    if i <= sampling
        output(i) = mean(data(1:i));     
    else
        output(i) = mean(data(i-sampling+1:i)); 
    end
end
        