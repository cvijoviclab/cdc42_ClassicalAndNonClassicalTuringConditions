function outputMatrix = readData(dataStr,Case,nuOfRepititions)
% Read data
pathStr = '../../Results/';
% Allocate memory for the outputMatrix
outputMatrix = zeros(5,5);
% Loop over all repititions, read data and gather in the output matrix 
for index = 1:nuOfRepititions
    data = readtable([pathStr,Case,num2str(index),dataStr]);
    temp = table2array(data);
    if index == 1
        outputMatrix = temp(2:end,2);            
    else
       outputMatrix = [outputMatrix, temp(2:end,2)]; 
    end
end

end