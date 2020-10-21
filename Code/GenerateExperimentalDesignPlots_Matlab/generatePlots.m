function [] = generatePlots(M,destinationFolder,variableVector,colourStr,upAndDownStr,opacity,extraOption)
%% Defining parameters and allocating memory
% Allocate the three output vectors which are going to write
% to the file
polTime = zeros(length(variableVector),1); % The 50% quantile
polTimeHigh = zeros(length(variableVector),1); % The 95% quantile
polTimeLow = zeros(length(variableVector),1);% The 5% quantile
    
for index = 1: length(variableVector)    
    % Calculate the quantiles 
    polTime(index,1) = quantile(M(index,:),0.50);% 50% quantile
    polTimeHigh(index,1) = quantile(M(index,:),0.90);% 90% quantile
    polTimeLow(index,1) = quantile(M(index,:),0.05);% 5% quantile   
end
%% Save the files to the destination files
% Now we save our three quantiles! 
%-------------------------------------------------------------------------------------------------
% 50%-QUANTILE
%-------------------------------------------------------------------------------------------------
if extraOption
    % We define a string needed for plotting
    plotStr = ['\addplot[forget plot,color=',colourStr,'] coordinates {'];
else
    % We define a string needed for plotting
    plotStr = ['\addplot[color=',colourStr,'] coordinates {'];    
end
% We define a string needed for closing the statement
closeStr = '};';
% Open file, print to it, then close!
fileID = fopen(destinationFolder{1},'w');
fprintf(fileID,'%s\n',plotStr);
for i = 1:length(variableVector)
    fprintf(fileID,'\t\t(%0.4f\t,\t%0.4f\t)\n',variableVector(i,1),polTime(i,1));
end
fprintf(fileID,'%s\n',closeStr);
fclose(fileID);
%-------------------------------------------------------------------------------------------------
% 90%-QUANTILE
%-------------------------------------------------------------------------------------------------
% We define a string needed for plotting
plotStr = ['\addplot[forget plot,densely dashed,color=',colourStr,',name path=Up',upAndDownStr,'] coordinates {'];
% We define a string needed for closing the statement
closeStr = '};';
% Open file, print to it, then close!
fileID = fopen(destinationFolder{2},'w');
fprintf(fileID,'%s\n',plotStr);
for i = 1:length(variableVector)    
    fprintf(fileID,'\t\t(%0.4f\t,\t%0.4f\t)\n',variableVector(i,1),polTimeHigh(i,1));
end
fprintf(fileID,'%s\n',closeStr);
fclose(fileID);
%-------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------
% 5%-QUANTILE
%-------------------------------------------------------------------------------------------------
% We define a string needed for plotting
plotStr = ['\addplot[forget plot,densely dashed,color=',colourStr,',name path=Down',upAndDownStr,'] coordinates {'];
% We define a string needed for closing the statement
closeStr = '};';
% Open file, print to it, then close!
fileID = fopen(destinationFolder{3},'w');
fprintf(fileID,'%s\n',plotStr);
for i = 1:length(variableVector)
    fprintf(fileID,'\t\t(%0.4f\t,\t%0.4f\t)\n',variableVector(i,1),polTimeLow(i,1));
end
fprintf(fileID,'%s\n',closeStr);
% We also add a string to fill in colour between the higher and lower
% quantile
fillStr = ['\addplot[',colourStr,'!50,opacity=',num2str(opacity),',forget plot] fill between[of=Up',upAndDownStr,' and Down',upAndDownStr,'];'];
fprintf(fileID,'%s\n',fillStr);
fclose(fileID);
%-------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------

end