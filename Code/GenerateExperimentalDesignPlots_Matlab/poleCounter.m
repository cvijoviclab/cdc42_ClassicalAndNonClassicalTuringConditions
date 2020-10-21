%nuOfPoles = xlsread('NumberOfPoles_Increasing_Gamma.xlsx');
index = 2;
if index ==1 % Classic
    nuOfPoles = xlsread('classic.xlsx');
    str = '\addplot[color=blue] coordinates {';
    fileID = fopen('classic.tex','w');
else % Non-classic
    nuOfPoles = xlsread('nonclassic.xlsx');    
    str = '\addplot[color=orange] coordinates {';
    fileID = fopen('nonclassic.tex','w');
end
fprintf(fileID,'%s\n',str);
[nuOfIter,c] = size(nuOfPoles);
for i = 1:nuOfIter
    fprintf(fileID,'(%d\t,\t%d)\n',nuOfPoles(i,1),nuOfPoles(i,2));     
end
fprintf(fileID,'};',str);
fclose(fileID);