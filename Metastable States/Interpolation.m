%%% This code is to interpolate the data, for each badger to get
%%% equidistant time

%%% Assumptions: travelling at the same speed, travelling in a straight
%%% line between points.

%%% tt = [0, 0.1, 0.2, ...]
%%% for each badger
%%% t = [ ]
%%% xx = interp1(t, x, tt,1)
%%% yy = interp1(t, y, tt, 1)

%%% Need columns with unique badger number (Number), easting (Easting),
%%% northing (Northing) and cumulative time (Cum_Time)
%%% read data in as table - then convert to numeric array

clear
clc

dataEDMD = readtable("file.csv"); % read in data
dataEDMD = [dataEDMD(:,"Number"),dataEDMD(:,"Easting"),dataEDMD(:,"Northing"),dataEDMD(:,"Cum_Time")];   

EDMD = table2array(dataEDMD);

clear dataEDMD % tidy workspace

%
IDchange = find(diff(EDMD(:,1))); %find when there are unique badgers
IDchange = [0 ; IDchange]; %add start
IDchange = [IDchange;length(EDMD)]; %add end

clc

interval = 20; % this will change based on how often you want to interpolate

for i = 1:1:length(IDchange)-1
    x = EDMD(IDchange(i)+1:IDchange(i+1),2);
    y = EDMD(IDchange(i)+1:IDchange(i+1),3); 
    t = EDMD(IDchange(i)+1:IDchange(i+1),4);
    tt = 0:interval:t(end); %interpolated for every 20 minutes
    xx = interp1(t, x, tt);
    yy = interp1(t, y, tt);
    Badger{i} = horzcat (tt',xx',yy'); %each badger interpolated put into individual cell
end

clear xx yy x y tt t 

%%% add column into each cell for each individual badger so it can be identified when put into one array
for i = 1:length(IDchange)-1 %change for how many individual badgers you have
    newrow = i * ones(length(Badger{i}),1);
    Badger{i} = [newrow ,Badger{i}];
end 

Badger_con = vertcat(Badger{:}); %put into one array

clear newrow i

% Extraction for EDMD in python

X1 = []; % first to second to last coordinates
Y1 = []; % second to last coordinates
for i = 1:length(IDchange)-1
    x = Badger{i}(1:end-1,3:4);
    x = x';
    X1 = [X1 , x];
    y = Badger{i}(2:end,3:4);
    y = y';
    Y1 = [Y1, y];
end

%% Save to put data into python

writematrix(transpose(X1),'X.csv')
writematrix(transpose(Y1),'Y.csv')
writematrix(Badger_con,'interpolate.csv')