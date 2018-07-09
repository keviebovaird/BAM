%%%%%%%%%%%%%%%%%BAM TEST PULSE ANALYSIS%%%%%%%%%%%%%%%
%This is the *final* version of the analysis
%Yay
%Kevie Bovaird, SLAM Lab, The Ohio State University
%July 2018

%% Set-up Nonsense
clear all; clc;

%% Rhythm Analysis
%% Load Files

results = cd; 
rhyResults = fopen('rhythmData.txt');
rhy = textscan(rhyResults, '%d %s %s %s');
fclose(rhyResults); 

%% Pre-Allocate Variables

rhySubNum = rhy{1};
rhySubInit = rhy{2};
rhyStim = rhy{3};
rhyResp = rhy{4};
rhyCorr = string();


%% Analyze Data

%write a function to turn stimulus name into binary
for rr= 1:length(rhyStim); 
     if contains (rhyStim{rr}, 'D');
         rhyCorr(rr) = "different";
     elseif contains (rhyStim{rr},'S');
         rhyCorr(rr) = "same";
     else
         rhyCorr(rr) = "Null";
     end
end

rhyCorr = rhyCorr';
rhyNew= [rhySubNum, rhyCorr, rhyResp]; %I forget if this is needed
tt= 1:length(rhyStim);
rhyScore= strcmp (rhyResp(tt), rhyCorr(tt));

%Sorting by Subject Number
rhySubNumSort = unique(rhySubNum, 'rows'); %this variable is every subject number I have without repeats
numSubjects = length(rhySubNumSort); %number of subjects. Useful to have
rhyCount= nan(1, numSubjects);

%Pull out Counted Score
for mm = 1:numSubjects
        rhyCount(mm) = sum(rhyScore(rhySubNumSort(mm)==rhySubNum));
end
rhyCount= rhyCount/20;

%the table. do I need a table? Eh. It's a table!
rhyScoreTable= table(rhySubNumSort, rhyCount');
rhyStimTable= table(rhySubNum, rhyStim, rhyScore);

%%%%%%%%%%%%%%%%%%Melody Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Load Files

melResults = fopen('melodyData.txt');
mel = textscan(melResults, '%d %s %s %s');
fclose(melResults); 

%% Pre-Allocate Variables

melSubNum = mel{1};
melSubInit = mel{2};
melStim = mel{3};
melResp = mel{4};
melCorr = string();


%% Analyze Data

%write a function to turn stimulus name into binary
for qq= 1:length(melStim); 
     if contains (melStim{qq}, 'D');
         melCorr(qq) = "different";
     elseif contains (melStim{qq},'S');
         melCorr(qq) = "same";
     else
         melCorr(qq) = "Null";
     end
end

melCorr = melCorr';
melNew= [melSubNum, melCorr, melResp]; %I forget if this is needed
ss= 1:length(melStim);
melScore= strcmp (melResp(ss), melCorr(ss));

%Sorting by Subject Number
melSubNumSort = unique(melSubNum, 'rows'); %this variable is every subject number I have without repeats
mnumSubjects = length(melSubNumSort); %number of subjects. Useful to have
melCount= nan(1, mnumSubjects);

%Pull out Counted Score
for nn = 1:mnumSubjects
        melCount(nn) = sum(melScore(melSubNumSort(nn)==melSubNum));
end
melCount= melCount/20;

%the table again
melScoreTable= table(melSubNumSort, melCount');
melStimTable= table(melSubNum, melStim, melScore);

%% Pulse Section
%%Let's read in the data and do some prep work, yes?
filename = 'pulseData.txt';

fid = fopen(filename);
j = 1;

NumOfBins=20;                       % change this to changebucket size 20 = 5%
xx= -NumOfBins/2:1:NumOfBins/2;

while true  %Reads file into DataCell Variable
        tline = fgetl(fid);
        if tline == -1
            break  % exits while when end of file is reached
        end
        tempHeaders = {};
        tempData = [];
        for i = 1:3
            [tempHeaders{i},tline] = strtok(tline); %tbh I forget why this works
        end
        tempData = str2num(tline);
        
        DataCell{j,1} = tempHeaders{1};    %Subject Number
        DataCell{j,2} = tempHeaders{2};    %Subject Initials
        DataCell{j,3} = tempHeaders{3};    %Pulse File
        DataCell{j,4} = tempData;          %Times in mseconds
        
        
        %Math Stuffs
        workablePulseTimes= cell2mat(DataCell(j,4)); %Let's prune off the last few seconds for accuracy
        prunePulses = find(workablePulseTimes<20001);
        DataCell{j,4} = workablePulseTimes(prunePulses);
        
        pulseString = DataCell{j,3};        % get the BPM from Pulse File Name
        [pnum,pulseString] = strtok(pulseString,'_');
        [bpm, pulseString] = strtok(pulseString,'_');
        BeatDelta_ms= 60000/(str2num(bpm));   %Beat in milliseconds
        BinSize_ms = BeatDelta_ms/NumOfBins;   %Bin Size = Beat in ms/Number of bins
        DataCell{j,5}=BeatDelta_ms;
        DataCell{j,6}=(DataCell{j,4})/BeatDelta_ms;   % Pulse Time/Beats in ms
        DataCell{j,7}=rem(DataCell{j,6},1);       % gets the decimal portion of the number 
        DataCell{j,8}=mean(DataCell{j,7});         % mean of the the pulse
        DataCell{j,9}=(DataCell{j,7}-DataCell{j,8})*360;
        DataCell{j,10}=std(DataCell{j,7});          % stdeviation of the pulse mean
        DataCell{j,11} = min(DataCell{j,6});   % How fast they got the beat
        DataCell{j,12} = length(DataCell{j,6});   %How many spacebar presses (FLAG if < 5?) as not enough data to represent skill??
        DataCell{j,13} = round((20000/DataCell{j,5})-0.1); %Total possible # of beats
        DataCell{j,14}= length(unique(round(DataCell{j,6}+(1-DataCell{j,8})))); %Total number of unique beats pressed
        DataCell{j,15}= DataCell{j,13}-DataCell{j,14}; %Misses
        DataCell{j,16}= length(DataCell{j,6})-(DataCell{j,14}); %False Alarms
        DataCell{j,17}= DataCell{j,14}/DataCell{j,13}; %Percentage of beats gotten
        
        %DoubleTime Nonsense
        DataCell{j,18}= BeatDelta_ms/2;             
        DataCell{j,19}=DataCell{j,4}/DataCell{j,18}; %PulseTime/Beats for Doubletime
        DataCell{j,20}=rem(DataCell{j,19},1);
        DataCell{j,21}=mean(DataCell{j,20});
        DataCell{j,22}= round(DataCell{j,6}); %The beat # for reg time
        DataCell{j,23}= diff(cell2mat(DataCell(j,22))); %math
        DataCell{j,24}= round(DataCell{j,19}); %the beat # for double time
        DataCell{j,25}= diff(cell2mat(DataCell(j,24))); %math
        %Let's call out those troublesome two-timers 
        DataCell{j,26}= (diff(cell2mat(DataCell(j,22))))==1;
        DataCell{j,27}= (diff(cell2mat(DataCell(j,24))))==1;
        DataCell{j,28}= sum(DataCell{j,26})<sum(DataCell{j,27}); %logical index of doubletimers (=1)
       
        %Generate Score
        DataCell{j,29}=BinSize_ms;
        DataCell{j,30}=round(((DataCell{j,7}-DataCell{j,8})*BeatDelta_ms)/BinSize_ms);
       %Replace Data with Corrected DoubleData
        if DataCell{j,28}==1
            DataCell{j,30}=round(((DataCell{j,20}-DataCell{j,21})*DataCell{j,18})/BinSize_ms);
        else
            DataCell{j,30};
        end
        
        PlotVector=DataCell{j,30};
        YCount = zeros(1,NumOfBins+1);    % Counts the number of values for each bin
        for k=1 : length(PlotVector)
            for m=1:NumOfBins+1
               if PlotVector(k) == xx(m) 
                    YCount(m)=YCount(m)+1;
               end
           end
        end
        DataCell{j,31} = length(DataCell{j,6});   %How many spacebar presses
        PulseScore=0;
        for m=1:NumOfBins+1
               PulseScore=PulseScore +YCount(m)*1/(abs(xx(m))+1);
        end
        DataCell{j,32} = PulseScore/DataCell{j,31};   %PulseScore/Spacebar presses (Avg Score space bar press )
        
        j = j+1;
end



%%Okay it's score time
subjectNumbers=str2double(DataCell(1:j-1,1));
subjectNumbersUnique=unique(str2double(DataCell(1:j-1,1)));
beatScore=[];
pulseScoreAv=[];
totalPulseScore=[];

% The Average % of Beats Gotten out of the total per subject
for ss=1:length(subjectNumbersUnique)
    %make this per subject
    currentSubject = subjectNumbersUnique(ss) == subjectNumbers;
    %Beat Score per subject
    beatScore(ss)=sum(cell2mat(DataCell(currentSubject,17)))/length(cell2mat(DataCell(currentSubject,17)));
    pulseScoreAv(ss)=mean(cell2mat(DataCell(currentSubject,32)));
    totalPulseScore(ss)= beatScore(ss)*0.1 + pulseScoreAv(ss)*0.9;
end

%%Let's make a table
MyPulseTable=table(subjectNumbersUnique,beatScore',pulseScoreAv',totalPulseScore');
MyPulseTable.Properties.VariableNames{1}='SubjectNo';
MyPulseTable.Properties.VariableNames{2}='PercentBeats';
MyPulseTable.Properties.VariableNames{3}='Consistency';
MyPulseTable.Properties.VariableNames{4}='PulseScore';


%% Total Analysis 

scoreBAM= [];

for zz= 1:length(subjectNumbersUnique)
    
    scoreBAM(zz)= (1/3)*(rhyCount(zz))+(1/3)*(melCount(zz))+(1/3)*(totalPulseScore(zz));

end
%Make a Table of EVERYTHING
scoreBAMTable= table(subjectNumbersUnique, rhyCount',melCount',totalPulseScore', scoreBAM');
scoreBAMTable.Properties.VariableNames{1}='SubjectNo';
scoreBAMTable.Properties.VariableNames{2}='RhythmScore';
scoreBAMTable.Properties.VariableNames{3}='MelodyScore';
scoreBAMTable.Properties.VariableNames{4}='PulseScore';
scoreBAMTable.Properties.VariableNames{5}='BAMScore';

%Write the Table to a File for R Analysis
writetable(scoreBAMTable, 'BAMScore.csv');
%% Code Graveyard

%Okay, lets find those hella low scores
%pulseScoresTot = cell2mat(DataCell(a:b,16));
%pulseScoresLow = pulseScoresTot < 30;
%pulseScoresInvestigate = pulseScoresTot(pulseScoresLow);
%Radian bullshit
%Okay, so take my pulse time stamps, divide by bpm, multiply by 360