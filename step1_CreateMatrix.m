% Before run this files:
% 1- Run cript "Simulate_Scenario.m" (to create the .mat files with the
% results for each Category
% 2- Save files in folder Data

%Restult: This script save a file with the mean and quantiles for Eisize,
%Hopitalization, ICU, Deaths for residents an Infected Staff
%This files is used in sptep 2_CreatePlots.m



for j = [0 1 3 5]    %This is because the files have the name Category0, Category1, Category3, and Category5

MeanEpisize = [];
MeanHosp    = [];
MeanICU     = [];
MeanStaff   = [];
MeanDeath   = []; 
QuanEpisize = [];
QuanHosp    = [];
QuanICU     = [];
QuanStaff   = [];
QuanDeath   = [];

for k=1:11  %11 is the number of R0
FileName    = ['./Data/Category' num2str(j) '_' num2str(k) '.mat'];
DATA        = matfile(FileName);

%Episize 
Episize     = DATA.vectorResR(length(DATA.vectorResR),:);
Hosp        = DATA.vectorResHT(length(DATA.vectorResHT),:);
ICU         = DATA.vectorResUT(length(DATA.vectorResUT),:);
Death       = DATA.vectorResD(length(DATA.vectorResD),:);
Staff       = DATA.vectorStaffR(length(DATA.vectorStaffR),:);

%resultsP=[resultsP;[convertCharsToStrings(F(k).name) MeanEpisize Quan_Episize Hosp Quan_Hosp ICU Quan_ICU Staff Quan_Staff]];

%MeanEpisize = [MeanEpisize [convertCharsToStrings(['R' num2str(k)])  mean(Episize)]];
%MeanHosp    = [MeanHosp    [convertCharsToStrings(['R' num2str(k)])  mean(Hosp)]];
%MeanICU     = [MeanICU     [convertCharsToStrings(['R' num2str(k)])  mean(ICU)]];
%MeanStaff   = [MeanStaff   [convertCharsToStrings(['R' num2str(k)])  mean(Staff)]];
%QuanEpisize = [QuanEpisize [convertCharsToStrings(['R' num2str(k)])  quantile(Episize,[0.25 0.975])]];
%QuanHosp    = [QuanHosp    [convertCharsToStrings(['R' num2str(k)])  quantile(Hosp,[0.25 0.975])]];
%QuanICU     = [QuanICU     [convertCharsToStrings(['R' num2str(k)])  quantile(ICU,[0.25 0.975])]];
%QuanStaff   = [QuanStaff   [convertCharsToStrings(['R' num2str(k)])  quantile(Staff,[0.25 0.975])]];

MeanEpisize = [MeanEpisize mean(Episize)];
MeanHosp    = [MeanHosp    mean(Hosp)];
MeanICU     = [MeanICU     mean(ICU)];
MeanDeath   = [MeanDeath   mean(Death)];
MeanStaff   = [MeanStaff   mean(Staff)];
QuanEpisize = [QuanEpisize quantile(Episize,[0.25 0.975])];
QuanHosp    = [QuanHosp    quantile(Hosp,[0.25 0.975])];
QuanICU     = [QuanICU     quantile(ICU,[0.25 0.975])];
QuanDeath   = [QuanDeath   quantile(Death,[0.25 0.975])];
QuanStaff   = [QuanStaff   quantile(Staff,[0.25 0.975])];
end

%Each row has the R0 value from R0_1 to R0_11, each column has the values
%saved avobe.
%Escenario0  = Baseline
%Escenario1  = Category1-2
%Escenario3  = Category3-4
%Escenario5  = Category5

save(['ResultsEscenario_' num2str(j) '.mat'], 'MeanEpisize', 'QuanEpisize', 'MeanHosp', 'QuanHosp', 'MeanICU', 'QuanICU', 'MeanStaff', 'QuanStaff', 'MeanDeath','QuanDeath');
end