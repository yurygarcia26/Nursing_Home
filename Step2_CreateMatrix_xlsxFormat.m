%In this code, a file is generated with the matrix that contains the average
%and the quantiles (95%)for the total of patients, hospitalized, recovered, etc., 
%The name of the file is MatrixLatexEscenario_#.mat

%scenarios 0 = Baselne, 1= Category1, 3 = Category 3, 5 = Category5 




resultsP = [];
Matrix   = [];

for j = 5    % keep j=# ,where # is the scenario number taht you want to see
             % j= [0, 1, 3, 5]  if you want to save all files at the same time

MeanEpisize = [];
MeanHosp    = [];
MeanICU     = [];
MeanStaff   = [];
MeanDeath   = [];
QuanEpisize = [];
QuanHosp    = [];
QuanICU     = [];
QuanStaff   = [];
Death       = [];

for k=1:11  %1:11 is the number of R0
FileName    = ['./Data/Category' num2str(j) '_' num2str(k) '.mat'];
DATA        = matfile(FileName);

%Episize 
Episize     = DATA.vectorResR(length(DATA.vectorResR),:);   %recovery people
Hosp        = DATA.vectorResHT(length(DATA.vectorResHT),:); %Total of Hospitalized people
ICU         = DATA.vectorResUT(length(DATA.vectorResUT),:);
Staff       = DATA.vectorStaffR(length(DATA.vectorStaffR),:);
Death       = DATA.vectorResD(length(DATA.vectorResD),:);

MeanEpisize = [MeanEpisize mean(Episize)];
MeanHosp    = [MeanHosp    mean(Hosp)];
MeanICU     = [MeanICU     mean(ICU)];
MeanDeath   = [MeanDeath   mean(Death)];
MeanStaff   = [MeanStaff   mean(Staff)];
QuanEpisize = [QuanEpisize quantile(Episize,[0.25 0.975])];
QuanHosp    = [QuanHosp    quantile(Hosp,[0.25 0.975])];
QuanICU     = [QuanICU     quantile(ICU,[0.25 0.975])];
QuanStaff   = [QuanStaff   quantile(Staff,[0.25 0.975])];
R0 = 2:0.2:4;

resultsP=[resultsP; [R0(k)  quantile(Episize,0.25) mean(Episize) quantile(Episize,0.975) quantile(Hosp,0.25) mean(Hosp) quantile(Hosp,0.975) quantile(ICU,0.25) mean(ICU) quantile(ICU, 0.975) quantile(Death,0.25) mean(Death) quantile(Death, 0.975)  quantile(Staff,0.25) mean(Staff) quantile(Staff, 0.975)]];
Matrix   =[Matrix; [R0(k)  mean(Episize) mean(Hosp) mean(ICU) mean(Death)]];

end
xlswrite("Category_5",Matrix);
xlswrite("ResultsCat5",resultsP);
%save(['MatrixLatexEscenario_' num2str(j) '.mat'], 'MeanEpisize', 'QuanEpisize', 'MeanHosp', 'QuanHosp', 'MeanICU', 'QuanICU', 'MeanStaff', 'QuanStaff');
end