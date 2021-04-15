
%% This script accompanies the work 
% B. Saraswat, S. Ansumali, M. K Prakash, Using high effective risk of Adult-Senior duo in multigenerational homes to prioritize COVID-19 vaccination, medRxiv 2021.04.14.21255468; doi: https://doi.org/10.1101/2021.04.14.21255468

%% Initialization block
clear
%% Define Infection, vaccination scenarios
%---Analysis performed over 180 days --- 
days=1:1:180; 

%---Four different COVID-19 second wave scenarios are considered. Different
%peak values up to 400,000 daily infections are discussed in the manuscript
%--- Infection Scenario 1---
TotalInfections = 1.5e5*(1+0.667*days/30);
TotalInfections(31:60)=2.5e5;
i=1:1:60;
TotalInfections(61:120)=2.5e5*(1-i/60);
TotalInfections(121:180)=0;
clear i;
%---end Infection Scneario 1---

%--- Infection Scenario 2---
% TotalInfections = 1.5e5*(1+1.667*days/30);
% i=1:1:60;
% TotalInfections(31:90)=4e5*(1-i/60);
% TotalInfections(91:180)=0;
% clear i;
% --- end Infection scenario 2 ---

%--- Infection Scenario 3 ---
% TotalInfections = 1.5e5*(1+1.667*days/30);
% i=1:1:150;
% TotalInfections(31:180)=4e5*(1-i/150);
% clear i;
% --- end Infection Scenario 3 ---

%--- Infection Scenario 4 ---
% TotalInfections = 1.5e5*(1+1.667*days/60);
% i=1:1:120;
% TotalInfections(61:180)=4e5*(1-i/120);
% clear i;
% --- end Infection Scenario 4 ---


Vaccinations = 3e6*(1+0*days/max(days)); %--- daily rate of vaccinations. 2 million or 3 million daily doses administered
SAR = 0.3; % --- secondary attack rate within a family
Factor.InfectionReductionVaccinated = 0.33; %---67% reduction in chance of infection 
Factor.TransmissionFromVaccinated = 0.3;  %---70% reduction in transmission
Factor.MortalityReductionVaccinated = 0.2;%--- 80% reduction in mortality
Factor.SeniorVaccineCoverage= 0.5;        %--- Fraction of Seniors (age > 60) already vaccinated before beginning a planned roll-down to younger groups
Factor.AdultVaccineCoverage=0.2;          %--- Fraction of Adults (age 20-60) already vaccinated before beginning this analysis
Strategy='Targeted'; %--- Targeted vaccination by preferentially vaccinating Adults from multigenerational homes
%Strategy='Uniform'; %--- Uniform vaccination for any one in the age groups 20-60

%% Prepare demographics - Population by age, Defining Seniors, family linkage, 

days=180;
%---Age bin defines the lower age. For example Age 30 means individuals in
%the age group 30-34
Age = [20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];
%--Population corrected for 2011 to 2021 census
Population = 1.125*[111424222 101413965 88594951 85140684 72438112 62318327 49069254 39146055 37663707 26454983 19208842 9232503 6220229 2383167 1446534 633297 605778];
%Mortality risk obtained from the Meta-analysis of Levin et al. European Journal of Epidemiology (2020) 35:1123?1138
%https://doi.org/10.1007/s10654-020-00698-1
Mortality  = 1/100*[0.0064	0.0115	0.02	0.04	0.07	0.13	0.22	0.42	0.77	1.43	2.64	4.71	8.13	15.48	27.62	27.62	27.62];
%---
Seniors.Age = Age(9:17);            %---Seniors are defined as those above 60 years
Seniors.Population=Population(9:17);
Seniors.MortalityRisk =Mortality(9:17);
Seniors.Susceptible=zeros(days,9);
Seniors.Susceptible(1,:)=Seniors.Population;
Seniors.Deaths  = zeros(days,9);
Seniors.Infected= zeros(days,9);
Seniors.Vaccinated= zeros(days,9);

%--Assigning Adult With a Senior in the Family
AdultsWF.Age = Age(1:8);  %--ages between 20-60 --
AdultsWF.Population=zeros(1,8);
%---assigning children below (first generation differs by 25 years of age
%as discussed in Manuscript)
for i=1:5
    AdultsWF.Population(i+3)=AdultsWF.Population(i+3)+Seniors.Population(i);
end
%--assigning grandchildren below (second generation differs by 50 years of age
%as discussed in Manuscript)
for i=3:9
    AdultsWF.Population(i-2)=AdultsWF.Population(i-2)+Seniors.Population(i);
end
AdultsWF.MortalityRisk =Mortality(1:8);
AdultsWF.Susceptible=zeros(days,8);
AdultsWF.Susceptible(1,:)=AdultsWF.Population;
AdultsWF.Infected = zeros(days,8);
AdultsWF.VacInfected = zeros(days,8);
AdultsWF.Deaths = zeros(days,8);
AdultsWF.Vaccinated = zeros(days,8);

%---Assigning Adults with No Seniors in the Family ---
AdultsNF.Age = Age(1:8); %-- ages between 20-60
AdultsNF.Population=Population(1:8)-AdultsWF.Population;
AdultsNF.MortalityRisk =Mortality(1:8);
AdultsNF.Susceptible=zeros(days,8);
AdultsNF.Susceptible(1,:)=AdultsNF.Population;
AdultsNF.Infected = zeros(days,8);
AdultsNF.VacInfected = zeros(days,8);
AdultsNF.Deaths = zeros(days,8);
AdultsNF.Vaccinated = zeros(days,8);
%% A first roll-down group to ages 45-60 started on 1 April 2021. By the time our analysis finished on 12 April 2021, vaccination coverage in this group is 20% (maximum in some states)
% -- We remove this fraction of population already vaccinated

for j=6:8
    AdultsWF.Vaccinated(1,j) = Factor.AdultVaccineCoverage * AdultsWF.Susceptible(1,j);
    AdultsNF.Vaccinated(1,j) = Factor.AdultVaccineCoverage * AdultsNF.Susceptible(1,j);

    AdultsWF.Susceptible(1,j) = (1-Factor.AdultVaccineCoverage) * AdultsWF.Susceptible(1,j);
    AdultsNF.Susceptible(1,j) = (1-Factor.AdultVaccineCoverage) * AdultsNF.Susceptible(1,j);
end

%% Assign primary infections to AdultsWF and AdultsWF using the total infections in the population

for day=2:180
    for j=1:8

        if strcmp(Strategy,'Targeted')
        if sum(AdultsWF.Susceptible(day-1,:))>Vaccinations(day)
            AdultsWF.Vaccinated(day,j)= Vaccinations(day)*AdultsWF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:)));
        elseif sum(AdultsWF.Susceptible(day-1,:))>10
            AdultsWF.Vaccinated(day,j)= AdultsWF.Susceptible(day-1,j);
            AdultsNF.Vaccinated(day,j)= (Vaccinations(day)-sum(AdultsWF.Susceptible(day-1,:)))*AdultsNF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)));
        else
            AdultsNF.Vaccinated(day,j)= Vaccinations(day)*AdultsNF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)));
        end
       elseif strcmp(Strategy,'Uniform')
           AdultsWF.Vaccinated(day,j)= Vaccinations(day)*AdultsWF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)));
           AdultsNF.Vaccinated(day,j)= Vaccinations(day)*AdultsNF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)));
        end
       
        AdultsWF.Infected(day,j)    = TotalInfections(day)*AdultsWF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)) + Factor.InfectionReductionVaccinated*(AdultsWF.Vaccinated(day,j) + AdultsNF.Vaccinated(day,j)) );
        AdultsNF.Infected(day,j)    = TotalInfections(day)*AdultsNF.Susceptible(day-1,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)) + Factor.InfectionReductionVaccinated*(AdultsWF.Vaccinated(day,j) + AdultsNF.Vaccinated(day,j)) );
        AdultsWF.VacInfected(day,j) = TotalInfections(day)*Factor.InfectionReductionVaccinated*AdultsWF.Vaccinated(day,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)) + Factor.InfectionReductionVaccinated*(AdultsWF.Vaccinated(day,j) + AdultsNF.Vaccinated(day,j)) );
        AdultsNF.VacInfected(day,j) = TotalInfections(day)*Factor.InfectionReductionVaccinated*AdultsNF.Vaccinated(day,j)/(sum(AdultsWF.Susceptible(day-1,:) + AdultsNF.Susceptible(day-1,j)) + Factor.InfectionReductionVaccinated*(AdultsWF.Vaccinated(day,j) + AdultsNF.Vaccinated(day,j)) );

        
        AdultsWF.Susceptible(day,j) = AdultsWF.Susceptible(day-1,j) - AdultsWF.Infected(day,j) - AdultsWF.Vaccinated(day,j);
        AdultsNF.Susceptible(day,j) = AdultsNF.Susceptible(day-1,j) - AdultsNF.Infected(day,j) - AdultsNF.Vaccinated(day,j);

        AdultsWF.Deaths(day,j) = AdultsWF.MortalityRisk(j) * AdultsWF.Infected(day,j);
        AdultsNF.Deaths(day,j) = AdultsNF.MortalityRisk(j) * AdultsNF.Infected(day,j);
        
    end
end

%% Count intergenerational infection spreads

for day=1:180
%--assign infections
Factor.InfectionReductionVaccinated = 0.3;
Factor.TransmissionFromVaccinated = 0.3; 
    for j=1:5
        Seniors.Infected(day,j) = Seniors.Infected(day,j) + SAR * (AdultsWF.Infected(day,j+3) + Factor.InfectionReductionVaccinated * AdultsWF.VacInfected(day,j+3) ); % from children
    end
    for j=3:9
        Seniors.Infected(day,j) = Seniors.Infected(day,j) + SAR * (AdultsWF.Infected(day,j-2) + Factor.InfectionReductionVaccinated * AdultsWF.VacInfected(day,j-2) ); % from Grandchildren
    end
%--estimate deaths among seniors    
    for j=1:9
    Seniors.Deaths(day,j)       = Seniors.Infected(day,j)* (Factor.SeniorVaccineCoverage * Factor.MortalityReductionVaccinated + (1-Factor.SeniorVaccineCoverage))*Seniors.MortalityRisk(j);
    end
end
%% Graphical results

fprintf('Total number of fatalities among Seniors is %d\n', sum(sum(Seniors.Deaths)));

%-- Reduction of susceptible population
x=[0 30 60 90 120 150 180];
area(1/1e6*[sum(AdultsWF.Susceptible,2) sum(AdultsNF.Susceptible,2)]); axis([0 180 0 700]); 
set(gca,'XTick',x)
pause

%-- Deaths among Seniors
plot(cumsum(sum(Seniors.Deaths,2))); axis([0 180 0 200e3]); 
set(gca,'XTick',x)
