% Tristan Dicke
% Utility Model
% 20 May 2020
%%
clear;clc;
% Model 1 %Uncomment for model selection with defined non-MOC and MOC models
% R=64;
% q=17;
% Model 2
% R=64;
% q=23;
% Model 3
R=61;
q=18;
SS=R+q;
% Model 1 %Likelihoods of models based on gravity and water levels, uncomment for analysis of interest
% load LikeWaterlvlotherSYmoc; %Likelihood based on Water Level change(MOC as truth Models)
% Like=LikeWaterlvlotherSYmoc;
% load LikeWaterlvlotherSY; % Likelihoods based of Gravity change(Non-MOCs as truth Models)
% Like2=LikeWaterlvlotherSY;
% LikeAppend=cat(4,Like2,Like);

% Model 2
% load LikeWaterotherSYM2moc;
% Like=LikeWaterotherSYM2moc;
% load LikeWaterotherSYM2;
% Like2=LikeWaterotherSYM2;
% LikeAppend=cat(4,Like2,Like);
% load LikeGravityotherSYM2moc;
% Like=LikeGravityotherSYM2moc;
% load LikeGravityotherSYM2;
% Like2=LikeGravityotherSYM2;
% LikeAppend=cat(4,Like2,Like);

%Model 3
% load LikeWaterotherSYM3moc;
% Like=LikeWaterotherSYM3moc;
% load LikeWaterotherSYM3;
% Like2=LikeWaterotherSYM3;
% LikeAppend=cat(4,Like2,Like);

% load LikeGravityPrismotherSYM3fixedmoc;
% Like=LikeGravityPrismotherSYM3fixedmoc;
% load LikeGravityPrismotherSYM3fixed;
% Like2=LikeGravityPrismotherSYM3fixed;
% LikeAppend=cat(4,Like2,Like);
% LikeAppend=cat(4,Like2,Like);
%%
% Utility Function from POI's
for i=1:SS
    if (APOI(i)<0.5) %APOI is the POI array from water_level_means_other.m that is the POI's of models
        UtilityTruth(i)=(-0.02*(APOI(i))+1);
    elseif (APOI(i)>=0.5)&&(APOI(i)<=1.5)
        UtilityTruth(i)=(-0.82*APOI(i)+1.4);
    else (APOI(i)>1.5)
        UtilityTruth(i)=(-0.02*(APOI(i))+0.2);
    end
    Utility1(i)=UtilityTruth(i);
end
% Utility Value of each model in ensemble for selected AG well location
%% Constant Likelihood of models (needed to determine the Value of data)
% for x=1:SS
%    for k=1:SS
%         if (x==k)
%             LikeConstant(:,:,k)=zeros(50,50);
%         else    
%             LikeConstant(:,:,k)=zeros(50,50)+(1/(SS-1)); 
%         end
%    end
%     LikeConstant1(:,:,:,x)=LikeConstant(:,:,:);
%  end
%
% LikeAppend=LikeConstant1;
%% Likelihood Weighted Utilities
for k=1:SS
    Likesq=squeeze(LikeAppend(:,:,1:SS,k));
    for i=1:50
        for j=1:50
            POIModels(i,j,k)=sum(squeeze(Likesq(i,j,:)).*reshape(UtilityTruth,[SS,1]),'OmitNan');
        end
    end
end
%Likelihood Weighted Utility Values
%% Likelihood Weighted Utility Errors
for k=1:SS
    for i=1:50
        for j=1:50
            Pointloc(i,j,k)=abs(POIModels(i,j,k)-UtilityTruth(k));  %Looping over all Utilities minus the truth utility of the model
        end
    end
end
Pointloc=Pointloc(:,:,1:SS);
utilityofothermodels=Pointloc(:,:,1:SS);

%% Mean and Standard Derviation of Errors
for i=1:50
    for j=1:50

        MeanErrorUtility(i,j)=mean(utilityofothermodels(i,j,:));
        StdErrorUtility(i,j)=std(utilityofothermodels(i,j,:));

        MeanErrormoc(i,j)=mean(utilityofothermodels(i,j,(R+1):SS));
        StdErrormoc(i,j)=std(utilityofothermodels(i,j,(R+1):SS));

    end
end
% Utility Error Statistics saved in this step
MeanUtility3GstepP=MeanErrorUtility; %First number and letter altered when saving for different analyses
StdUtility3GstepP=StdErrorUtility;
MeanMOCU3GstepP=MeanErrormoc;
StdMOCU3GstepP=StdErrormoc;

%Data used to make graphics *simple subtraction is added to value of oberservation plots

