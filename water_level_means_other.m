%% Name: Tristan Dicke
%Date: 7 April 2020
%Project: Data Collected for Heads- Gravity and Water Levels

clc;clear; % clearing workscape so extra variables are not over written 
%%
gam=(6.67430e-11)*1e5; %Gravitational Constant
p=0.1;  %p/Sy=porosity/Specific yield of changing cells 
d=997; %d=density water
load base_top_elev.csv; %Elevation Profile
%% Model Selection (Uncomment the s and q values)
% s=number of non MOC models
% q=number of MOC models
% w=total number of models in ensemble from first script
%Model 1
% s=64; 
% q=17;
%Model 2
% s=64;
% q=23;
% w=s+q;
%Model 3
s=61;
q=18;
w=s+q;
%% Uploading MOC into 3D Matrix
%projectdir = 'C:\Users\Tristan Dicke\Documents\Gravity Research Papers\MOC_s';
projectdir = 'C:\Users\Tristan Dicke\Documents\Gravity Research Papers\ConstantSYModel3Moc';
%projectdir=folder of head outputs for ensemble generation for MOC models
dinfo = dir( fullfile(projectdir, '*.csv') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );
precipsInt = cell(num_files, 1);
for K = 1 : num_files
  this_file = filenames{K};
  x = load(this_file);
  precipsInt{K} = x(2:end,:);
  %precipsInt{1,K}(isnan(mydata{1,K}))=0;
end
Headprofilesmocs2 = cat(3,precipsInt{:,1});
%Headprofilesmocs2=MOC Heads of models in ensemble(listed(NTNA),(YTNA),(YTYA))
%% Uploading Other Models
%projectdir = 'C:\Users\Tristan Dicke\Documents\Gravity Research Papers\AcceptedModels';
projectdir = 'C:\Users\Tristan Dicke\Documents\Gravity Research Papers\ConstantSYModel3';
%projectdir=folder of head outputs for ensemble generation for non-MOC models
dinfo = dir( fullfile(projectdir, '*.csv') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );
precipsInt = cell(num_files, 1);
for K = 1 : num_files
  this_file = filenames{K};
  x = load(this_file);
  precipsInt{K} = x(2:end,:);
  %precipsInt{1,K}(isnan(mydata{1,K}))=0;
end
Headprofilesmocs = cat(3,precipsInt{:,1});
%Headprofilesmocs=Non-MOC Heads of models in ensemble(listed(NTNA),(YTNA),(YTYA))
%% Set up of Arrays for initial Analysis
n=s; %Array indexs for non-MOCs
ytya_array=[1:n];
ytya_array=(ytya_array*3-2+2);
ntna_array=[1:n];
ntna_array=(ntna_array*3-2);
ytna_array=[1:n];
ytna_array=(ytna_array*3-2+1);


m=q; %Array indexs for MOCs
mntna_array=[1:m];
mntna_array=(mntna_array*3-2);
mytya_array=[1:m];
mytya_array=(mytya_array*3-2+2);
mytna_array=[1:m];
mytna_array=(mytna_array*3-2+1);

ind=ytya_array;
ind2=mytya_array;
nind=ytna_array;
nind2=mytna_array;
% Drawdown within Town Well (Can change the location for POI's at different points within model)
POI=Headprofilesmocs(22,39,ind); %Point location of interest in model-town well(22,39)
POI_1=Headprofilesmocs(22,39,nind);
POI2=Headprofilesmocs2(22,39,ind2);
POI_2=Headprofilesmocs2(22,39,nind2);
POI=reshape(POI,1,s);
POI2=reshape(POI2,1,q);
POI_1=reshape(POI_1,1,s);
POI_2=reshape(POI_2,1,q);
POI=POI_1-POI;
POI2=POI_2-POI2;
APOI=[POI POI2]; % POI for all models (non-mocs and mocs)
APOI=reshape(APOI,w,1); 

%% Non-MOC Initial Water Level/Gravity chnage Arrays
 % Setting Up Loops for Conditions
    n=s; % arrays for truth models(can be an s or q depending on models investigated)
    ntna_array=[1:n];
    ntna_array=ntna_array*3-2;
    ytna_array=[1:n];
    ytna_array=(ytna_array*3-2+1);
    ytya_array=[1:n];
    ytya_array=(ytya_array*3-2+2);
    drho=997*(0.1); %m
  %Locations of Measurement  
    Xc=[500:1000:49500];%  Instrument location X-array
    Yc=[500:1000:49500];% Instrument location Y-array
    Xc1=[0:1000:49000];Xc2=[1000:1000:50000]; %X1 and X2 dimension arrays for Prism Formula
    Yc1=[0:1000:49500];Yc2=[1000:1000:50000]; %Y1 and Y2 dimension arrays for Prism Formula
    
    for k=1:s  % Number of non-MOC model looped(can be an s or q depending on models investigated)
        id=ntna_array(k); %indexing for ntna head profiles
        id2=ytna_array(k); %indexing for ytna head profiles
        for ii=1:50 % Looping over x values in domain
            for iii=1:50 % Looping over y values in domina
                    Xo=Xc(ii);  %Position of measurement X
                    Yo=Yc(iii); %Position of measuremnt Y
                    Zo=base_top_elev(ii,iii);    %Position of measurment Z
                  for i=1:50 %Looping over domain for gravity of each X point
                     for j=1:50 %Looping over domain for gravity of each Y point
                    %Uncomment following for Point Mass Approximation if applicable and comment out one of the X and Y loops
%                     Zcom=(Headprofilesmocs(:,:,id)+Headprofilesmocs(:,:,id2))./2;      
%                     dm=1000*1000*(Headprofilesmocs(:,:,id)-Headprofilesmocs(:,:,id2))*(p)*d; %x,y,z dimensions of each cell
%                     Zoo=repmat(Zo,50,50);
%                     Xoo=repmat(Xo,50,50);
%                     Yoo=repmat(Yo,50,50);
%                     Xco=repmat(Xc,50,1);
%                     Yco=repmat(Yc,50,1)';
%                     r=sqrt(((Xoo-Xco).^2)+((Yoo-Yco).^2)+((Zoo-Zcom).^2));
%                     gz=(sqrt((Zoo-Zcom).^2).*gam.*dm)./(r.^3);
%                     gzn=reshape(gz,1,2500);
                    %Gravity Response Of a Rectanular Prism
                        dx1=(Xc1(i)-Xo); % X distance from measurement location to  X1
                        dx2=(Xc2(i)-Xo); % X distance from measurement location to  X2
                        dy1=(Yc1(j)-Yo); % Y distance from measurement location to  Y1
                        dy2=(Yc2(j)-Yo); % X distance from measurement location to  Y2
                        dz1=(Zo-Headprofilesmocs(i,j,id)); % Z distance from measurement location to  Z1
                        dz2=(Zo-Headprofilesmocs(i,j,id2)); % Z distance from measurement location to  Z2                        R111=sqrt(dx1.^2+dy1.^2+dz1.^2);
                        R112=sqrt(dx2.^2+dy1.^2+dz1.^2); % R distance to point on prism
                        R121=sqrt(dx1.^2+dy2.^2+dz1.^2); % R distance to point on prism
                        R122=sqrt(dx2.^2+dy2.^2+dz1.^2); % R distance to point on prism
                        R211=sqrt(dx1.^2+dy1.^2+dz2.^2); % R distance to point on prism
                        R212=sqrt(dx2.^2+dy1.^2+dz2.^2); % R distance to point on prism
                        R221=sqrt(dx1.^2+dy2.^2+dz2.^2); % R distance to point on prism
                        R222=sqrt(dx2.^2+dy2.^2+dz2.^2); % R distance to point on prism

                        g111=-[dz1.*atan((dx1.*dy1)./(dz1.*R111))-dx1.*log(R111+dy1)-dy1.*log(R111+dx1)]; %Z component due to gravity without constants incorporated at point
                        g112=+[dz1.*atan((dx2.*dy1)./(dz1.*R112))-dx2.*log(R112+dy1)-dy1.*log(R112+dx2)]; %Z component due to gravity without constants incorporated at point
                        g121=+[dz1.*atan((dx1.*dy2)./(dz1.*R121))-dx1.*log(R121+dy2)-dy2.*log(R121+dx1)]; %Z component due to gravity without constants incorporated at point
                        g122=-[dz1.*atan((dx2.*dy2)./(dz1.*R122))-dx2.*log(R122+dy2)-dy2.*log(R122+dx2)]; %Z component due to gravity without constants incorporated at point

                        g211=+[dz2.*atan((dx1.*dy1)./(dz2.*R211))-dx1.*log(R211+dy1)-dy1.*log(R211+dx1)]; %Z component due to gravity without constants incorporated at point
                        g212=-[dz2.*atan((dx2.*dy1)./(dz2.*R212))-dx2.*log(R212+dy1)-dy1.*log(R212+dx2)]; %Z component due to gravity without constants incorporated at point
                        g221=-[dz2.*atan((dx1.*dy2)./(dz2.*R221))-dx1.*log(R221+dy2)-dy2.*log(R221+dx1)]; %Z component due to gravity without constants incorporated at point
                        g222=+[dz2.*atan((dx2.*dy2)./(dz2.*R222))-dx2.*log(R222+dy2)-dy2.*log(R222+dx2)]; %Z component due to gravity without constants incorporated at point 
                        gz(i,j)=drho.*gam.*(g111+g112+g121+g122+g211+g212+g221+g222);  %Gravity change due to each cell at a defined measurement location
                    end 
                  end
                  gzn=reshape(gz,1,2500); 
                  gzt(ii,iii,k)=sum(gzn);  %Gravity change at a single measurement location 
            end
        end
          
                 
               
                            
               % Uncomment following equations for water level change at each location         
%                gzt(:,:,k)=Headprofilesmocs(:,:,id)-Headprofilesmocs(:,:,id2); %water level change between ntna and ytna
%              % Headprofilesmocs will have a 2 added to the end to specify MOCS
    end

 %% Data for MOCs(GZT of MOCs)
    Xc=[500:1000:49500];Yc=[500:1000:49500];
    Xc1=[0:1000:49000];Xc2=[1000:1000:50000];
    Yc1=[0:1000:49500];Yc2=[1000:1000:50000];
    % Setting Up Loops for Condtitions
    m=q;
    mntna_array=[1:m];
    mntna_array=mntna_array*3-2;
    mytna_array=[1:m];
    mytna_array=(mytna_array*3-2+1);
    mytya_array=[1:m];
    mytya_array=(mytya_array*3-2+2);
    % Model Run Loop

    for k=1:q % Same process as comments describe above but only for MOCs
    id=mntna_array(k);
    id2=mytna_array(k);
    for ii=1:50
            for iii=1:50
                    Xo=Xc(ii); 
                    Yo=Yc(iii);
                    Zo=base_top_elev(ii,iii);    
                  for i=1:50
                     for j=1:50
                        dx1=(Xc1(i)-Xo);
                        dx2=(Xc2(i)-Xo);
                        dy1=(Yc1(j)-Yo);
                        dy2=(Yc2(j)-Yo);
                        dz1=(Zo-Headprofilesmocs2(i,j,id));
                        dz2=(Zo-Headprofilesmocs2(i,j,id2));
                        R111=sqrt(dx1.^2+dy1.^2+dz1.^2);
                        R112=sqrt(dx2.^2+dy1.^2+dz1.^2);
                        R121=sqrt(dx1.^2+dy2.^2+dz1.^2);
                        R122=sqrt(dx2.^2+dy2.^2+dz1.^2);
                        R211=sqrt(dx1.^2+dy1.^2+dz2.^2);
                        R212=sqrt(dx2.^2+dy1.^2+dz2.^2);
                        R221=sqrt(dx1.^2+dy2.^2+dz2.^2);
                        R222=sqrt(dx2.^2+dy2.^2+dz2.^2);

                        g111=-[dz1.*atan((dx1.*dy1)./(dz1.*R111))-dx1.*log(R111+dy1)-dy1.*log(R111+dx1)];
                        g112=+[dz1.*atan((dx2.*dy1)./(dz1.*R112))-dx2.*log(R112+dy1)-dy1.*log(R112+dx2)];
                        g121=+[dz1.*atan((dx1.*dy2)./(dz1.*R121))-dx1.*log(R121+dy2)-dy2.*log(R121+dx1)];
                        g122=-[dz1.*atan((dx2.*dy2)./(dz1.*R122))-dx2.*log(R122+dy2)-dy2.*log(R122+dx2)];

                        g211=+[dz2.*atan((dx1.*dy1)./(dz2.*R211))-dx1.*log(R211+dy1)-dy1.*log(R211+dx1)];
                        g212=-[dz2.*atan((dx2.*dy1)./(dz2.*R212))-dx2.*log(R212+dy1)-dy1.*log(R212+dx2)];
                        g221=-[dz2.*atan((dx1.*dy2)./(dz2.*R221))-dx1.*log(R221+dy2)-dy2.*log(R221+dx1)];
                        g222=+[dz2.*atan((dx2.*dy2)./(dz2.*R222))-dx2.*log(R222+dy2)-dy2.*log(R222+dx2)];
                        gz(i,j)=drho.*gam.*(g111+g112+g121+g122+g211+g212+g221+g222);  
                    end 
                  end
                  gzn=reshape(gz,1,2500);
                  gztm(ii,iii,k)=sum(gzn);  
                %Point Mass Approximation
%     for i=1:50
%         for j=1:50
%             Xo=Xc(i);
%             Yo=Yc(j);
%             Zo=base_top_elev(i,j);
%             Zcom=(Headprofilesmocs2(:,:,id)+Headprofilesmocs2(:,:,id2))/2;
%             dm=1000*1000*(Headprofilesmocs2(:,:,id)-Headprofilesmocs2(:,:,id2))*(p)*d;
%             Zoo=repmat(Zo,50,50);
%             Xoo=repmat(Xo,50,50);
%             Yoo=repmat(Yo,50,50);
%             Xco=repmat(Xc,50,1);
%             Yco=repmat(Yc,50,1)';
%             r=sqrt(((Xoo-Xco).^2)+((Yoo-Yco).^2)+((Zoo-Zcom).^2));
%             gz=(sqrt((Zoo-Zcom).^2).*gam.*dm)./(r.^3);
%             gzn=reshape(gz,1,2500); 
%             gztm(j,i,k)=sum(gzn);   
            % Water Level Change for MOCs
%            gztm(:,:,k)=Headprofilesmocs2(:,:,id)-Headprofilesmocs2(:,:,id2); %water level change
        end
    end

    end
     gztfinal=cat(3,gzt,gztm);       %Matrix of gravity change or water level change values measured at every location in the model ensemble( Used to create likelihoods)

%% Truth Model Selection(Model that the data is being compared to)
t=q; % Selection of Subset of models (non-MOC-s,MOC-q) (remember to switch the index of Headprofilesmocs varaible for Non-MOC models)
tntna_array=[1:t];
tntna_array=tntna_array*3-2;
tytna_array=[1:t];
tytna_array=(tytna_array*3-2+1);
tytya_array=[1:t];
tytya_array=(tytya_array*3-2+2);
Xc=[500:1000:49500];Yc=[500:1000:49500];
 Xc1=[0:1000:49000];Xc2=[1000:1000:50000];
 Yc1=[0:1000:49500];Yc2=[1000:1000:50000];
for tmm=1:q % Same process as comments describe above until gzt_ft
        tm1=tntna_array(tmm);
        tm2=tytna_array(tmm);
        tm3=tytya_array(tmm);
         for ii=1:50
            for iii=1:50
                    Xo=Xc(ii); 
                    Yo=Yc(iii);
                    Zo=base_top_elev(ii,iii);    
                  for i=1:50
                     for j=1:50
                        dx1=(Xc1(i)-Xo);
                        dx2=(Xc2(i)-Xo);
                        dy1=(Yc1(j)-Yo);
                        dy2=(Yc2(j)-Yo);
                        dz1=(Zo-Headprofilesmocs2(i,j,id2));
                        dz2=(Zo-Headprofilesmocs2(i,j,id));
                        R111=sqrt(dx1.^2+dy1.^2+dz1.^2);
                        R112=sqrt(dx2.^2+dy1.^2+dz1.^2);
                        R121=sqrt(dx1.^2+dy2.^2+dz1.^2);
                        R122=sqrt(dx2.^2+dy2.^2+dz1.^2);
                        R211=sqrt(dx1.^2+dy1.^2+dz2.^2);
                        R212=sqrt(dx2.^2+dy1.^2+dz2.^2);
                        R221=sqrt(dx1.^2+dy2.^2+dz2.^2);
                        R222=sqrt(dx2.^2+dy2.^2+dz2.^2);

                        g111=-[dz1.*atan((dx1.*dy1)./(dz1.*R111))-dx1.*log(R111+dy1)-dy1.*log(R111+dx1)];
                        g112=+[dz1.*atan((dx2.*dy1)./(dz1.*R112))-dx2.*log(R112+dy1)-dy1.*log(R112+dx2)];
                        g121=+[dz1.*atan((dx1.*dy2)./(dz1.*R121))-dx1.*log(R121+dy2)-dy2.*log(R121+dx1)];
                        g122=-[dz1.*atan((dx2.*dy2)./(dz1.*R122))-dx2.*log(R122+dy2)-dy2.*log(R122+dx2)];

                        g211=+[dz2.*atan((dx1.*dy1)./(dz2.*R211))-dx1.*log(R211+dy1)-dy1.*log(R211+dx1)];
                        g212=-[dz2.*atan((dx2.*dy1)./(dz2.*R212))-dx2.*log(R212+dy1)-dy1.*log(R212+dx2)];
                        g221=-[dz2.*atan((dx1.*dy2)./(dz2.*R221))-dx1.*log(R221+dy2)-dy2.*log(R221+dx1)];
                        g222=+[dz2.*atan((dx2.*dy2)./(dz2.*R222))-dx2.*log(R222+dy2)-dy2.*log(R222+dx2)];
                        gz(i,j)=drho.*gam.*(g111+g112+g121+g122+g211+g212+g221+g222);  
                    end 
                  end
                  gzn=reshape(gz,1,2500);
                  gzt_truth(ii,iii)=sum(gzn);  
            end
         end
         
            %Point Mass Approximation
%         for i=1:50
%             for j=1:50
%                 Xo=Xc(i);
%                 Yo=Yc(j);
%                 Zo=base_top_elev(i,j);
%                 Zcom=(Headprofilesmocs(:,:,tm1)+Headprofilesmocs(:,:,tm2))./2;   % Models Index within folder
%                 dm=(1000*1000).*(Headprofilesmocs(:,:,tm1)-Headprofilesmocs(:,:,tm2)).*(p).*d;  % Models Index within folder
%                 Zoo=repmat(Zo,50,50);
%                 Xoo=repmat(Xo,50,50);
%                 Yoo=repmat(Yo,50,50);
%                 Xco=repmat(Xc,50,1);
%                 Yco=repmat(Yc,50,1)';
%                 r=sqrt(((Xoo-Xco).^2)+((Yoo-Yco).^2)+((Zoo-Zcom).^2));
%                 gz=(sqrt((Zoo-Zcom).^2).*gam.*dm)./(r.^3);
%                 gzn=reshape(gz,1,2500);
% %                 gzn(gzn <0.015) = 0;
%                 gzt_truth(j,i)=sum(gzn);
                %Comment out above and uncomment below to use water levels
% %                   gzt_truth=Headprofilesmocs2(:,:,tm1)-Headprofilesmocs2(:,:,tm2);  %Water Level Change
%             end
%         end
        TruePOI1=Headprofilesmocs2(22,39,tm2)-Headprofilesmocs2(22,39,tm3); %True POI of truth models( delete the 2 for non-moc models)
        u=tmm; %Model Number
 % GZT of Other models

    gzt_ft=repmat(gzt_truth,1,1,w); 
    for i=1:w
        m1 = gztfinal(:,:,i);
        MisM=((m1-gzt_ft(:,:,i)).^2); %Mismatch calculation of each model with the truth model
        m2{i} = MisM;     
    end
    MisM2 = cat(3,m2{:}); %Mismatch of each models cells
    MisM2(MisM2==0) = NaN;  %Elimatinating the truth model from ensemble

    for a=1:50
        for b=1:50
          MISM=MisM2(a,b,:);    %(Point location (#(1-50),#(1-50),set of models)
          MISM=reshape(MISM,1,w); 
          MISMofPoint=sum(1./MISM,'omitnan'); % Setting the common metric of each point
          LikeH(a,b,:)=((1./MISM)./(MISMofPoint)); % Calculation of Likelihood of each model at each point
          CheckSum(a,b)=sum(LikeH(a,b,:),3,'omitnan'); % Insurance metric to check that likelihoods add to 1
        end
    end
    LikeGravityPrismotherSYM3fixedmoc(:,:,:,tmm)=LikeH; % Save likelihoods for each truth model run 
%% POI Predictions 
for a=1:50
    for b=1:50
       cc1= LikeH(a,b,1:w); % Likelihood of current truth model at (x,y) location
       cc= squeeze(cc1); % shrink the matrix so likelohoods can weight predicted values of each model
       bb= APOI.*cc; % Applied Weighting of each POI with likelihood 
       POIweighted(a,b)=sum(bb,'omitnan'); % Summed to a single value for each (x,y) measurement location for set truth model
    end
end

POIGravityPrismotherSYPredictionM3fixedmoc(:,:,tmm)=POIweighted; % Saved Likelihood POI values for each truth model

%% Likelihood Weighted POI Error
TruePOI2=repmat(TruePOI1,50,50);        %Truth POI of Truth model (replicated over domain)

POIerrorGravityPrismSYPredictionM3fixedmoc(:,:,tmm)=abs(POIGravityPrismotherSYPredictionM3fixedmoc(:,:,tmm)-TruePOI2); %POI error from truth model

 end

