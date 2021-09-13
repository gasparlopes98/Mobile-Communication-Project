clc
clear all, close all

viewer = siteviewer("Buildings","\leiriafinal.osm");
 
%======= Frequências================
freq_N= 1830e6; %[1830e6, 2620e6];
freq_N2=2620e6;
freq_SE=791e6;  %[1855e6,791e6,2145e6];
freq_SE2=1855e6;
freq_SE3=2145e6;
freq_SO= 1805e6; %[1805e6,2110e6];
freq_SO2 = 2110e6;


%=== Site 1 === 
site1.lats = 39.74314151310942;
site1.lons = -8.809610052851166;
site1.freq = [freq_N, freq_SO,freq_SE2]; % 
ant1_N = design(reflectorSpherical,freq_N);
ant1_N.Tilt = [-92 90]; 
ant1_N.TiltAxis = [1 0 0; 0 1 0];
ant1_SO = design(reflectorSpherical,freq_SO);
ant1_SO.Tilt = [-92 90 120]; %10
ant1_SO.TiltAxis = [1 0 0;0 1 0; 0 0 1]; % ; 0 0 1
ant1_SE2 = design(reflectorSpherical,freq_SE2);
ant1_SE2.Tilt = [-92 90 -120]; %10
ant1_SE2.TiltAxis = [1 0 0;0 1 0; 0 0 1];
site1.antena = [ant1_N,ant1_SO,ant1_SE2]; % 

 
%=== Site 2 === 
site2.lats = 39.7423600;
site2.lons = -8.8066526;
site2.freq = [freq_N,freq_N2]; % 800e6, 810e6
ant2_N = design(reflectorSpherical,freq_N);
ant2_N.Tilt = [-91 90];
ant2_N.TiltAxis = [1 0 0; 0 1 0]; 
ant2_N2 = design(reflectorSpherical,freq_N2);
ant2_N2.Tilt = [-91 90];
ant2_N2.TiltAxis = [1 0 0; 0 1 0]; 
site2.antena = [ant2_N,ant2_N2];

 

%=== Site 3 === 
site3.lats = 39.748609;
site3.lons = -8.804629;
site3.freq = [freq_SE,freq_SE2,freq_SE3,freq_SO,freq_SO2]; %791e6,815e6
ant3_SE = design(reflectorSpherical,freq_SE);
ant3_SE.Tilt = [-92 90 -120]; 
ant3_SE.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
ant3_SE2 = design(reflectorSpherical,freq_SE2);
ant3_SE2.Tilt = [-92 90 -120]; 
ant3_SE2.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
ant3_SE3 = design(reflectorSpherical,freq_SE3);
ant3_SE3.Tilt = [-92 90 -120];
ant3_SE3.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
ant3_SO = design(reflectorSpherical,freq_SO);
ant3_SO.Tilt = [-92 90 120];
ant3_SO.TiltAxis = [1 0 0; 0 1 0; 0 0 1]; 
ant3_SO2 = design(reflectorSpherical,freq_SO2);
ant3_SO2.Tilt = [-92 90 120];
ant3_SO2.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
site3.antena = [ant3_SE,ant3_SE2,ant3_SE3,ant3_SO,ant3_SO2];

 

%=== Site 4 === 
site4.lats = 39.744955;
site4.lons = -8.811599;
site4.freq = [freq_SE,freq_SE2,freq_SE3];
ant4_SE = design(reflectorSpherical,freq_SE);
ant4_SE.Tilt = [-92 90 -120]; 
ant4_SE.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
ant4_SE2 = design(reflectorSpherical,freq_SE2);
ant4_SE2.Tilt = [-92 90 -120]; 
ant4_SE2.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
ant4_SE3 = design(reflectorSpherical,freq_SE3);
ant4_SE3.Tilt = [-92 90 -120];
ant4_SE3.TiltAxis = [1 0 0; 0 1 0; 0 0 1];
site4.antena = [ant4_SE,ant4_SE2,ant4_SE3];

 

%=== Site 5 === 
site5.lats = 39.74697956695493;
site5.lons = -8.807343960520848;
site5.freq = [freq_N2, freq_SO2,freq_SE,freq_SE3];
ant5_N2 = design(reflectorSpherical,freq_N2);
ant5_N2.Tilt = [-92 90]; 
ant5_N2.TiltAxis = [1 0 0; 0 1 0];
ant5_SO2 = design(reflectorSpherical,freq_SO2);
ant5_SO2.Tilt = [-92 90 120]; %10
ant5_SO2.TiltAxis = [1 0 0;0 1 0; 0 0 1]; % ; 0 0 1
ant5_SE = design(reflectorSpherical,freq_SE);
ant5_SE.Tilt = [-92 90 -120]; %10
ant5_SE.TiltAxis = [1 0 0;0 1 0; 0 0 1];
ant5_SE3 = design(reflectorSpherical,freq_SE3);
ant5_SE3.Tilt = [-92 90 -120]; %10
ant5_SE3.TiltAxis = [1 0 0;0 1 0; 0 0 1];
site5.antena = [ant5_N2,ant5_SO2,ant5_SE,ant5_SE3];


tx1= txsite('Name','Site 1',...
               'Latitude', site1.lats,...
               'Longitude', site1.lons,...
                 'Antenna', site1.antena,...
           "AntennaHeight", 22,...
    "TransmitterFrequency", site1.freq,...
        "TransmitterPower",50,...
        'SystemLoss',2);  
tx2= txsite('Name','Site 2',...
               'Latitude', site2.lats,...
               'Longitude', site2.lons,...
                 'Antenna', site2.antena,...
           "AntennaHeight", 31,...
    "TransmitterFrequency", site2.freq,...
        "TransmitterPower",50,...
        'SystemLoss',2);
tx3= txsite('Name','Site 3',...
                'Latitude', site3.lats,...
               'Longitude', site3.lons,...
                 'Antenna', site3.antena,...
           "AntennaHeight", 25.5,...
    "TransmitterFrequency", site3.freq,...
        "TransmitterPower",50,...
        'SystemLoss',2);       
tx4= txsite('Name','Site 4',...
              'Latitude', site4.lats,...
               'Longitude', site4.lons,...
                 'Antenna', site4.antena,...
           "AntennaHeight", 20,...
    "TransmitterFrequency", site4.freq,...
        "TransmitterPower",50,...
        'SystemLoss',2);    
tx5= txsite('Name','Site 5',...
              'Latitude', site5.lats,...
               'Longitude', site5.lons,...
                 'Antenna', site5.antena,...
           "AntennaHeight", 30,...
    "TransmitterFrequency", site5.freq,...
        "TransmitterPower",50,...
        'SystemLoss',2);   

txs=[tx1, tx2, tx3, tx4, tx5];
%show(txs)

%========== Cobertura Longley-Rice ============
% 
% coverage(txs,'SignalStrengths',-100:5:-5,...
%         'MaxRange', 750,...
%         'Resolution', 3)

%========== SINR =============== 
%sinr(txs)

%========== Padrões das Antenas ===========
pattern(tx1(1),site1.freq(1))  
pattern(tx1(2),site1.freq(2)) 
pattern(tx1(3),site1.freq(3))

pattern(tx2(1),site2.freq(1)) 
pattern(tx2(2),site2.freq(2)) 

pattern(tx3(1),site3.freq(1)) 
pattern(tx3(2),site3.freq(2)) 
pattern(tx3(3),site3.freq(3)) 
pattern(tx3(4),site3.freq(4)) 
pattern(tx3(5),site3.freq(5)) 

pattern(tx4(1),site4.freq(1)) 
pattern(tx4(2),site4.freq(2))
pattern(tx4(3),site4.freq(3)) 

pattern(tx5(1),site5.freq(1))
pattern(tx5(2),site5.freq(2))
pattern(tx5(3),site5.freq(3))
pattern(tx5(4),site5.freq(4))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulações RXS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%viewer = siteviewer("Buildings","C:\Users\Mateo\Documents\IPL\6S\mov\leiriafinal.osm");

%=============== BLOCO =====================================================

inilat = 39.74392260879267;
nextlat = inilat;
inilon = -8.811286734948416;
names = ["ini"];
lats = [inilat];
lons = [inilon];
%1 lat = 111.32 km -> 111320m
%1 lon = 40075 km * cos( latitude ) / 360
for i = 1:10
        for j = 1:15
            
            if i == 1
                nextlat = 0.0004 + lats(end);
            lats(i,j) = nextlat;
            nextlon = 0.0004 + lons(end);
            lons(i,j) = nextlon;
            
            end
            names(i,j) = string(i)+string(j);
        end
    %nextlat = -0.0001 + lats(end);
    lons(i,:) = lons(1,:);
    %lats(i) = nextlat;
    lats(i,:) = lats(1,:) -0.0005*(i-1);
end

lats = reshape(lats,[1,15*10]);
lons = reshape(lons,[1,15*10]);
names = reshape(names,[1,15*10]);

%lats = [39.74452169295092];
%lons = [-8.808020439303462];

%=============== RUA =====================================================

% for i =0:25
%     names(i+1) = string(i*10);
% end
% 
% 
% lats =   [39.744623363789, 39.74470690000001, 39.74478733171882,...
%         39.744863638647274, 39.74493946052537, 39.74502195431312,...
%         39.745104366660065, 39.745182219988564, 39.74526007323098,...
%         39.74532864588302, 39.74539412496557, 39.74545651049144,...
%         39.745517864796106, 39.74558128137187, 39.745639026484504,...
%         39.7456921313345, 39.74574523614322, 39.74579974812201,...
%         39.74586007096292, 39.74593018975343, 39.74599824615689,...
%         39.74606836480727, 39.74613848338624, 39.74620138381538,...
%         39.74626479976203, 39.7463313091068];
%     
% lons = [-8.809767302444596, -8.809719722715101, -8.809667419638565,...
%       -8.809602376069185, -8.809543353896194, -8.809497085792202,...
%       -8.80945247252369, -8.809393463927295, -8.809335125879151,...
%       -8.809259353470146, -8.809178887202888, -8.809093727065939,...
%       -8.809010578587992, -8.808925418454562, -8.808837576120515,...
%       -8.808743028251861, -8.808647809832907, -8.808555797717258,...
%       -8.808469296475998, -8.808394194624109, -8.808319092774155,...
%       -8.808246002575547, -8.808172241827481, -8.808087752243557,...
%       -8.808008627080213, -8.807928831364636];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rxs = rxsite('Name', names,...
      'Antenna',pifa, 'Latitude',lats,...
       'Longitude',lons);
show(rxs)

pm = propagationModel("raytracing", ...
                        'Method', 'sbr');
pm.MaxNumReflections = 4;

[~,numrxs] = size(rxs);

num = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
allpowrt(:,1) = num;
allpowlr(:,1) = num;
bestpow(1:3,1:numrxs) = 0;
for i = 1:numrxs
    ssPerfect = sigstrength(rxs(i),txs,pm);
    allpowrt(:,i+1) = ssPerfect;
    ssPerfectlr = sigstrength(rxs(i),txs);
    allpowlr(:,i+1) = ssPerfectlr;
    [maxpowrt,indrt] = max(ssPerfect);
    [maxpowlr,indlr] = max(ssPerfectlr);
    if maxpowrt >= maxpowlr
        maxpow = maxpowrt;
        ind = indrt;
        method = 1;%'ray-tracing';
    else
        maxpow = maxpowlr;
        ind = indlr;
        method = 2;%'longley-rice';
    end
    bestpow(:,i) = [maxpow;ind;method];
end
%ssPerfect = sigstrength(rxs,txs,pm);

raytrace(txs, rxs, pm) 
