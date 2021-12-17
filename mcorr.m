function mcorr(rootxml,roottemp,DAYDIR)

% Multichannel xcorr 
% dipendenze: ReadMSEEDFast.m
%             findwav_mseed.m
%
% Legge QuAKEML catalog e i relativi files mseed 
% Estrae i template e fa la crosscorrelazione
% Scrive i risultati di ogni canale/fase nella dir WORK 
% Al termine crea una struttura matlab nella dir OUT
% GLI EVENTI DEVONO AVERE PICK P e S 
%
% per preparare i templates usare il codice python alla fine del codice
%
% -----------------------------------------
% usage:    mcorr(rootxml,roottemp,DAYDIR)
% -----------------------------------------
%
% Davide Piccinini Dec 2020 mailto:davide.piccinini@ingv.it
% mcorr ver. 1.1

% 04/02/2021
% + Velocizzata l'associazione degli eventi
% + cambiato il modo di lettura dei dati mseed (da ReadMSEEDFast a rdmseed)
% + incluse tutte le funzioni associate (findwav e rdmseed)
% 17/02/2021
% + aggiunto in findwav_mseed il controllo sulla maxCC dei doppioni
% + aggiunge il valore di prewin ai tempi degli slaves
% 18/02/2021
% + aggiunto un controllo sulla PREF e SREF 
% 17/12/2021
% + upgraded decodeQML to correctly get only real picks


%% SETUP WORKING DIR
if exist('WORK','dir')==0
!mkdir WORK
!mkdir OUT
end

%% SETUP PARAMETERS
Plwin=  0.5; % sec of P
Slwin=  0.8; % sec of S
prewin= 0.1;  % add this before P and S

CORRE=1;  % time correction for Plwin sets Plwin eq to S-P 

NPHP=4;   % Numero minimo di fasi P
NPHS=2;   % numero minimo di fasi S

FILT=[3 15]; % estremi del filtro

THR=0.65;  % soglia di CC

GARBAGE=1;% clear WORK directory each template
% warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LEGGE L'EVENTO
%KK=textread([rootxml '.xml'],'%s','delimiter','\n');
fprintf('Getting template picks from %s\n',rootxml)
filexml=[rootxml '.xml'];
X=decodeQML(filexml);
%keyboard
%% esce se il n di fasi P ed S non è sufficiente
namex=fieldnames(X);
np=0; ns=0;
for k=1:length(namex)-1
    XU=getfield(X,char(namex(k)));
    if XU.P > 0; np=np+1; end
    if XU.S > 0; ns=ns+1; end
end
if np < NPHP || ns < NPHS ; fprintf('Template has insufficient number of P and S [%d %d / %d %d]\n',np,ns,NPHP,NPHS);return; end 
%% LEGGE IL TEMPLATE MSEED
D=dir([roottemp '.*.ms']);
if isempty(D)==1
    return
end

nn=0;
for ii=1:length(D)
        nn=nn+1;
        S(nn)=ReadMSEED([D(ii).folder '/' D(ii).name ]);
end
% keyboard
for j=1:length(S)
    S(j).station=strrep(S(j).station,' ',''); %% accounting for 3char station code
    trace=S(j).data;
    trace=trace-mean(trace);
    trace=tapera(trace,50);
    S(j).datafilt = filtrax(trace,FILT(1),FILT(2),S(j).sampleRate);
end

n=0;
FN=fieldnames(X);
for k=1:length(S)
    %pause
    NM=S(k).station;
    NM=strrep(NM,' ',''); %% accounting for 3char station code
    for jj=1:length(FN)-1
        %fprintf('%s %s\n',NM,char(FN(jj)));
        if strcmp(char(FN(jj)),NM)==1
            %fprintf('%s %s\n',NM,char(FN(jj)));
            PT=X.(NM).P;
            Punc=X.(NM).Punc;
            ST=X.(NM).S;
            Sunc=X.(NM).Sunc;
            if CORRE==1 && isempty (PT) ==0 && isempty (ST) ==0
                %% check P win and resize to SP time
                SPT=86400*(ST-PT);
                if Plwin > SPT
                    Plwin=SPT;
                    fprintf('Resizing Pwin length @%s to %4.2f s\n',NM,Plwin)
                end
            end
            if isempty (PT) == 0
                n=n+1;
                [~, ps]=min(abs(S(k).matlabTimeVector-PT));
                S(k).PTIME=PT;
                S(k).Punc=Punc;               
                try
                    S(k).ptemplate=(S(k).datafilt(round(ps-(prewin*S(k).sampleRate)):round(ps+(Plwin*S(k).sampleRate))));
                    if  ST > 0
                        [~, ss]=min(abs(S(k).matlabTimeVector-ST));
                        S(k).STIME=ST;
                        S(k).Sunc=Sunc;
                        S(k).stemplate=(S(k).datafilt(ss-(prewin*S(k).sampleRate):ss+(Slwin*S(k).sampleRate)));
                    else
                        S(k).STIME=[];
                        S(k).Sunc =[];
                        S(k).stemplate=[];
                    end                
                NS(n)=S(k);
                catch
                    fprintf('Something goes wrong with the pick [Maybe S outside the template?] @%s\n',NM)
                    n=n-1;
                end
            end
        end
    end
end
%keyboard
%% CHECK AGAIN P AND S NUMBER
np=0;
ns=0;
for k=1:length(NS)
    if NS(k).PTIME > 0
        np=np+1;
    end
    if NS(k).STIME > 0
        ns=ns+1;
    end
end
if np < NPHP || ns < NPHS; fprintf('Insufficient phases --- STOP ---\n'); return; end

% keyboard

for k=1:length(NS)   %% ADDING ORIGIN TIME
    NS(k).OTime=char(X.HEADER(1));
    NS(k).MOTime=cell2mat(X.HEADER(2));
    NS(k).LAT=cell2mat(X.HEADER(3));
    NS(k).LON=cell2mat(X.HEADER(4));
    NS(k).LON=cell2mat(X.HEADER(5));
end
OTime=NS(1).OTime;
% MOTime=NS(1).MOTime;

fprintf('Number of stations : %3.0f\n',length(NS))

T0=clock;
%keyboard
%% %%%% XCORR CORE CALL %%%%%%%%
parfor k=1:length(NS)
%for k=1:length(NS)   % uncomment for non parallel
    STA=NS(k).station;
    d=dir([DAYDIR '/*' STA '.*.*HZ*']);
    if length(d) > 0
        DAY=ReadMSEED([d(1).folder '/' d(1).name]);
        trace=DAY.data;
        if length(trace) > 1000 %% skip dayfiles with no data
        trace=trace-mean(trace);
        trace=tapera(trace,50);
        DAY.datafilt = filtrax(trace,FILT(1),FILT(2),DAY.sampleRate);
        fileday=DAY;
        findwav_mseed(fileday,NS(k),THR,'P');
            if isempty(getfield(NS(k),'STIME'))==0
                findwav_mseed(fileday,NS(k),THR,'S');
            end
        end
    else
        fprintf('No data in dayfiles for station %s\n',STA)
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ET=etime(clock,T0);
fprintf('Total time used for correlation: %8.2f sec\n',ET)
%% %%%%%%%%%%%%%%%%%%%%
%% ASSOCIA LE P E LE S
%% ====================
T0=clock;
%tic

T = struct2table(NS);
SS=sortrows(T,'PTIME');
NS=table2struct(SS);

% keyboard


%parfor k=1:length(NS)
for k=1:length(NS);  % uncomment for non parallel
    STA=NS(k).station;
    PTEMP(k)=NS(k).PTIME;
    try
    Dp=dir(['WORK/' OTime '.' STA '*.P.' DAYDIR '.out']);
    fprintf('Reading %s\n',['WORK/' OTime '.' STA '*.P.' DAYDIR '.out'])

    [PT, CORR(k).AP, CORR(k).CCP]=textread(['WORK/' Dp(1).name],'%s %f %f');
    CORR(k).PT=datenum(PT,'yyyy-mm-ddTHH:MM:SS.FFF');
    catch ME
    fprintf('No P time for Station %s\n',STA);
    CORR(k).AP=NaN;
    CORR(k).CCP=NaN;
    CORR(k).PT=NaN;
    PTEMP(k)=NaN;
    end
    
    try
    STEMP(k)=NS(k).STIME;
    Ds=dir(['WORK/' OTime '.' STA '*.S.' DAYDIR '.out']);
    fprintf('Reading %s\n',['WORK/' OTime '.' STA '*.S.' DAYDIR '.out'])
    [ST, CORR(k).AS, CORR(k).CCS]=textread(['WORK/' Ds(1).name],'%s %f %f');
    CORR(k).ST=datenum(ST,'yyyy-mm-ddTHH:MM:SS.FFF');  
    catch ME
    fprintf('No S time for Station %s\n',STA);
    STEMP(k)=NaN;
    CORR(k).AS=NaN;
    CORR(k).CCS=NaN;
    CORR(k).ST=NaN;
    end    
end

%keyboard
%% STABILISCE LA STAZIONE DI RIFERIMENTO PER IL CALCOLO DELLE TT
SEC=1/86400;
PREF=1; %STAZIONE DI RIFERIMENTO PER LE P E' LA PRIMA STAZIONE
if isnan(CORR(PREF).CCP)==1
    for k=1:length(CORR)
        n(k)=length(CORR(k).CCP);
    end
    is=find(n>1);
    PREF=min(is);
end

%% SREF POTREBBE NON ESSERE LA 1° STAZ PERCHE' POTREBBE NON AVERE LA P
for k=1:length(NS)
    if isempty(NS(k).STIME)
        SREF(k)=0;
    else
        SREF(k)=NS(k).STIME;
    end
end
is=find(SREF > 0);
SREF=min(is);  
for k=1:length(CORR)
       n(k)=length(CORR(k).CCS);
end
is=find(n > 1);
if min(is) ~= SREF
    SREF=min(is);
end

fprintf('REFERENCE STATION FOR P WAVES IS %s \n',NS(PREF).station);
fprintf('REFERENCE STATION FOR S WAVES IS %s \n',NS(SREF).station);
    
%keyboard
%%
TTP=(PTEMP-PTEMP(PREF))./SEC; % TEMPI DI ARRIVO P DEL TEMPLATE RELATIVI RISPETTO ALLA REF
TTS=(STEMP-STEMP(SREF))./SEC; % TEMPI DI ARRIVO S DEL TEMPLATE RELATIVI RISPETTO ALLA REF

for k=1:length(NS)
    CORR(k).PCHECK=CORR(k).PT-(TTP(k)*SEC);  %% toglie a tutti i match la travel time dell'evento per le P
    CORR(k).SCHECK=CORR(k).ST-(TTS(k)*SEC);  %% toglie a tutti i match la travel time dell'evento per le S
end
% keyboard
%% ALT CODE FOR EXTRACT EVENTS IMPROVED SPEED
%% Pwaves
fprintf('Collecting P picks\n')
THR=0.2*SEC;
clear POUT
clear ia ib
A=CORR(PREF).PCHECK;
for k=1:length(NS)
    B=CORR(k).PCHECK;
    [LIA, LIB]=ismembertol(A(1:end),B(1:end),THR/(max(abs([A(1:end); B(1:end)]))));
    ia(:,k)=(LIA);
    ib(:,k)=(LIB);
end
%%
S=sum(ia');
ii=find(S>=NPHP); %% CERCA TUTTI QUELLI CHE HANNO ALMENO NPHP FASI P
%%
for NEV=1:length(ii)
    n=0;
    TEMP=[];
    for nsta=1:length(NS)
        if (ib(ii(NEV),nsta))~=0
            n=n+1;
            TEMP.sta(n)={NS(nsta).station};    %% real station
            TEMP.Ptim(n)=CORR(nsta).PT(ib(ii(NEV),nsta));  %% real matched P time
            TEMP.Pcc(n) =CORR(nsta).CCP(ib(ii(NEV),nsta)); %% related cc
        end
    end
    POUT(NEV)=TEMP;
end

if isempty(NEV)==1; fprintf('No P-events found\n'); return; end  %% SE NON TROVA P ESCE
%% S waves
fprintf('Collecting S picks\n')
THR=0.3*SEC;
NEV=0;
clear SOUT
clear ia ib
A=CORR(SREF).SCHECK;
%keyboard
for k=1:length(NS)
    B=CORR(k).SCHECK;
    [LIA, LIB]=ismembertol(A(1:end),B(1:end),THR/(max(abs([A(1:end); B(1:end)]))));
    ia(:,k)=(LIA);
    ib(:,k)=(LIB);
end
S=sum(ia');
ii=find(S>=NPHS); %% CERCA TUTTI QUELLI CHE HANNO ALMENO NPHP FASI S

for NEV=1:length(ii)
    n=0;
    TEMP=[];
    for nsta=1:length(NS)
        if (ib(ii(NEV),nsta))~=0
            n=n+1;
            TEMP.sta(n)={NS(nsta).station};    %% real station
            TEMP.Stim(n)=CORR(nsta).ST(ib(ii(NEV),nsta));  %% real matched P time
            TEMP.Scc(n) =CORR(nsta).CCS(ib(ii(NEV),nsta)); %% related cc
        end
    end
    SOUT(NEV)=TEMP;
end
%keyboard
if isempty(NEV)==1; fprintf('No S-events found\n'); return; end

%keyboard
%%
%% MIX P and S matches according to template travel times  ||
fprintf('Combining %8.0d P-events & %8.0d S-events \n',length(POUT),length(SOUT))

%% alt code

SPTIMES=STEMP-PTEMP;  % S-P del template
for k=1:length(NS)
    Ts(k)={NS(k).station};
end
SRCH=0.2/86400;
%keyboard
% INIZIALIZZA SLAVE
SL.Psta=  [];
SL.Ptim=  [];
SL.Pcc=   [];
SL.Ssta=  [];
SL.Stim=  [];
SL.Scc=   [];
SL.TemplateOTime=[];
SLAVE=repmat(SL,1,1000); %% maximum number of matches per day

for uu=1:length(SOUT)
    Spipo(uu)=SOUT(uu).Stim(1);
end

% Offset is the prewin defined at the begin of the program added 17/02/2021
offset=prewin/86400;
%keyboard
L=0;
for j=1:length(POUT)
    Pt=POUT(j).Ptim;
    Ps=POUT(j).sta;
    [i, ~]=ismember(Ts,Ps); %% individua fra tutte solo le staz del possibile match
    STEO=Pt+SPTIMES(i);     %% aggiunge alle P del match S-P del template
        
    %% aggiunto il 4/2/2021 limita la ricerca delle S all'intervallo P -> Steorica+10 secondi
    MAXS=min(STEO)+(10/86400);
    % if MAXS > max(Spipo); MAXS=max(Spipo); end
    yy=find(Spipo > Pt(PREF) & Spipo < MAXS); % indici 
    SOUTTEMP=SOUT(yy);                                   % SOUTTEMP contiene solo le S degli indici
    %%
    for k=1:length(SOUTTEMP)
        St=SOUTTEMP(k).Stim;
        Ss=SOUTTEMP(k).sta;
        NumbS=ismember(Ps,Ss);  %% INDICI DELLE STAZIONI COINCIDENTI
        IS=ismembertol(STEO(NumbS),St,SRCH/(max(abs([STEO(NumbS) St(1:end)])))); % CERCA GLI ARRIVI S INTORNO ALL'ARRIVO DEL TEMPLATE
        if IS==1
            L=L+1;
            SLAVE(L).Psta= POUT(j).sta;
            SLAVE(L).Ptim= POUT(j).Ptim + offset;
            SLAVE(L).Pcc=  POUT(j).Pcc;
            SLAVE(L).Ssta= SOUTTEMP(k).sta;
            SLAVE(L).Stim= SOUTTEMP(k).Stim + offset;
            SLAVE(L).Scc=  SOUTTEMP(k).Scc;
            SLAVE(L).TemplateOTime=NS(1).MOTime;
        end
    end
end
SLAVE=SLAVE(1:L); %% RESIZE SLAVE

fprintf('===============================\n');
fprintf('# total matches found= %4.0f   \n',length(SLAVE));
fprintf('===============================\n');

ET=etime(clock,T0);
fprintf('Total time used for collect %4.0f Events: %8.2f sec \n',length(SLAVE),ET)

save(['OUT/' OTime '.' DAYDIR '.mat'],'SLAVE');

%% CLEAR OUT WORKFILES 
if GARBAGE==1
    fprintf('purging working directory...\n') 
    !rm WORK/*
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% END OF MAIN %%%%%

%% BEGIN OF FUNCTIONS


function c = filtrax(S,L,H,Fs);
%c = filtrax(S,L,H,Fs) Filtra un vettore S, usando un filtro di tipo Butterworth
%del 4 ordine.
%S=segnale in entrata
%H=limite superiore del filtro (Hz)
%L=limite inferiore del filtro (Hz)
%Fs=numero di campioni/secondo (cps)
%Davide Piccinini -Dicembre 1999-

if nargin < 4
    error('Usage: c=filtrax(S,L,H,Fs)');
end
if L & H ~=0
    [b,a]=butter(4,[L H]*2/Fs);
    %[H,w]=freqz(b,a,1024);
    c=filtfilt (b,a,S);
else
    c=S;
end

% 
function vout=tapera(vin,n)
%usage : vout=tapera(vin,n)
%Applica un taper di dimensioni n campioni
%agli estremi del vettore vin.
%Davide Piccinini - Aprile 2000


[r,c]=size(vin);

if r < c
   vin=vin';
end

%t=taper (n,1);
t=hanning(n*2);
t=t;
ti=t(1:n);
tf=t(n+1:n*2);
l=length(vin);
vout=vin;
vin=vin;
vout(1:n)=vin(1:n).*ti;
vout(l-(n-1):l)=vin(l-(n-1):l).*tf;

%%
function X=decodeQML(filexml)

% usage : X=decodeQML(filexml)
% read QuakeML files

% revision dec 17 2021
% Davide Piccinini

PRT=1;
%%
fid=fopen(filexml);
CC=textscan(fid,'%s','delimiter','\n');
fclose(fid);
KK=CC{1};

X=[];
n=0;

% DECODE LOCATION
for i=1:length(KK)
    LIN=char(KK(i));
    if findstr(LIN,'<origin publicID')==1
        OTime=char(KK(i+2));
        OTime=OTime(8:end-12);
        MOTime=datenum(OTime,'yyyy-mm-ddTHH:MM:SS.FFF');
    end
    if findstr(LIN,'<latitude')==1
        OLAT=char(KK(i+1)); LAT=OLAT(8:13);LAT=str2double(LAT);
    end
    if findstr(LIN,'<longitude')==1
        OLON=char(KK(i+1)); LON=OLON(8:13);LON=str2double(LON);
    end
    if findstr(LIN,'<depth>')==1        
        ODEP=char(KK(i+1)); fi=findstr(ODEP,'</value>');
        DEP=ODEP(8:fi-1);DEP=str2double(DEP)/1000;
    end
    if findstr(LIN,'<mag>')==1
        OMAG=char(KK(i+1));MAG=OMAG(8:10);MAG=str2double(MAG);
    end     
end
%%
% DECODE ARRIVALS
n=0;
for i=1:length(KK)
    LIN=char(KK(i));
    if findstr(LIN,'arrival publicID')==2
        n=n+1;
        PickID(n)=(i+1);
    end
end
fprintf('Number of arrivals: %d\n',n);
%%
n=0;
for i=1:length(KK)
    LIN=char(KK(i));
    if findstr(LIN,'arrival publicID')==2
        n=n+1;
        INI(n)=i;
    end
    if findstr(LIN,'</arrival>')==1
        FIN(n)=i;
    end
end
%%
n=0;
for i=1:length(FIN)
    BLK=KK(INI(i):FIN(i));
    for j=1:length(BLK)
        CC=char(BLK(j));
        if findstr(CC,'<timeWeight>')==1
            fin=findstr(CC,'</timeW');
            n=n+1;
            TW(n)=str2double(CC(13:fin-1));
        end
    end
end
%%
n=0;
for k=1:length(TW);
    if TW(k) >0
        BLK=KK(INI(k):FIN(k));
        for j=1:length(BLK)
            CC=char(BLK(j));
            if findstr(CC,'<pickID')==1
                f=findstr(CC,'=');
                ID=CC(f+1:end-9);
                n=n+1;
                PICKID(n)=str2double(ID);
            end
        end
    end
end
        
%%
for k=1:length(PICKID)
    PID=PICKID(k);
    PID=num2str(PID);
    for j=1:length(KK)
        LIN=char(KK(j));
        UU=findstr(LIN,PID);
        YY=findstr(LIN,'pick publicID');
        if YY==2 & isempty(UU)==0
            BLKI(k)=j;
            BLKF(k)=j+16;
            % pause
        end
    end
end
%%
nn=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(PICKID)   %%% RIMETTERE A POSTO QUESTA PARTE prende anche le fasi automatiche e NON DEVE FARLO  
    BLK=char(KK(BLKI(k):BLKF(k)));
    [r,c]=size(BLK);
    for jj=1:r
        LL=char(BLK(jj,:));
        i=strfind(LL,'<evaluationMode');
        if isempty(i)==0
            if strcmp(LL(17),'m')==1                
                MAN=1;
                nn=nn+1;
            else
                MAN=0;
            end
        end
    end
    if MAN==1    
        for jj=1:r
            LL=char(BLK(jj,:));
            i=strfind(LL,'stationCode');
            if isempty(i)==0
                staz=LL(i+13:i+13+3);
                if strcmp(staz(end),'"')==1
                    staz=staz(1:end-1);
                end
                if isfield(X,staz)==0
                    X.(staz)  =[];
                    X.(staz).P=[];
                    X.(staz).Punc=[];
                    X.(staz).S=[];
                    X.(staz).Sunc=[];
                end
                %pause
            end
        end    
        TIMEPICK=BLK(3,8:32);
        TIME    =datenum(TIMEPICK,'yyyy-mm-ddTHH:MM:SS.FFF');
        UNC     =str2double(BLK(4,14:16));
        PHASE   =BLK(9,12);
        if strcmp(PHASE,'P')
            X.(staz).P=TIME;
            X.(staz).Punc=UNC;
        end
        if strcmp(PHASE,'S')
            X.(staz).S=TIME;
            X.(staz).Sunc=UNC;
        end
    end
end
fprintf('Picks used : %d\n',nn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

if isempty(X)==0
    X.HEADER={OTime,MOTime,LAT,LON,DEP,MAG};
end

if PRT==1
    NM=fieldnames(X);
    NTOT=length(NM);
    for k=1:NTOT-1
        PT=X.(char(NM(k))).P;
        ST=X.(char(NM(k))).S;
        fprintf('%4s %s %s \n', char(NM(k)),datestr(PT,'yyyy-mm-ddTHH:MM:SS.FFF'),datestr(ST,'HH:MM:SS.FFF'))
    end
end


function OUT=ReadMSEED(filein)

% Read mseed files
% OUT=ReadMSEED(filein)
% OUT is a structure

        %keyboard
        %fprintf('Reading %s\n',filein)
        AA=rdmseed(filein);
        OUT.network=AA(1).NetworkCode;
        OUT.station=AA(1).StationIdentifierCode;
        OUT.channel=AA(1).ChannelIdentifier;
        OUT.dataquality='';
        OUT.location='';
        OUT.type='';
        OUT.sampleRate=AA(1).SampleRate;
        OUT.sampleType='i';
        DIFF=AA(1).RecordStartTimeMATLAB - round(AA(1).RecordStartTimeMATLAB);
        DIFF=abs(DIFF/86400);
        if abs(DIFF) < 1/OUT.sampleRate/86400; 
            AA(1).RecordStartTimeMATLAB=round(AA(1).RecordStartTimeMATLAB);            
        end
        OUT.dateTimeString=datestr(AA(1).RecordStartTimeMATLAB,'yyyy/mm/dd HH:MM:SS.FFF');    
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'yyyy'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'mm'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'dd'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'HH'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'MM'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'SS.FFF'));
        DL=0;
        for k=1:length(AA)
            DL=DL+length(AA(k).d);
            LEN(k)=AA(k).NumberSamples;
        end
        START=AA(1).t(1); 
        END  =AA(end).t(end);
        OUT.sampleCount=DL;
        OUT.numberOfSample=DL;
        OUT.matlabTimeVector=linspace(START,END,DL)';
        OUT.data=vertcat(AA.d);        

function []=findwav_mseed(fileday,filesamp,THR,PHASE);

%% usage: []=findwav_mseed(fileday,filesamp,THR,PHASE);
%% Based on Matching filter algorithm
%
%  filesamp = file sac che contiene la f.o. campione
%  fileday  = file sac nel quale bisogna effettuare la ricerca
%  fvar     = 1 se il filesamp ? una variabile o un file fisico [optional]
% 
% Controllare la durata del match (DUR)
% la larghezza del filtro 
% e il valore della CC
% 
% Crea in uscita un file "wavl.out"
%

% implementato riducendo considerevolmente
% la finestra di analisi

% Aggiunto il dimensionamento iniziale della variabile MC

% Mod. 16/07/08
% riordinata la parte della scrittura
% commentato alcune i blocchi principali
% aggiunto il sort delle var per i files di uscita

% Mod 31/10/2020 
% i files sono entrambi in memoria 


OUTFILESAC=0;

%keyboard
DAY=fileday;
%DAYNAME=datestr(DAY.matlabTimeVector(1),'yyyymmdd');
DAYNAME=datestr(datenum(DAY.dateTimeString,'yyyy/mm/dd HH:MM:SS.FFF'),'yyyymmdd');


if strcmpi(PHASE,'P')
    ONAME=datestr(filesamp.PTIME,'yyyymmddHHMMSS.FFF');
else
    ONAME=datestr(filesamp.STIME,'yyyymmddHHMMSS.FFF');
end
NAMEOUT=[filesamp.OTime '.' filesamp.station '.' filesamp.channel '.' PHASE '.' DAYNAME '.out'];

%if strcmp(filesamp.station,'SEI')==1
%keyboard
%end

if strcmpi(PHASE,'P')==1
    LW=length(filesamp.ptemplate);
    if (LW/2)-floor(LW/2) > 0   %% SE DISPARI
        WAVLET=filesamp.ptemplate;
        WAVLET=WAVLET(1:end-1);
        LW=length(WAVLET);
    else
        WAVLET=filesamp.ptemplate;
        LW=length(WAVLET);
    end
else
    LW=length(filesamp.stemplate);    
    if (LW/2)-floor(LW/2) > 0   %% SE DISPARI
        WAVLET=filesamp.stemplate;
        WAVLET=WAVLET(1:end-1);
        LW=length(WAVLET);
    else
        WAVLET=filesamp.stemplate;
        LW=length(WAVLET);
    end
end

if (LW/2)-floor(LW/2) > 0; fprintf('%d\n');end

LD=length(fileday.datafilt);
LCOR=(LW*2)-1;

SHIFT=floor(LW/2);
NWIN=floor((LD-LW)/SHIFT)+1;
%% inizializza MC !! <<-- questo fa guardagnare un sacco di tempo
%MC=zeros(LCOR,NWIN);
CCMAX=zeros(1,NWIN);
MAXTRAX=[];
%% INIZIALIZZA OUTSAC e STARTORI
OUTSAC=zeros(LW,NWIN);
STARTORI=zeros(1,NWIN);

L=length(DAY.datafilt);
NW=ceil(L/LW);

% disp(sprintf('*** CC THRESHOLD = %3.2f',THR)) 
fprintf('Matching started @ %s\n',DAY.station) 
t0 = clock;
N=0;
for J=1:NWIN;
    INI=(J-1)*SHIFT+1;
    FIN=INI+LW-1;
    TRAX=DAY.datafilt(INI:FIN);    
    [C,LAGS]=xcorr(WAVLET,TRAX,'coeff');
    [CMAX,JMAX]=max((C));  %% ABS o non ABS?
    if CMAX >= THR %% supera la soglia della CC --> FILTER MATCHED
%        keyboard
%        fprintf('%8.0f / %8.0f\n',J,NWIN);

        TMAX=round(LW-JMAX);
        I1=INI+TMAX; if I1 < 1; I1=1; end
        I2=I1 +LW-1;                
        TIM(J)=I1;      % PRIMO CAMPIONE DEL MATCH
        if TIM(J) > 0
            try
                MAXTRAX(J)=max(DAY.datafilt(I1:I2)); %% allineato sullo stack
            catch
                fprintf('%s %s -->> %d %d %d --- %d / %d \n',filesamp.station,PHASE,I1,I2,length(DAY.datafilt),J,NWIN)
                MAXTRAX(J)=max(DAY.datafilt(I1:length(DAY.datafilt)));
            end
            CCMAX  (J)=CMAX;
            N=N+1;
            STARTORI(N)=I1-1;
            
        end
    end
    % MC(:,J)=C; %% matrice delle CC
end

% disp(sprintf('...done in %d secs!',round(etime(clock,t0))));
fprintf('Corr @ %s for %s done in %5.0f secs!\n',DAY.station,PHASE,round(etime(clock,t0)));
%disp(' ');
%keyboard
%return
%% 
if isempty (MAXTRAX)==0
    %% ESTRAE SOLO I DATI  DHE HANNO MATCHATO IL FILTRO %%
    IDAMP=find(MAXTRAX > 0);% indice in base alle ampiezze
    CCOUT=CCMAX(find(CCMAX > 0));  %MAX della CC
    ITIM =TIM(find(TIM > 0));
    IDAMP=IDAMP';
    MT   =MAXTRAX(IDAMP)';
    %keyboard
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PER EVITARE CHE LO STESSO EVENTO SIA STATO MATCHATO 2 VOLTE
    %% CONSECUTIVAMENTE, CONTROLLO LE AMPIEZZE E TOLGO LE COPPIE DI VALORI
    %% UGUALI
    [B,I,J]=unique(MT,'legacy');
    %% ROUTINE PER SCEGLIERE I IN FUNZIONE DEL VALORE DI CC PIU' ALTO
    for k=1:max(J);
        i=find(J==k);
        if length(i)>1
            [v w]=max(CCOUT(i));
            igood=i(w);
            I(k)=igood;
%            k
%            pause
        end
    end
    %%
    NMT=MT(I);
    NIDAMP=IDAMP(I);
    NCCOUT=CCOUT(I);
    NTIM   =ITIM(I);
    NSTARTORI=STARTORI(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SORT DELLE USCITE IN FUNZIONE DEL TEMPO    
    [SNTIM,J] = sort(NTIM);
    SNMT      = NMT(J);
    SNIDAMP   = NIDAMP(J);
    SNCCOUT   = NCCOUT(J);
    SNSTARTORI= NSTARTORI(J);
    if SNSTARTORI(1)==0; SNSTARTORI(1)=1;end %% SE IL PRIMO CAMPIONE E' 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fout=fopen(['WORK/' NAMEOUT],'w');
    for i=1:length(SNCCOUT);
        OUTTIME=datestr(DAY.matlabTimeVector(SNSTARTORI(i)),'yyyy-mm-ddTHH:MM:SS.FFF');
        fprintf(fout,'%s %8.6e %5.3f\n',OUTTIME,SNMT(i), SNCCOUT(i));%
    end
    fclose (fout);
%    disp(sprintf('%2d %s matches found @ %s',length(SNMT),PHASE,DAY.station));
end
return


function varargout = rdmseed(varargin)
%RDMSEED Read miniSEED format file.
%	X = RDMSEED(F) reads file F and returns a M-by-1 structure X containing
%	M blocks ("data records") of a miniSEED file with headers, blockettes, 
%	and data in dedicated fields, in particular, for each data block X(i):
%		         t: time vector (DATENUM format)
%		         d: data vector (double)
%		BLOCKETTES: existing blockettes (substructures)
%
%	Known blockettes are 100, 500, 1000, 1001 and 2000. Others will be
%	ignored with a warning message.
%
%	X = RDMSEED(F,ENCODINGFORMAT,WORDORDER,RECORDLENGTH), when file F does 
%	not include the Blockette 1000 (like Seismic Handler outputs), specifies:
%		- ENCODINGFORMAT: FDSN code (see below); default is 10 = Steim-1;
%		- WORDORDER: 1 = big-endian (default), 0 = little-endian;
%		- RECORDLENGTH: must be a power of 2, at least 256 (default is 4096).
%	If the file contains Blockette 1000 (which is mandatory in the SEED 
%	convention...), these 3 arguments are ignored except with 'force' option.
%
%	X = RDMSEED without input argument opens user interface to select the 
%	file from disk.
%
%	[X,I] = RDMSEED(...) returns a N-by-1 structure I with N the detected 
%	number of different channels, and the following fields:
%	    ChannelFullName: channel name,
%	        XBlockIndex: channel's vector index into X,
%	         ClockDrift: vector of time interval errors, in seconds,
%	                     between each data block (relative to sampling
%	                     period). This can be compared to "Max Clock Drift"
%	                     value of a Blockette 52.
%	                        = 0 in perfect case
%	                        < 0 tends to overlapping
%	                        > 0 tends to gapping
%	  OverlapBlockIndex: index of blocks (into X) having a significant 
%	                     overlap with previous block (less than 0.5
%	                     sampling period).
%	        OverlapTime: time vector of overlapped blocks (DATENUM format).
%	      GapBlockIndex: index of blocks (into X) having a significant gap
%	                     with next block (more than 0.5 sampling period).
%	            GapTime: time vector of gapped blocks (DATENUM format).
%
%	RDMSEED(...) without output arguments plots the imported signal by 
%	concatenating all the data records, in one single plot if single channel
%	is detected, or subplots for multiplexed file (limited to 10 channels).
%	Gaps are shown with red stars, overlaps with green circles.
%
%	[...] = RDMSEED(F,...,'be') forces big-endian reading (overwrites the
%	automatic detection of endianness coding, which fails in some cases).
%
%	[...] = RDMSEED(F,...,'notc') disables time correction.
%
%	[...] = RDMSEED(F,...,'nullhead') ignores null header (some files may
%	start with a series of null bytes).
%
%	[...] = RDMSEED(F,...,'plot') forces the plot with output arguments.
%
%	[...] = RDMSEED(F,...,'v') uses verbose mode (displays additional 
%	information and warnings when necessary). Use 'vv' for extras, 'vvv'
%	for debuging.
%
%	Some instructions for usage of the returned structure:
%	
%	- to get concatenated time and data vectors from a single-channel file:
%		X = rdmseed(f,'plot');
%		t = cat(1,X.t);
%		d = cat(1,X.d);
%
%	- to get the list of channels in a multiplexed file:
%		[X,I] = rdmseed(f);
%		char(I.ChannelFullName)
%
%	- to extract the station component n from a multiplexed file:
%		[X,I] = rdmseed(f);
%		k = I(n).XBlockIndex;
%		plot(cat(1,X(k).t),cat(1,X(k).d))
%		datetick('x')
%		title(I(n).ChannelFullName)
%
%	Known encoding formats are the following FDSN codes:
%		 0: ASCII
%		 1: 16-bit integer
%		 2: 24-bit integer
%		 3: 32-bit integer
%		 4: IEEE float32
%		 5: IEEE float64
%		10: Steim-1
%		11: Steim-2
%		12: GEOSCOPE 24-bit (untested)
%		13: GEOSCOPE 16/3-bit gain ranged
%		14: GEOSCOPE 16/4-bit gain ranged
%		19: Steim-3 (alpha and untested)
%
%	See also MKMSEED to export data in miniSEED format.
%
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%		Institut de Physique du Globe de Paris
%	Created: 2010-09-17
%	Updated: 2018-08-09
%
%	Acknowledgments:
%		Ljupco Jordanovski, Jean-Marie Saurel, Mohamed Boubacar, Jonathan Berger,
%		Shahid Ullah, Wayne Crawford, Constanza Pardo, Sylvie Barbier,
%		Robert Chase, Arnaud Lemarchand, Alexandre Nercessian.
%
%		Special thanks to Martin Mityska who also inspired me with his ingenious
%		ReadMSEEDFast.m function.
%
%	References:
%		IRIS (2010), SEED Reference Manual: SEED Format Version 2.4, May 2010,
%		  IFDSN/IRIS/USGS, http://www.iris.edu
%		Trabant C. (2010), libmseed: the Mini-SEED library, IRIS DMC.
%		Steim J.M. (1994), 'Steim' Compression, Quanterra Inc.

%	History:
%
%		[2018-08-09]
%			- MAJOR CODE UPDATE: now processes the binary data in memory 
%			  after a global file reading.
%			- removes all global variables.
%		[2017-11-21]
%			- adds option 'nullhead' to bypass null bytes header.
%		[2015-01-05]
%			- fixes a bug when a data block has 0 sample declared in header
%			  but some data in the record (STEIM-1/2 coding).
%		[2014-06-29]
%			- 24-bit uncompressed format tested (bug correction), thanks to
%			  Arnaud Lemarchand.
%		[2014-05-31]
%			- applies the time correction to StartTime and X.t (if needed).
%			- new option 'notc' to disable time correction.
%			- Geoscope 16/4 format passed real data archive tests.
%			- fixes a problem when plotting multiplexed channels (thanks to 
%			  Robert Chase).
%		[2014-03-14]
%			- Improved endianness automatic detection (see comments).
%			- Accepts mixed little/big endian encoding in a single file.
%			- minor fixes.
%		[2013-10-25]
%			- Due to obsolete syntax of bitcmp(0,N) in R2013b, replaces all
%			  by: 2^N-1 (which is much faster...)
%		[2013-02-15]
%			- Tests also DayOfYear in header to determine automatically 
%			  little-endian coding of the file.
%			- Adds option 'be' to force big-endian reading (overwrites
%			  automatic detection).
%		[2012-12-21]
%			- Adds a verbose mode
%		[2012-04-21]
%			- Correct bug with Steim + little-endian coding
%			  (thanks to Shahid Ullah)
%		[2012-03-21]
%			- Adds IDs for warning messages
%		[2011-11-10]
%			- Correct bug with multiple channel name length (thanks to
%			  Jonathan Berger)
%		[2011-10-27]
%			- Add LocationIdentifier to X.ChannelFullName
%		[2011-10-24]
%			- Validation of IEEE double encoding (with PQL)
%			- Import/plot data even with file integrity problem (like PQL)
%		[2011-07-21]
%			- Validation of ASCII encoding format (logs)
%			- Blockettes are now stored in substructures below a single
%			  field X.BLOCKETTES
%			- Add import of blockettes 500 and 2000
%			- Accept multi-channel files with various data coding
%		[2010-10-16]
%			- Alpha-version of Steim-3 decoding...
%			- Extend output parameters with channel detection
%			- Add gaps and overlaps on plots
%			- Add possibility to force the plot
%		[2010-10-02]
%			- Add the input formats for GEOSCOPE multiplexed old data files
%			- Additional output argument with gap and overlap analysis
%			- Create a plot when no output argument are specified
%			- Optimize script coding (30 times faster STEIM decoding!)
%		[2010-09-28]
%			- Correction of a problem with STEIM-1 nibble 3 decoding (one 
%			  32-bit difference)
%			- Add reading of files without blockette 1000 with additional
%			  input arguments (like Seismic Handler output files).
%			- Uses warning() function instead of fprintf().
%
%	Copyright (c) 2018, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

if nargin > 6
	error('Too many input arguments.')
end

% default input arguments
makeplot = 0;	% make plot flag
verbose = 0;	% verbose flag/level 
forcebe = 0;	% force big-endian
ef = 10;		% encoding format default
wo = 1;			% word order default
rl = 2^12;		% record length default
force = 0;		% force input argument over blockette 1000 (UNDOCUMENTED)
notc = 0;		% force no time correction (over ActivityFlags)
nullhead = 0;	% allow null bytes before header

if nargin < 1
	[filename,pathname] = uigetfile('*','Please select a miniSEED file...');
	f = fullfile(pathname,filename);
else
	f = varargin{1};
end

if ~ischar(f) || ~exist(f,'file')
	error('File %s does not exist.',f);
end

if nargin > 1
	verbose = any(strcmpi(varargin,'v')) + 2*any(strcmpi(varargin,'vv')) ...
	          + 3*any(strcmpi(varargin,'vvv'));
	makeplot = any(strcmpi(varargin,'plot'));
	forcebe = any(strcmpi(varargin,'be'));
	notc = any(strcmpi(varargin,'notc'));
	force = any(strcmpi(varargin,'force'));
	nullhead = any(strcmpi(varargin,'nullhead'));
end
nargs = (makeplot>0) + (verbose>0) + (forcebe>0) + (notc>0) + (force>0) ...
	 + (nullhead>0);


if nargin > (1 + nargs)
	ef = varargin{2};
	if ~isnumeric(ef) || ~any(ef==[0:5,10:19,30:33])
		error('Argument ENCODINGFORMAT must be a valid FDSN code value.');
	end
end

if nargin > (2 + nargs)
	wo = varargin{3};
	if ~isnumeric(wo) || (wo ~= 0 && wo ~= 1)
		error('Argument WORDORDER must be 0 or 1.');
	end
end

if nargin > (3 + nargs)
	rl = varargin{4};
	if ~isnumeric(rl) || rl < 256 || rem(log(rl)/log(2),1) ~= 0
		error('Argument RECORDLENGTH must be a power of 2 and greater or equal to 256.');
	end
end

if nargout == 0
	makeplot = 1;
end

% sensible limits for multiplexed files
max_channels = 20;	% absolute max number of channels to plot
max_channel_label = 6;	% max. number of channels for y-labels

% file is opened in Big-Endian encoding (this is encouraged by SEED)
fid = fopen(f,'rb','ieee-be');
le = 0;
offset = 0;

% --- tests if the header is mini-SEED
% the 7th character must be one of the "data header/quality indicator", usually 'D'
header = fread(fid,20,'*char');
if ~ismember(header(7),'DRMQ')
	if ismember(header(7),'VAST')
		error('File seems to be a SEED Volume. Cannot read it.');
	else
		if header(1)==0
			if nullhead
				if verbose
					fprintf('Null header option: bypassing...');
				end
				c = 0;
				fseek(fid,0,'bof');
				while c==0
					c = fread(fid,1,'*char');
					offset = offset + 1;
				end
				if verbose
					fprintf(' %d null bytes.\n',offset);
				end
				header = fread(fid,6,'*char');
				if ~ismember(header(6),'DRMQ')
					error('File is not in mini-SEED format. Cannot read it.');
				else
					offset = offset - 1;
				end
			else
				error('File starts with null bytes... if you believe it is still a miniseed file, try the ''nullhead'' option.');
			end
		else
			error('File is not in mini-SEED format. Cannot read it.');
		end
	end
end

i = 1;

% --- main loop that reads data records until the end of the file
while offset >= 0
	[X(i),offset] = read_data_record(f,fid,offset,le,ef,wo,rl,forcebe,verbose,notc,force);
	i = i + 1;
end

fclose(fid);

if nargout > 0
	varargout{1} = X;
end

% --- analyses data
if makeplot || nargout > 1

	% test if the file is multiplexed or a single channel
	un = unique(cellstr(char(X.ChannelFullName)));
	nc = numel(un);
	for i = 1:nc
		k = find(strcmp(cellstr(char(X.ChannelFullName)),un{i}));
		I(i).ChannelFullName = X(k(1)).ChannelFullName;
		I(i).XBlockIndex = k;
		I(i).ClockDrift = ([diff(cat(1,X(k).RecordStartTimeMATLAB));NaN]*86400 - cat(1,X(k).NumberSamples)./cat(1,X(k).SampleRate))./cat(1,X(k).NumberSamples);
		I(i).OverlapBlockIndex = k(find(I(i).ClockDrift.*cat(1,X(k).NumberSamples).*cat(1,X(k).SampleRate) < -.5) + 1);
		I(i).OverlapTime = cat(1,X(I(i).OverlapBlockIndex).RecordStartTimeMATLAB);
		I(i).GapBlockIndex = k(find(I(i).ClockDrift.*cat(1,X(k).NumberSamples).*cat(1,X(k).SampleRate) > .5) + 1);
		I(i).GapTime = cat(1,X(I(i).GapBlockIndex).RecordStartTimeMATLAB);
	end
end
if nargout > 1
	varargout{2} = I;
end

% --- plots the data
if makeplot

	figure
	
	xlim = [min(cat(1,X.t)),max(cat(1,X.t))];

	% test if all data records have the same length
	rl = unique(cat(1,X.DataRecordSize));
	if numel(rl) == 1
		rl_text = sprintf('%d bytes',rl);
	else
		rl_text = sprintf('%d-%d bytes',min(rl),max(rl));
	end
	
	% test if all data records have the same sampling rate
	sr = unique(cat(1,X.SampleRate));
	if numel(sr) == 1
		sr_text = sprintf('%g Hz',sr);
	else
		sr_text = sprintf('%d # samp. rates',numel(sr));
	end
	
	% test if all data records have the same encoding format
	ef = unique(cellstr(cat(1,X.EncodingFormatName)));
	if numel(ef) == 1
		ef_text = sprintf('%s',ef{:});
	else
		ef_text = sprintf('%d different encod. formats',numel(ef));
	end
			
	if nc == 1
		plot(cat(1,X.t),cat(1,X.d))
		hold on
		for i = 1:length(I.GapBlockIndex)
			plot(I.GapTime(i),X(I.GapBlockIndex(i)).d(1),'*r')
		end
		for i = 1:length(I.OverlapBlockIndex)
			plot(I.OverlapTime(i),X(I.OverlapBlockIndex(i)).d(1),'og')
		end
		hold off
		set(gca,'XLim',xlim)
		datetick('x','keeplimits')
		grid on
		xlabel(sprintf('Time\n(%s to %s)',datestr(xlim(1)),datestr(xlim(2))))
		ylabel('Counts')
		title(sprintf('mini-SEED file "%s"\n%s (%d rec. @ %s - %g samp. @ %s - %s)', ...
			f,un{1},length(X),rl_text,numel(cat(1,X.d)),sr_text,ef_text),'Interpreter','none')
	else
		% plot is done only for real data channels...
		if nc > max_channels
			warning('Plot has been limited to %d channels (over %d). See help to manage multiplexed file.', ...
				max_channels,nc);
			nc = max_channels;
		end
		for i = 1:nc
			subplot(nc*2,1,i*2 + (-1:0))
			k = I(i).XBlockIndex;
			if ~any(strcmp('ASCII',cellstr(cat(1,X(k).EncodingFormatName))))
				plot(cat(1,X(k).t),cat(1,X(k).d))
				hold on
				for ii = 1:length(I(i).GapBlockIndex)
					if ~isempty(X(I(i).GapBlockIndex(ii)).d)
						plot(I(i).GapTime(ii),X(I(i).GapBlockIndex(ii)).d,'r')
					else
						plot(repmat(I(i).GapTime(ii),1,2),ylim,'r')
					end
				end
				for ii = 1:length(I(i).OverlapBlockIndex)
					if ~isempty(X(I(i).OverlapBlockIndex(ii)).d)
						plot(I(i).OverlapTime(ii),X(I(i).OverlapBlockIndex(ii)).d,'g')
					else
						plot(repmat(I(i).OverlapTime(ii),1,2),ylim,'g')
					end
				end
				hold off
			end
			set(gca,'XLim',xlim,'FontSize',8)
			h = ylabel(un{i},'Interpreter','none');
			if nc > max_channel_label
				set(gca,'YTick',[])
				set(h,'Rotation',0,'HorizontalAlignment','right','FontSize',8)
			end
			datetick('x','keeplimits')
			set(gca,'XTickLabel',[])
			grid on
			if i == 1
				title(sprintf('mini-SEED file "%s"\n%d channels (%d rec. @ %s - %g data - %s - %s)', ...
					f,length(un),length(X),rl_text,numel(cat(1,X(k).d)),sr_text,ef_text),'Interpreter','none')
			end
			if i == nc
				datetick('x','keeplimits')
				xlabel(sprintf('Time\n(%s to %s)',datestr(xlim(1)),datestr(xlim(2))))
			end
		end
		v = version;
		if str2double(v(1))>=7
			linkaxes(findobj(gcf,'type','axes'),'x')
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,offset] = read_data_record(f,fid,offset,le,ef,wo,rl,forcebe,verbose,notc,force)
% read_data_record uses global variables f, fid, offset, le, ef, wo, rl, 
%	and verbose. It reads a data record and returns a structure D.


fseek(fid,offset,'bof');

% --- read fixed section of Data Header (48 bytes)
D.SequenceNumber        = fread(fid,6,'*char')';
D.DataQualityIndicator  = fread(fid,1,'*char');
D.ReservedByte          = fread(fid,1,'*char');
D.StationIdentifierCode = fread(fid,5,'*char')';
D.LocationIdentifier    = fread(fid,2,'*char')';
D.ChannelIdentifier	    = fread(fid,3,'*char')';
D.NetworkCode           = fread(fid,2,'*char')';
D.ChannelFullName = sprintf('%s:%s:%s:%s',deblank(D.NetworkCode), ...
	deblank(D.StationIdentifierCode),deblank(D.LocationIdentifier), ...
	deblank(D.ChannelIdentifier));

% Start Time decoding
[D.RecordStartTime,swapflag] = readbtime(fid,forcebe);
D.RecordStartTimeISO = sprintf('%4d-%03d %02d:%02d:%07.4f',D.RecordStartTime);

if swapflag
	if le
		machinefmt = 'ieee-be';
		le = 0;
	else
		machinefmt = 'ieee-le';
		le = 1;
	end
	position = ftell(fid);
	fclose(fid);
	fid = fopen(f,'rb',machinefmt);
	fseek(fid,position,'bof');
	if verbose > 0
		warning('RDMSEED:DataIntegrity', ...
			'Sequence # %s: need to switch file encoding to %s...\n', ...
			D.SequenceNumber,machinefmt);
	end
end

D.NumberSamples	        = fread(fid,1,'uint16');

% Sample Rate decoding
SampleRateFactor        = fread(fid,1,'int16');
SampleRateMultiplier    = fread(fid,1,'int16');
if SampleRateFactor > 0
	if SampleRateMultiplier >= 0
		D.SampleRate = SampleRateFactor*SampleRateMultiplier;
	else
		D.SampleRate = -1*SampleRateFactor/SampleRateMultiplier;
	end
else
	if SampleRateMultiplier >= 0
		D.SampleRate = -1*SampleRateMultiplier/SampleRateFactor;
	else
		D.SampleRate = 1/(SampleRateFactor*SampleRateMultiplier);
	end
end

D.ActivityFlags          = fread(fid,1,'uint8');
D.IOFlags                = fread(fid,1,'uint8');
D.DataQualityFlags       = fread(fid,1,'uint8');
D.NumberBlockettesFollow = fread(fid,1,'uint8');
D.TimeCorrection         = fread(fid,1,'int32');	% Time correction in 0.0001 s
D.OffsetBeginData        = fread(fid,1,'uint16');
D.OffsetFirstBlockette   = fread(fid,1,'uint16');

% --- read the blockettes
OffsetNextBlockette = D.OffsetFirstBlockette;

D.BLOCKETTES = [];
b2000 = 0;	% Number of Blockette 2000

for i = 1:D.NumberBlockettesFollow
	fseek(fid,offset + OffsetNextBlockette,'bof');
	BlocketteType = fread(fid,1,'uint16');
	
	switch BlocketteType
		
		case 1000
			% BLOCKETTE 1000 = Data Only SEED (8 bytes)
			OffsetNextBlockette = fread(fid,1,'uint16');
			D.BLOCKETTES.B1000.EncodingFormat = fread(fid,1,'uint8');
			D.BLOCKETTES.B1000.WordOrder = fread(fid,1,'uint8');
			D.BLOCKETTES.B1000.DataRecordLength = fread(fid,1,'uint8');
			D.BLOCKETTES.B1000.Reserved = fread(fid,1,'uint8');
			
		case 1001
			% BLOCKETTE 1001 = Data Extension (8 bytes)
			OffsetNextBlockette = fread(fid,1,'uint16');
			D.BLOCKETTES.B1001.TimingQuality = fread(fid,1,'uint8');
			D.BLOCKETTES.B1001.Micro_sec = fread(fid,1,'int8');
			D.BLOCKETTES.B1001.Reserved = fread(fid,1,'uint8');
			D.BLOCKETTES.B1001.FrameCount = fread(fid,1,'uint8');
			
		case 100
			% BLOCKETTE 100 = Sample Rate (12 bytes)
			OffsetNextBlockette = fread(fid,1,'uint16');
			D.BLOCKETTES.B100.ActualSampleRate = fread(fid,1,'float32');
			D.BLOCKETTES.B100.Flags = fread(fid,1,'uint8');
			D.BLOCKETTES.B100.Reserved = fread(fid,1,'uint8');
		
		case 500
			% BLOCKETTE 500 = Timing (200 bytes)
			OffsetNextBlockette = fread(fid,1,'uint16');
			D.BLOCKETTES.B500.VCOCorrection = fread(fid,1,'float32');
			D.BLOCKETTES.B500.TimeOfException = readbtime(fid,forcebe);
			D.BLOCKETTES.B500.MicroSec = fread(fid,1,'int8');
			D.BLOCKETTES.B500.ReceptionQuality = fread(fid,1,'uint8');
			D.BLOCKETTES.B500.ExceptionCount = fread(fid,1,'uint16');
			D.BLOCKETTES.B500.ExceptionType = fread(fid,16,'*char')';
			D.BLOCKETTES.B500.ClockModel = fread(fid,32,'*char')';
			D.BLOCKETTES.B500.ClockStatus = fread(fid,128,'*char')';
		
		case 2000
			% BLOCKETTE 2000 = Opaque Data (variable length)
			b2000 = b2000 + 1;
			OffsetNextBlockette = fread(fid,1,'uint16');
			BlocketteLength = fread(fid,1,'uint16');
			OffsetOpaqueData = fread(fid,1,'uint16');
			D.BLOCKETTES.B2000(b2000).RecordNumber = fread(fid,1,'uint32');
			D.BLOCKETTES.B2000(b2000).DataWordOrder = fread(fid,1,'uint8');
			D.BLOCKETTES.B2000(b2000).Flags = fread(fid,1,'uint8');
			NumberHeaderFields = fread(fid,1,'uint8');
			HeaderFields = splitfield(fread(fid,OffsetOpaqueData-15,'*char')','~');
			D.BLOCKETTES.B2000(b2000).HeaderFields = HeaderFields(1:NumberHeaderFields);
			% Opaque data are stored as a single char string, but must be
			% decoded using appropriate format (e.g., Quanterra Q330)
			D.BLOCKETTES.B2000(b2000).OpaqueData = fread(fid,BlocketteLength-OffsetOpaqueData,'*char')';
		
		otherwise
			OffsetNextBlockette = fread(fid,1,'uint16');

			if verbose > 0
				warning('RDMSEED:UnknownBlockette', ...
					'Unknown Blockette number %d (%s)!\n', ...
					BlocketteType,D.ChannelFullName);
			end
	end
end

% --- read the data stream
fseek(fid,offset + D.OffsetBeginData,'bof');

if ~force && isfield(D.BLOCKETTES,'B1000')
	EncodingFormat = D.BLOCKETTES.B1000.EncodingFormat;
	WordOrder = D.BLOCKETTES.B1000.WordOrder;
	D.DataRecordSize = 2^D.BLOCKETTES.B1000.DataRecordLength;
else
	EncodingFormat = ef;
	WordOrder = wo;
	D.DataRecordSize = rl;
end

uncoded = 0;

D.d = NaN;
D.t = NaN;

switch EncodingFormat
	
	case 0
		% --- decoding format: ASCII text
		D.EncodingFormatName = {'ASCII'};
		D.d = fread(fid,D.DataRecordSize - D.OffsetBeginData,'*char')';

	case 1
		% --- decoding format: 16-bit integers
		D.EncodingFormatName = {'INT16'};
		dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/2),'*int16');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		D.d = dd(1:D.NumberSamples);
		
	case 2
		% --- decoding format: 24-bit integers
		D.EncodingFormatName = {'INT24'};
		dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/3),'bit24=>int32');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		D.d = dd(1:D.NumberSamples);
		
	case 3
		% --- decoding format: 32-bit integers
		D.EncodingFormatName = {'INT32'};
		dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/4),'*int32');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		D.d = dd(1:D.NumberSamples);
		
	case 4
		% --- decoding format: IEEE floating point
		D.EncodingFormatName = {'FLOAT32'};
		dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/4),'*float');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		D.d = dd(1:D.NumberSamples);
		
	case 5
		% --- decoding format: IEEE double precision floating point
		D.EncodingFormatName = {'FLOAT64'};
		dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/8),'*double');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		D.d = dd(1:D.NumberSamples);

	case {10,11,19}
		% --- decoding formats: STEIM-1 and STEIM-2 compression
		%    (c) Joseph M. Steim, Quanterra Inc., 1994
		steim = find(EncodingFormat==[10,11,19]);
		D.EncodingFormatName = {sprintf('STEIM%d',steim)};
		
		% Steim compression decoding strategy optimized for Matlab
		% -- by F. Beauducel, October 2010 --
		%
		%	1. loads all data into a single 16xM uint32 array
		%	2. gets all nibbles from the first row splitted into 2-bit values
		%	3. for each possible nibble value, selects (find) and decodes
		%	   (bitsplit) all the corresponding words, and stores results
		%	   in a 4xN (STEIM1) or 7xN (STEIM2) array previously filled with
		%	   NaN's. For STEIM2 with nibbles 2 or 3, decodes also dnib values
		%	   (first 2-bit of the word)
		%	5. reduces this array with non-NaN values only
		%	6. integrates with cumsum
		%
		% This method is about 30 times faster than a 'C-like' loops coding...
		
		frame32 = fread(fid,[16,(D.DataRecordSize - D.OffsetBeginData)/64],'*uint32');
		if xor(~WordOrder,le)
			frame32 = swapbytes(frame32);
		end
		
		% specific processes for STEIM-3
		if steim == 3
			% first bit = 1 means second differences
			SecondDiff = bitshift(frame32(1,:),-31);
			% checks for "squeezed flag"... and replaces frame32(1,:)
			squeezed = bitand(bitshift(frame32(1,:),-24),127);
			k = find(bitget(squeezed,7));
			if ~isempty(k)
				moredata24 = bitand(frame32(1,k),16777215);
				k = find(squeezed == 80);	% upper nibble 8-bit = 0x50
				if ~isempty(k)
					frame32(1,k) = hex2dec('15555555');
				end
				k = find(squeezed == 96);	% upper nibble 8-bit = 0x60
				if ~isempty(k)
					frame32(1,k) = hex2dec('2aaaaaaa');
				end
				k = find(squeezed == 112);	% upper nibble 8-bit = 0x70
				if ~isempty(k)
					frame32(1,k) = hex2dec('3fffffff');
				end
			end
		end

		% nibbles is an array of the same size as frame32...
		nibbles = bitand(bitshift(repmat(frame32(1,:),16,1),repmat(-30:2:0,size(frame32,2),1)'),3);
		x0 = bitsign(frame32(2,1),32);	% forward integration constant
		xn = bitsign(frame32(3,1),32);	% reverse integration constant
		
		switch steim
			
		case 1
			% STEIM-1: 3 cases following the nibbles
			ddd = NaN*ones(4,numel(frame32));	% initiates array with NaN
			k = find(nibbles == 1);			% nibble = 1 : four 8-bit differences
			if ~isempty(k)
				ddd(1:4,k) = bitsplit(frame32(k),32,8);
			end
			k = find(nibbles == 2);			% nibble = 2 : two 16-bit differences
			if ~isempty(k)
				ddd(1:2,k) = bitsplit(frame32(k),32,16);
			end
			k = find(nibbles == 3);			% nibble = 3 : one 32-bit difference
			if ~isempty(k)
				ddd(1,k) = bitsign(frame32(k),32);
			end

		case 2	
			% STEIM-2: 7 cases following the nibbles and dnib
			ddd = NaN*ones(7,numel(frame32));	% initiates array with NaN
			k = find(nibbles == 1);			% nibble = 1 : four 8-bit differences
			if ~isempty(k)
				ddd(1:4,k) = bitsplit(frame32(k),32,8);
			end
			k = find(nibbles == 2);			% nibble = 2 : must look in dnib
			if ~isempty(k)
				dnib = bitshift(frame32(k),-30);
				kk = k(dnib == 1);		% dnib = 1 : one 30-bit difference
				if ~isempty(kk)
					ddd(1,kk) = bitsign(frame32(kk),30);
				end
				kk = k(dnib == 2);		% dnib = 2 : two 15-bit differences
				if ~isempty(kk)
					ddd(1:2,kk) = bitsplit(frame32(kk),30,15);
				end
				kk = k(dnib == 3);		% dnib = 3 : three 10-bit differences
				if ~isempty(kk)
					ddd(1:3,kk) = bitsplit(frame32(kk),30,10);
				end
			end
			k = find(nibbles == 3);				% nibble = 3 : must look in dnib
			if ~isempty(k)
				dnib = bitshift(frame32(k),-30);
				kk = k(dnib == 0);		% dnib = 0 : five 6-bit difference
				if ~isempty(kk)
					ddd(1:5,kk) = bitsplit(frame32(kk),30,6);
				end
				kk = k(dnib == 1);		% dnib = 1 : six 5-bit differences
				if ~isempty(kk)
					ddd(1:6,kk) = bitsplit(frame32(kk),30,5);
				end
				kk = k(dnib == 2);		% dnib = 2 : seven 4-bit differences (28 bits!)
				if ~isempty(kk)
					ddd(1:7,kk) = bitsplit(frame32(kk),28,4);
				end
			end
			
		case 3	% *** STEIM-3 DECODING IS ALPHA AND UNTESTED ***
			% STEIM-3: 7 cases following the nibbles
			ddd = NaN*ones(9,numel(frame32));	% initiates array with NaN
			k = find(nibbles == 0);				% nibble = 0 : two 16-bit differences
			if ~isempty(k)
				ddd(1:2,k) = bitsplit(frame32(k),32,16);
			end
			k = find(nibbles == 1);				% nibble = 1 : four 8-bit differences
			if ~isempty(k)
				ddd(1:4,k) = bitsplit(frame32(k),32,8);
			end
			k = find(nibbles == 2);				% nibble = 2 : must look even dnib
			if ~isempty(k)
				dnib2 = bitshift(frame32(k(2:2:end)),-30);
				w60 = bitand(frame32(k(2:2:end)),1073741823) ...
					+ bitshift(bitand(frame32(k(1:2:end)),1073741823),30);	% concatenates two 30-bit words
				kk = find(dnib2 == 0);		% dnib = 0: five 12-bit differences (60 bits)
				if ~isempty(kk)
					ddd(1:5,k(2*kk)) = bitsplit(w60,60,12);
				end
				kk = find(dnib2 == 1);		% dnib = 1: three 20-bit differences (60 bits)
				if ~isempty(kk)
					ddd(1:3,k(2*kk)) = bitsplit(w60,60,20);
				end
			end
			k = find(nibbles == 3);				% nibble = 3 : must look 3rd bit
			if ~isempty(k)
				dnib = bitshift(frame32(k),-27);
				kk = k(dnib == 24);		% dnib = 11000 : nine 3-bit differences (27 bits)
				if ~isempty(kk)
					ddd(1:9,kk) = bitsplit(frame32(kk),27,3);
				end
				kk = k(dnib == 25);		% dnib = 11001 : Not A Difference
				if ~isempty(kk)
					ddd(1,kk) = bitsign(frame32(kk),27);
				end
				kk = k(dnib > 27);		% dnib = 111.. : 29-bit sample (29 bits)
				if ~isempty(kk)
					ddd(1,kk) = bitsign(frame32(kk),29);
				end
			end
		end
		
		% Little-endian coding: needs to swap bytes
		if ~WordOrder
			ddd = flipud(ddd);
		end
		dd = ddd(~isnan(ddd));		% reduces initial array ddd: dd is non-NaN values of ddd
		
		% controls the number of samples
		if numel(dd) ~= D.NumberSamples
			if verbose > 1
				warning('RDMSEED:DataIntegrity','Problem in %s sequence # %s [%d-%03d %02d:%02d:%07.4f]: number of samples in header (%d) does not equal data (%d).\n', ...
					D.EncodingFormatName{:},D.SequenceNumber,D.RecordStartTimeISO,D.NumberSamples,numel(dd));
			end
			if numel(dd) < D.NumberSamples
				D.NumberSamples = numel(dd);
			end
		end

		% rebuilds the data vector by integrating the differences
		D.d = cumsum([x0;dd(2:D.NumberSamples)]);
		
		% controls data integrity...
		if D.d(end) ~= xn
			warning('RDMSEED:DataIntegrity','Problem in %s sequence # %s [%s]: data integrity check failed, last_data=%d, Xn=%d.\n', ...
				D.EncodingFormatName{:},D.SequenceNumber,D.RecordStartTimeISO,D.d(end),xn);
		end

		if D.NumberSamples == 0
			D.d = nan(0,1);
		end
		
		% for debug purpose...
		if verbose > 2
			D.dd = dd;
			D.nibbles = nibbles;
			D.x0 = x0;
			D.xn = xn;
		end

	case 12
		% --- decoding format: GEOSCOPE multiplexed 24-bit integer
		D.EncodingFormatName = {'GEOSCOPE24'};
		dd = fread(fid,(D.DataRecordSize - D.OffsetBeginData)/3,'bit24=>double');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		D.d = dd(1:D.NumberSamples);
		
	case {13,14}
		% --- decoding format: GEOSCOPE multiplexed 16/3 and 16/4 bit gain ranged
		%	(13): 16/3-bit (bit 15 is unused)
		%	(14): 16/4-bit
		%	bits 15-12 = 3 or 4-bit gain exponent (positive) 
		%	bits 11-0 = 12-bit mantissa (positive)
		%	=> data = (mantissa - 2048) / 2^gain
		geoscope = 7 + 8*(EncodingFormat==14); % mask for gain exponent
		D.EncodingFormatName = {sprintf('GEOSCOPE16-%d',EncodingFormat-10)};
		dd = fread(fid,(D.DataRecordSize - D.OffsetBeginData)/2,'*uint16');
		if xor(~WordOrder,le)
			dd = swapbytes(dd);
		end
		dd = (double(bitand(dd,2^12-1))-2^11)./2.^double(bitand(bitshift(dd,-12),geoscope));
		D.d = dd(1:D.NumberSamples);
		
	case 15
		% --- decoding format: US National Network compression
		D.EncodingFormatName = {'USNN'};
		uncoded = 1;
		
	case 16
		% --- decoding format: CDSN 16-bit gain ranged
		D.EncodingFormatName = {'CDSN'};
		uncoded = 1;
		
	case 17
		% --- decoding format: Graefenberg 16-bit gain ranged
		D.EncodingFormatName = {'GRAEFENBERG'};
		uncoded = 1;
		
	case 18
		% --- decoding format: IPG - Strasbourg 16-bit gain ranged
		D.EncodingFormatName = {'IPGS'};
		uncoded = 1;
		
	case 30
		% --- decoding format: SRO format
		D.EncodingFormatName = {'SRO'};
		uncoded = 1;
		
	case 31
		% --- decoding format: HGLP format
		D.EncodingFormatName = {'HGLP'};
		uncoded = 1;
		
	case 32
		% --- decoding format: DWWSSN gain ranged format
		D.EncodingFormatName = {'DWWSSN'};
		uncoded = 1;
		
	case 33
		% --- decoding format: RSTN 16-bit gain ranged
		D.EncodingFormatName = {'RSTN'};
		uncoded = 1;
		
	otherwise
		D.EncodingFormatName = {sprintf('** Unknown (%d) **',EncodingFormat)};
		uncoded = 1;
		
end

if uncoded
	error('Sorry, the encoding format "%s" is not yet implemented.',D.EncodingFormatName);
end

% Applies time correction (if needed)
D.RecordStartTimeMATLAB = datenum(double([D.RecordStartTime(1),0,D.RecordStartTime(2:5)])) ...
	+ (~notc & bitand(D.ActivityFlags,2) == 0)*D.TimeCorrection/1e4/86400;
tv = datevec(D.RecordStartTimeMATLAB);
doy = datenum(tv(1:3)) - datenum(tv(1),1,0);
D.RecordStartTime = [tv(1),doy,tv(4:5),round(tv(6)*1e4)/1e4];
D.RecordStartTimeISO = sprintf('%4d-%03d %02d:%02d:%07.4f',D.RecordStartTime);

D.t = D.RecordStartTimeMATLAB;

% makes the time vector and applies time correction (if needed)
if EncodingFormat > 0
	D.t = D.t + (0:(D.NumberSamples-1))'/(D.SampleRate*86400);
end


offset = ftell(fid);
fread(fid,1,'char');	% this is to force EOF=1 on last record.
if feof(fid)
	offset = -1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = splitfield(s,d)
% splitfield(S) splits string S of D-character separated field names

C = textscan(s,'%s','Delimiter',d);
c = C{1};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,swapflag] = readbtime(fid,forcebe)
% readbtime reads BTIME structure from current opened file and returns
%	D = [YEAR,DAY,HOUR,MINUTE,SECONDS]

Year		= fread(fid,1,'*uint16');
DayOfYear	= fread(fid,1,'*uint16');
Hours		= fread(fid,1,'uint8');
Minutes		= fread(fid,1,'uint8');
Seconds		= fread(fid,1,'uint8');
fseek(fid,1,0);	% skip 1 byte (unused)
Seconds0001	= fread(fid,1,'*uint16');

% Automatic detection of little/big-endian encoding
% -- by F. Beauducel, March 2014 --
%
% If the 2-byte day is >= 512, the file is not opened in the correct
% endianness. If the day is 1 or 256, there is a possible byte-swap and we
% need to check also the year; but we need to consider what is a valid year:
% - years from 1801 to 2047 are OK (swapbytes >= 2312)
% - years from 2048 to 2055 are OK (swapbytes <= 1800)
% - year 2056 is ambiguous (swapbytes = 2056)
% - years from 2057 to 2311 are OK (swapbytes >= 2312)
% - year 1799 is ambiguous (swapbytes = 1799)
% - year 1800 is suspicious (swapbytes = 2055)
%
% Thus, the only cases for which we are 'sure' there is a byte-swap, are:
% - day >= 512
% - (day == 1 or day == 256) and (year < 1799 or year > 2311)
%
% Note: in IRIS libmseed, the test is only year>2050 or year<1920.
if ~forcebe && (DayOfYear >= 512 || (ismember(DayOfYear,[1,256]) && (Year > 2311 || Year < 1799)))
	swapflag = 1;
	Year = swapbytes(Year);
	DayOfYear = swapbytes(DayOfYear);
	Seconds0001 = swapbytes(Seconds0001);
else
	swapflag = 0;
end
d = [double(Year),double(DayOfYear),Hours,Minutes,Seconds + double(Seconds0001)/1e4];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = bitsplit(x,b,n)
% bitsplit(X,B,N) splits the B-bit number X into signed N-bit array
%	X must be unsigned integer class
%	N ranges from 1 to B
%	B is a multiple of N

sign = repmat((b:-n:n)',1,size(x,1));
x = repmat(x',b/n,1);
d = double(bitand(bitshift(x,flipud(sign-b)),2^n-1)) ...
	- double(bitget(x,sign))*2^n;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = bitsign(x,n)
% bitsign(X,N) returns signed double value from unsigned N-bit number X.
% This is equivalent to bitsplit(X,N,N), but the formula is simplified so
% it is much more efficient

d = double(bitand(x,2^n-1)) - double(bitget(x,n)).*2^n;



%% PYTHON CODE 

% 
% #!/usr/bin/env python3
% # -*- coding: utf-8 -*-
% """
% Created on Thu Nov 12 15:25:09 2020
% 
% @author: piccinini
% """
% import os
% from obspy.clients.fdsn import Client
% from obspy import UTCDateTime
% client = Client("http://verna.pi.ingv.it:8080")
% 
% CSC_INI="2016-03-31T00:00:00"  # Data Inizio Cascinelle
% 
% 
% latmin= 42.859
% lonmin= 11.710
% latmax= 42.891
% lonmax= 11.756
% 
% # ESTRAGGO IL CATALOGO DA VERNA
% client = Client('http://verna.pi.ingv.it:8080')
% ST = UTCDateTime(CSC_INI)
% ED = UTCDateTime(CSC_INI)+(86400)
% try:
%   cat = client.get_events(starttime=ST, endtime=ED,includearrivals=True)
%   print(str(len(cat)) + ' eventi scaricati')
% except:
%   print('No events in catalog')
% 
% 
% n=0
% for ev in cat:
% #    if ev.event_type=='earthquake':
% #    print(ev)
%     if ev.origins[0].evaluation_mode=='manual' : 
%         
%         n=n+1
%         TOUT=ev.origins[0].time
%         OUT=str(TOUT)[0:-4] + '.xml'
%         ev.write('TEMPLATES/'+OUT,format='QUAKEML')
%         WVF=client.get_waveforms('TV','*','*','?HZ',starttime=TOUT-30,endtime=TOUT+30)
%         WVF.merge(fill_value='interpolate')
%         for fil in WVF:
%           fil.write('TEMPLATES/' +str(TOUT)[0:-4]+ '.' +fil.stats.station+ '.ms',format='MSEED')
%         WVF=client.get_waveforms('IV','ARCI','*','?HZ',starttime=TOUT-30,endtime=TOUT+30)
%         WVF.merge(fill_value='interpolate')
%         for fil in WVF:
%           fil.write('TEMPLATES/' +str(TOUT)[0:-4]+ '.' +fil.stats.station+ '.ms',format='MSEED')
%         WVF=client.get_waveforms('IV','MCIV','*','?HZ',starttime=TOUT-30,endtime=TOUT+30)
%         WVF.merge(fill_value='interpolate')
%         for fil in WVF:
%           fil.write('TEMPLATES/' +str(TOUT)[0:-4]+ '.' +fil.stats.station+ '.ms',format='MSEED')
% 
% print (str(n) + ' eventi scritti')
% 
% 
% #CSC_INI="2016-04-02T00:00:00"  # Data Inizio Cascinelle
% #client = Client('http://verna.pi.ingv.it:8080')
% #ST = UTCDateTime(CSC_INI)
% #ED = UTCDateTime(CSC_INI)+(86400)
% 
% 
% OUTDIR=CSC_INI[0:4]+CSC_INI[5:7]+CSC_INI[8:10]
% try:
%   os.mkdir(OUTDIR)
% except:
%   print('Dir esiste')
% WVF=client.get_waveforms('TV','*','*','?HZ',starttime=ST,endtime=ED)
% WVF.merge(fill_value='interpolate')
% for sta in WVF:
%     #sta.merge(fill_value='interpolate')
%     sta.write(OUTDIR+'/'+CSC_INI+'.'+sta.stats.station+'.'+sta.stats.network+'.'+sta.stats.channel,format='MSEED')
% WVF=client.get_waveforms('IV','ARCI','*','?HZ',starttime=ST,endtime=ED)
% WVF.merge(fill_value='interpolate')
% for sta in WVF:
%     #sta.merge(fill_value='interpolate')
%     sta.write(OUTDIR+'/'+CSC_INI+'.'+sta.stats.station+'.'+sta.stats.network+'.'+sta.stats.channel,format='MSEED')
% WVF=client.get_waveforms('IV','MCIV','*','?HZ',starttime=ST,endtime=ED)
% WVF.merge(fill_value='interpolate')
% for sta in WVF:
%     #sta.merge(fill_value='interpolate')
%     sta.write(OUTDIR+'/'+CSC_INI+'.'+sta.stats.station+'.'+sta.stats.network+'.'+sta.stats.channel,format='MSEED')
%             


