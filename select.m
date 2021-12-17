clear all

% SELECT AND WRITE OUT IN NLLOC CORRELATION OUTPUT
D=dir('*.mat');

% D=dir('2009-01-23T23:15:46.940.20090123.mat')    
%%
CCTRH=[.65 ]

for JJ=1:length(CCTRH)
    L=0;
    THR=CCTRH(JJ);
    %THR=0.65
    NP = 4;
    NS = 3;
    
    for k=1:length(D)
        S=load(D(k).name);
        S.SLAVE=S.NEWSTRUCT;
        N=length(S.SLAVE);
        if N > 0;
            for j=1:N
                NUMP=length(S.SLAVE(j).Pcc);
                NUMS=length(S.SLAVE(j).Scc);
                MP= mean(S.SLAVE(j).Pcc);
                MS= mean(S.SLAVE(j).Scc);
                MTP=(S.SLAVE(j).Ptim(1));
                MTS=(S.SLAVE(j).Stim(1));
                
                if MP > THR && MS > THR && NUMP >= NP && NUMS >= NS
                    L=L+1;
                    SLA(L).Psta= S.SLAVE(j).Psta;
                    SLA(L).Ptim= S.SLAVE(j).Ptim;
                    SLA(L).Pcc=  S.SLAVE(j).Pcc;
                    SLA(L).Ssta= S.SLAVE(j).Ssta;
                    SLA(L).Stim= S.SLAVE(j).Stim;
                    SLA(L).Scc=  S.SLAVE(j).Scc;
                    SLA(L).TemplateOTime=S.SLAVE(j).TemplateOTime;
                    SLA(L).MTP=MTP;
                    SLA(L).MTS=MTS;
                    SLA(L).MP=MP;
                    SLA(L).MS=MS;
                end
            end
        end
    end
    
    fprintf('Total # of events with CC > %3.2f & NP >= %d & NS >= %d : %d\n',THR,NP,NS,length(SLA))
    
    T=struct2table(SLA);
    S=sortrows(T,8);
    TIMEP=S{:,8};
    % DTIMEP=diff(TIMEP)/86400;
    TOL=2/86400;
    [C,ID,IJ] = uniquetol(TIMEP,TOL/max(abs(TIMEP)));
    
    NMAX=max(IJ);
    TOERASE=[];
    
    for j=1:NMAX
        NSTAP=[];NSTAS=[];CCP=[];CCS=[];
        II=find(IJ==j);
        for k=1:length(II)
            NSTAP(k)=length(T{II(k),'Psta'}{:});
            NSTAS(k)=length(T{II(k),'Ssta'}{:});
            CCP  (k)=(T{II(k),'MP'});
            CCS  (k)=(T{II(k),'MS'});
        end
        Z=[NSTAP ; NSTAS; CCP; CCS];
        SCORE=sum(Z);
        [w v]=max(SCORE);
        CHOICE=II(v);
        
        OUT=setdiff(II,CHOICE);
        LO=length(OUT);
        LE=length(TOERASE);
        TOERASE(LE+1:LE+length(OUT))=OUT ;
        %pause
    end
    
    [r c]=size(S);
    if exist('TOERASE','var') == 1
        S([TOERASE],:)=[]; % TOGLIE GLI INDICI CHE HA SCARTATO
        fprintf('%4.0f double events removed\n',length(TOERASE))
    else
        fprintf('No doubles found, %4d events keeped\n',r);
    end
    
    SLAVE=table2struct(S);
    [r,c]=size(S);
    fprintf('Total slaves collected: %d\n',r)
    RR(JJ)=r
    
    clear SLA TOERASE II
    
end

%%
for k=1:length(SLAVE); OTIME(k)=SLAVE(k).TemplateOTime;TIME(k)=SLAVE(k).Ptim(1);CCP(k)=SLAVE(k).MP; end

U=unique(OTIME);
clear TU
clear CC

for k=1:length(U);
    i=find(OTIME==U(k));
    TU(k,2:length(i)+1)=NaN;
    TU(k,1)=TIME(k);
    CC(k,2:length(i)+1)=NaN;
    CC(k,1)=CCP(k);
    
    for j=1:length(i)
        TU(k,j+1)=TIME(i(j));
        CC(k,j+1)=CCP(i(j));
    end
    i=[];
end

[r c]=size(TU)

figure;
subplot(211)
for k=1:r
    ii=find(TU(k,:) > 0);
    X=k.*ones(size(TU(k,ii)));
    scatter(TU(k,ii(2:end)),X(2:end),30,CC(k,ii(2:end)),'s','filled'); hold on
end
datetick('x','mmdd','keeplimits')
box on
grid on


subplot(212); plot(TIME,1:length(TIME),'sk','MarkerFaceColor','r','MarkerSize',4);datetick('x','mmdd','keeplimits');grid on
xlabel('Time')
ylabel('Number of events')

%%
N=17
fid=fopen('lista.template','w'); 
for k=1:length(TU(N,:));
    if TU(N,k) > 0
    fprintf(fid,'%s\n',datestr(TU(N,k),'yyyy-mm-ddTHH:MM:SS.FFF'));
    end
end; 
fclose(fid);


%% % WRITE OUT NNLINLOC PHS FILE

mkdir('NLL_OBS')
STHR=sprintf('%3.2f',THR);

for k=1:length(SLAVE)
    SLAVE(k).Psta;
    NAME=datestr(SLAVE(k).Ptim(1),'yyyymmdd_HHMMSS.FFF');
    fout=fopen(['NLL_OBS/' NAME '.' STHR '.xml'],'w');
    for j=1:length(SLAVE(k).Psta);
        PSTA=SLAVE(k).Psta(j);
        PCH ='Z';
        PHS ='P';
        TIME=datestr(SLAVE(k).Ptim(j),'yyyymmdd HHMM SS.FFF');
        WW  =(1-SLAVE(k).Pcc(j))/3;
        fprintf(fout,'%-6s ?    %s    ? %s      ? %s0 GAU %9.2e -1.00e+00 -1.00e+00 -1.00e+00\n',...
            char(PSTA),PCH,PHS,TIME,WW);
    end
    for j=1:length(SLAVE(k).Ssta);
        SSTA=SLAVE(k).Ssta(j);
        SCH ='Z';
        PHS ='S';
        TIME=datestr(SLAVE(k).Stim(j),'yyyymmdd HHMM SS.FFF');
        WW  =(1-SLAVE(k).Scc(j))/3;
        fprintf(fout,'%-6s ?    %s    ? %s      ? %s0 GAU %9.2e -1.00e+00 -1.00e+00 -1.00e+00\n',...
            char(SSTA),SCH,PHS,TIME,WW);
    end
    fclose(fout);
end

%ZIPFILE=sprintf('zip FIRENZUOLA_%1d_%1d_%4.2f.zip *.xml',NP,NS,THR);
%eval(ZIPFILE)









