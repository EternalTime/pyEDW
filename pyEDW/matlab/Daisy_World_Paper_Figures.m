%   Parameter setup
%   theta (Parameter Vector)

f       =   0.88;
Ag      =   0.3;
Ab      =   0.1;
Aw      =   0.6;
Topt    =   300;
Q       =   0.1;
gG      =   1;
gD      =   0.2;
tauS    =   3;
tauE    =   5;
dT      =   30;

dt = 0.1;

theta0 = [ ...
            f , ...             %    1:     f
            1-Ag , ...          %    2:     1 - A_G
            Ab-Ag, ...          %    3:     A_B - A_G
            Aw-Ag, ...          %    4:     A_W - A_G
            Q, ...              %    5:     Q
            gD/gG, ...          %    6:     gamma_D / gamma_G
            1/(gG*tauE), ...    %    7:     1 / gamma_G tau_E
            1/(gG*tauS), ...    %    8:     1 / gamma_G tau_S
            8*(Topt/dT)^4, ...  %    9:     8 (T_opt / Delta T)^4
            0.05, ...           %    10:    delta
            0 ...               %    11:    lambda
        ];

transPink = [245 169 184]/255;

dTs = linspace(0,80,129);
dTs(1) = [];
delTs = dTs/Topt;
limT = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899]/theta0(9)^0.25;
cx = repmat(fliplr(linspace(0,1,50))',[1 3]);cy = [flipud(cx);cx];
cx = [cx;zeros(50,3)];cz = flipud(cx);

myCmap = repmat([1 1 1],[100,1]).*cx+repmat(SkyBlue,[100,1]).*cy +...
    repmat(transPink,[100 1]).*cz;


%% Limits of variables


fig = figure("Position",[0 0 800 600]);
ax1 = axes(fig,"Position",[0.1 0.1 0.38 0.38]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;
ax2 = axes(fig,"Position",[0.55 0.1 0.38 0.38]);
hold(ax2,"on");box(ax2,"on");ax2.LineWidth = 1.5;
ax3 = axes(fig,"Position",[0.1 0.55 0.38 0.38]);
hold(ax3,"on");box(ax3,"on");ax3.LineWidth = 1.5;
ax4 = axes(fig,"Position",[0.55 0.55 0.38 0.38]);
hold(ax4,"on");box(ax4,"on");ax4.LineWidth = 1.5;
Nstart = 7;
delT = dTs/Topt;
for Nd = Nstart:128
    load("data_"+num2str(Nd,"%03.0f")+".mat")
    f1 = squeeze(data(:,:,1));f1 = f1(:);
    f2 = squeeze(data(:,:,2));f2 = f2(:);
    L = squeeze(data(:,:,4));L = L(:);
    T = squeeze(data(:,:,3));T = T(:);
    idx = (f1>0.001)&(f2>0.001);
    sum(idx)

    f1 = f1(idx);
    f2 = f2(idx);
    L = L(idx);
    T = T(idx);

    minF = 0.8;maxF = 1.1;
    f1min = minF*min(f1);f1max = maxF*max(f1);
    %df1 =2*std(f1);f1min = max(0,mean(f1)-df1);f1max = min(f,mean(f1)+df1);
    f2min = minF*min(f2);f2max = maxF*max(f2);
    %df2 =1.5*std(f2);f2min = max(0,mean(f2)-df2);f2max = min(f,mean(f2)+df2);
    %Tmin = minF*min(T);Tmax = maxF*max(T);
    dT1 =3.5*std(T);Tmin = mean(T)-dT1;Tmax = mean(T)+dT1;
    Lmin = minF*min(L);Lmax = maxF*max(L)-0.2;
    %dL1 =3*std(L);Lmin = mean(L)-dL1;Lmax = mean(L)+dL1;
    if Nd ==Nstart
        plt1 = plot(ax1,delT(Nd)*[1 1],[f1min, f1max],'.',"Color",RochesterBlue,"MarkerSize",10);
        plt2 = plot(ax2,delT(Nd)*[1 1],[f2min, f2max],'.',"Color",ForagerPink,"MarkerSize",10);
        plt3 = plot(ax3,delT(Nd)*[1 1],[Tmin, Tmax],'.',"Color",ResourceGreen,"MarkerSize",10);
        plt4 = plot(ax4,delT(Nd)*[1 1],[Lmin, Lmax],'.',"Color",DartmouthGreen,"MarkerSize",10);
    else
        plt1.XData = [plt1.XData delT(Nd)*[1 1]];
        plt2.XData = [plt2.XData delT(Nd)*[1 1]];
        plt3.XData = [plt3.XData delT(Nd)*[1 1]];
        plt4.XData = [plt4.XData delT(Nd)*[1 1]];
        plt1.YData = [plt1.YData f1min f1max];
        plt2.YData = [plt2.YData f2min f2max];
        plt3.YData = [plt3.YData Tmin Tmax];
        plt4.YData = [plt4.YData Lmin Lmax];
        drawnow
        pause(eps);
    end
end

%% DATA COLLECTION: Extracting Efficacy and Viability
% conditioned on average luminosoty

load("data_000.mat");
eLmu = sum(data(:,:,4),2)/500;
eT = data(:,:,3);
eTmu = sum(data(:,:,3),2)/500;
Efficacy = zeros([400,128]);
pB = zeros([400,128]);
pB1 = zeros([400,128]);
pB2 = zeros([400,128]);
f1 = zeros([400,128]);
f2 = zeros([400,128]);
avgT = zeros([2,400,128]);
for ii = 1:128
    tic
    fname = "data_"+num2str(ii,"%03.0f")+".mat";
    load(fname,"data");
    idx = (data(:,:,1)+data(:,:,2))>10^-3;
    Tmu = sum(data(:,:,3).*idx,2)./max(3,sum(idx,2));
    eTmu = sum(eT.*idx,2)./max(3,sum(idx,2));
    avgT(:,:,ii) = [Tmu, eTmu]';
    Efficacy(:,ii) = Tmu./eTmu-1;
    pB(:,ii) = mean(idx,2);
    pB1(:,ii) = mean(data(:,:,1)>10^-3,2);
    pB2(:,ii) = mean(data(:,:,2)>10^-3,2);
    f1(:,ii) = mean(data(:,:,1).*(data(:,:,1)>10^-3),2);
    f2(:,ii) = mean(data(:,:,2).*(data(:,:,2)>10^-3),2);
    toc
end
eT = repmat(mean(eT,2),[1 128]);
Efficacy(isnan(Efficacy)) = 0;
aT = squeeze(avgT(1,:,:));
aT(aT == 0) = eT(aT == 0);
%
kernel = [1 7 1]';kernel = kernel/sum(kernel(:));
for ii = 1:2
    % only vertical smoothing
    Efficacy = conv2(padarray(Efficacy,[1 0],"replicate"),kernel,"valid");
    pB = conv2(padarray(pB,[1 0],"replicate"),kernel,"valid");
    aT = conv2(padarray(aT,[1 0],"replicate"),kernel,"valid");
end
kernel = [1 2 1;2 24 2;1 2 1];kernel = kernel/sum(kernel(:));
for ii = 1:0
    Efficacy = conv2(padarray(Efficacy,[1 1],"replicate"),kernel,"valid");
    pB = conv2(padarray(pB,[1 1],"replicate"),kernel,"valid");
end
Efficacy(pB<10^-3)=nan;
f1(pB1<10^-3) = 0;
f2(pB2<10^-3) = 0;
idx = (f1==0)&(f2==0);

aT(idx) = eT(idx);



save("data_temperaturesV.mat","avgT","eT","aT","Efficacy","Ls","Ts","delTs","pB","pB1","pB2","f1","f2")

%%
%% DATA COLLECTION: Extracting Biotic Characteristics

load("data_000.mat");
eLmu = sum(data(:,:,4),2)/500;
eT = data(:,:,3);
eT = repmat(mean(eT,2),[1 128]);
eTmu = sum(data(:,:,3),2)/500;
pB = zeros([400,128]);
pB1 = zeros([400,128]);
pB2 = zeros([400,128]);
f1 = zeros([400,128]);
f2 = zeros([400,128]);
f1f2 = zeros([400,128]);
f1pf2 = zeros([400,128]);
S = zeros([400,128]);
aT = eT;
Rmat = zeros(400,128,4,4);
for ii = 1:128
    tic
    fname = "data_"+num2str(ii,"%03.0f")+".mat";
    load(fname,"data");
    idx = (data(:,:,1)+data(:,:,2))>2*10^-3;
    aT(:,ii) = sum(data(:,:,3).*idx,2)./max(3,sum(idx,2));
    pB(:,ii) = mean(idx,2);
    pB1(:,ii) = mean(data(:,:,1)>10^-3,2);
    pB2(:,ii) = mean(data(:,:,2)>10^-3,2);
    f1(:,ii) = mean(data(:,:,1).*(data(:,:,1)>10^-3),2);
    f2(:,ii) = mean(data(:,:,2).*(data(:,:,2)>10^-3),2);
    f1s(:,ii) = mean(data(:,:,1).^2.*(data(:,:,1)>10^-3),2);
    f2s(:,ii) = mean(data(:,:,2).^2.*(data(:,:,2)>10^-3),2);
    f1f2(:,ii) = mean(data(:,:,1).*data(:,:,2).*(sum(data,3)>10^-3),2);
    f1pf2(:,ii) = mean((data(:,:,1)+data(:,:,2)).*(sum(data,3)>10^-3),2);
    S(:,ii) = -mean((data(:,:,1).*log2(data(:,:,1))+data(:,:,2).*log2(data(:,:,2))).*((data(:,:,1)>10^-3)&(data(:,:,2)>10^-3)),2);
    for jj = 1:400
        Rmat(jj,ii,:,:) = cov(squeeze(data(jj,:,:)));
    end
    toc
end
idx = (aT<min(eT(:))) | ((eT>1) & (aT<0.86));
aT(idx) = eT(idx);
%
kernel = [1 2 1];kernel = kernel/sum(kernel(:));
for ii = 1:0
    % only vertical smoothing
    pB = conv2(padarray(pB,[0 1],"replicate"),kernel,"valid");
    pB1 = conv2(padarray(pB1,[0 1],"replicate"),kernel,"valid");
    pB2 = conv2(padarray(pB2,[0 1],"replicate"),kernel,"valid");
    f1 = conv2(padarray(f1,[0 1],"replicate"),kernel,"valid");
    f2 = conv2(padarray(f2,[0 1],"replicate"),kernel,"valid");
    f1s = conv2(padarray(f1s,[0 1],"replicate"),kernel,"valid");
    f2s = conv2(padarray(f2s,[0 1],"replicate"),kernel,"valid");
    f1f2 = conv2(padarray(f1f2,[0 1],"replicate"),kernel,"valid");
    f1pf2 = conv2(padarray(f1pf2,[0 1],"replicate"),kernel,"valid");
    aT = conv2(padarray(aT,[0 1],"replicate"),kernel,"valid");
    S = conv2(padarray(S,[0 1],"replicate"),kernel,"valid");
end
kernel = [1 2 1;2 24 2;1 2 1];kernel = kernel/sum(kernel(:));
for ii = 1:0
    pB = conv2(padarray(pB,[1 1],"replicate"),kernel,"valid");
    pB1 = conv2(padarray(pB1,[1 1],"replicate"),kernel,"valid");
    pB2 = conv2(padarray(pB2,[1 1],"replicate"),kernel,"valid");
    f1 = conv2(padarray(f1,[1 1],"replicate"),kernel,"valid");
    f2 = conv2(padarray(f2,[1 1],"replicate"),kernel,"valid");
    f1s = conv2(padarray(f1s,[1 1],"replicate"),kernel,"valid");
    f2s = conv2(padarray(f2s,[1 1],"replicate"),kernel,"valid");
    f1f2 = conv2(padarray(f1f2,[1 1],"replicate"),kernel,"valid");
    f1pf2 = conv2(padarray(f1pf2,[1 1],"replicate"),kernel,"valid");
    aT = conv2(padarray(aT,[1 1],"replicate"),kernel,"valid");
    S = conv2(padarray(S,[1 1],"replicate"),kernel,"valid");
end


save("data_temperaturesV.mat","avgT","eT","aT","Ls","Ts","delTs","pB","pB1","pB2","f1","f2","f1s","f2s","f1f2","f1pf2","S","Rmat")

%% difference in probability of survival between phenotypes
close all 
load("data_efficacy.mat");
fig = figure("Position",[0  0 600,600]);
ax1 = axes(fig,"Position",[0.11 0.15 0.75 0.76]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;

im1 = surf(ax1,delTs,Ls(1:(end)),10*f1.*f2,"EdgeColor","none","FaceAlpha",1);
colormap(ax1,myCmap);
cb = colorbar(ax1,"Position",[sum(ax1.Position([1 3]))+0.02 ax1.Position(2) 0.025 ax1.Position(4)],...
        "LineWidth", 1.5,"TickLabelInterpreter","latex","FontSize",12);
cb.Label.String = "$|p(f_W>0)-p(f_B>0)|$";
cb.Label.FontSize = 18;
cb.Label.Interpreter = 'latex';
set(ax1,"CLim",[0 1]);
cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95]));
cts(end) = [];

%[~,ct] = contour3(ax1,delTs,Ls,(1+pB)*max(Efficacy(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);     
%[~,ct1] = contour3(ax1,delTs,Ls,(1+pB1),cts,"LineWidth",0.8,"Color",RochesterBlue);     
%[~,ct2] = contour3(ax1,delTs,Ls,(1+pB2),cts,"LineWidth",0.8,"Color",ForagerPink);   
contour3(ax1,delTs,Ls,1+f1/f,1+logspace(-5,0,20),"Color",[1 0 0],"LineWidth",2);
contour3(ax1,delTs,Ls,1+f2/f,1+logspace(-5,0,20),"Color",[0 0 1],"LineWidth",2);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);
plot3(ax1,[0 5],[1 1],[1 1],"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax1,delTs,p3u,ones(size(delTs)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax1,delTs,p3l,ones(size(delTs)),"--","Color",[0 0 0],"LineWidth",1.5)
set(ax1,"YDir","normal","XLim",[0,max(delTs)+0.001],...
    "YLim",[min(Ls)-0.001,Ls(end-1)+0.001],"LineWidth",1.5,"FontSize",16)
xlabel(ax1,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax1,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")


%%
close all
fig = figure("Position",[0 0 900 900]);
ax = axes(fig,"Position",[0.1 0.1 0.85 0.85]);
delT = dTs/Topt;
plot(ax,[0 3],[0 3],"-","Color",[0 0 0],"LineWidth",2.0)
hold(ax,"on");box(ax,"on");
pltB = plot(ax,[1+delT(1) 1+delT(1) 1-delT(1) 1-delT(1)],...
    [1-delT(1) 1+delT(1) 1+delT(1) 1-delT(1)],"--","Color",[0 0 0],"LineWidth",2.0);
set(ax,"XLim",[0.5 1.4],"YLim",[0.5 1.4],"LineWidth",1.5)

for kk = [1 40 80 128]
load("data_"+num2str(kk,"%03.0f")+".mat")
T = data(:,:,3);
L = data(:,:,4);
Te = L.^0.25;
if kk == 1
    plt = plot(ax,Te(1,:),T(1,:),'.',"Color",[0 0 0],"MarkerSize",10);
else
    plt.XData = Te(1,:);
    plt.YData = T(1,:);
    pltB.XData = [1+delT(kk) 1+delT(kk) 1-delT(kk) 1-delT(kk)];
    pltB.YData = [1-delT(kk) 1+delT(kk) 1+delT(kk) 1-delT(kk)];
end
xlabel(ax,num2str(kk,"%03.0f"),"FontSize",25)
drawnow
pause(eps)
for nn = 2:400
    plt.XData = Te(nn,:);
    plt.YData = T(nn,:);
    drawnow
    pause(eps)
end
hold("off")
end


%% DATA COLLECTION : extract all the entropies of AE systems


H = zeros(128,400,16);
getH = @(p) -sum(p(p(:)>0).*(log2(p(p(:)>0))));
for kk = 1:128
    tic
    load("data_"+num2str(kk,"%03.0f")+".mat");
    for ll = 1:400
        idx = squeeze(data(ll,:,1).*data(ll,:,2))>0.0000001;
        if sum(idx)>1
            Nbins = 1+ceil(sqrt(sum(idx)));
            f1 = squeeze(data(ll,idx,1));
            f2 = squeeze(data(ll,idx,2));
            T = squeeze(data(ll,idx,3));
            L = squeeze(data(ll,idx,4));

            x1 = linspace(min(f1),max(f1),Nbins-1);
            x2 = linspace(min(f2),max(f2),Nbins-1);
            x3 = linspace(min(T),max(T),Nbins-1);
            x4 = linspace(min(L),max(L),Nbins-1);

            pAE = zeros(Nbins,Nbins,Nbins,Nbins)/Nbins^4;
            for ii = 1:length(f1)
                idx1 = 1+sum(f1(ii)>x1);
                idx2 = 1+sum(f2(ii)>x2);
                idx3 = 1+sum(T(ii)>x3);
                idx4 = 1+sum(L(ii)>x4);
                pAE(idx1,idx2,idx3,idx4) =  pAE(idx1,idx2,idx3,idx4) +1;
            end
            pAE = pAE/sum(pAE(:));
            
            H(kk,ll,1) = getH(squeeze(sum(pAE,[2 3 4])));
            H(kk,ll,2) = getH(squeeze(sum(pAE,[1 3 4])));
            H(kk,ll,3) = getH(squeeze(sum(pAE,[1 2 4])));
            H(kk,ll,4) = getH(squeeze(sum(pAE,[1 2 3])));

            H(kk,ll,5) = getH(squeeze(sum(pAE,[1 2])));
            H(kk,ll,6) = getH(squeeze(sum(pAE,[1 3])));
            H(kk,ll,7) = getH(squeeze(sum(pAE,[1 4])));
            H(kk,ll,8) = getH(squeeze(sum(pAE,[2 3])));
            H(kk,ll,9) = getH(squeeze(sum(pAE,[2 4])));
            H(kk,ll,10) = getH(squeeze(sum(pAE,[3 4])));

            H(kk,ll,11) = getH(squeeze(sum(pAE,1)));
            H(kk,ll,12) = getH(squeeze(sum(pAE,2)));
            H(kk,ll,13) = getH(squeeze(sum(pAE,3)));
            H(kk,ll,14) = getH(squeeze(sum(pAE,4)));

            H(kk,ll,15) = getH(pAE);
            H(kk,ll,16) = Nbins;
            
        end
    end
    [kk toc]
end
save("data_entropy3.mat","H","Ls","Ts","Nbins","delTs")
  
%% DATA COLLECTION: Extract entropies of E
load("data_000.mat")
Nbins = 1+ceil(sqrt(500));

H = zeros(400,4);
getH = @(p) -sum(p(p(:)>0).*(log2(p(p(:)>0))));
for ll = 1:400
    T = squeeze(data(ll,:,3));
    L = squeeze(data(ll,:,4));

    x3 = linspace(min(T),max(T),Nbins-1);
    x4 = linspace(min(L),max(L),Nbins-1);

    pE = zeros(Nbins,Nbins);
    for ii = 1:length(f1)
        idx3 = 1+sum(T(ii)>x3);
        idx4 = 1+sum(L(ii)>x4);
        pE(idx3,idx4) =  pE(idx3,idx4) +1;
    end
    pE = pE/sum(pE(:));

    H(ll,1) = getH(squeeze(sum(pE,2)));
    H(ll,2) = getH(squeeze(sum(pE,1)));

    H(ll,3) = getH(pE);
    H(ll,4) = Nbins;
end
save("data_entropyE.mat","H","Ls","Ts","Nbins","delTs")
%% E:E correlations drop with appearance of Agent
 
load("data_efficacy.mat");
load("data_entropyE.mat");
IE = H(:,1) + H(:,2) - H(:,3);
for kk =1:2
IE = conv(padarray(IE,[3 0],"replicate"),[1 1 3 5 3 1 1]/15,"valid"); 
end
load("data_entropy3.mat","H");

idx = (pB1>10^-3)&(pB2>10^-3);
%IAE = H(:,:,10) + H(:,:,5) - H(:,:,15);
IEE = H(:,:,3) + H(:,:,4) - H(:,:,5);

dIE = 1-IEE'./repmat(IE,[1 128]);
dIE(not(idx)) = 0;
kernel = [1 2 1;2 4 2;1 2 1];kernel = kernel/sum(kernel(:));
for kk =1:5
    dIE = conv2(padarray(dIE,[1 1],"replicate"),kernel,"valid");
end
dIE(not(idx)) = nan;

close all
fig = figure("Position",[0 0 650 600]);
ax = axes(fig,"Position",[0.13 0.13*6.5/6 0.73 0.73*6.5/6]);
surf(ax,delTs,Ls,dIE,"EdgeColor","none","FaceAlpha",1.0)
colormap(ax,myCmap);view(ax,[0 90]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
cb = colorbar(ax,"Position",[sum(ax.Position([1 3]))+0.02 ax.Position(2) 0.025 ax.Position(4)]);
set(ax,"CLim",[0 max(dIE(:))*1.05])
cb.Label.String = "$-\Delta I_E/I_E$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 16; 
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.5;

%surf(ax,delTs,Ls,100*Efficacy,"FaceColor",[0.8 0 0],"EdgeColor","none","FaceAlpha",0.5)
%surf(ax,delTs,Ls,repmat(IE,[1 128]),"FaceColor",[0.8 0 0],"EdgeColor","none","FaceAlpha",0.5)

cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95])*max(dIE(:)));
cts(end) = [];
[~,ct] = contour3(ax,delTs,Ls,(1+f1+f2)*max(dIE(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);    
%contour(ax1,delTs,Ts,Efficacy,25,"Color",[0 0 0]);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);

plot3(ax,[0 5],[1 1],[1 1]*max(dIE(:)),"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax,delTs,p3u,ones(size(delTs))*max(dIE(:)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax,delTs,p3l,ones(size(delTs))*max(dIE(:)),"--","Color",[0 0 0],"LineWidth",1.5)


%surf(ax,delTs,Ls,IEE',"FaceColor",[0 0 0.8],"EdgeColor","none","FaceAlpha",0.5)
set(ax,"XLim",[0,max(delTs)-0.01],"YLim",[min(Ls)+0.1 max(Ls)-0.2],"FontSize",16)
xlabel(ax,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")
grid(ax,"off")
%% A:A correlations 
 
load("data_efficacy.mat");

idx = pB>10^-2;
IAA = (H(:,:,1) + H(:,:,2) - H(:,:,10))';
kernel = [1 2 1;2 4 2;1 2 1];kernel = kernel/sum(kernel(:));
for kk = 1:3
    IAA = conv2(padarray(IAA,[1 1],"replicate"),kernel,"valid");
end
IAA(not(idx)) = nan;

close all
fig = figure("Position",[0 0 650 600]);
ax = axes(fig,"Position",[0.13 0.13*6.5/6 0.73 0.73*6.5/6]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
surf(ax,delTs,Ls,IAA,"EdgeColor","none","FaceAlpha",1.0)
colormap(ax,myCmap)

cb = colorbar(ax,"Position",[sum(ax.Position([1 3]))+0.02 ax.Position(2) 0.025 ax.Position(4)]);
set(ax,"CLim",[0 max(IAA(:))*1.05])
cb.Label.String = "$I(A:A)$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 16; 
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.5;

cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95])*max(IAA(:)));
cts(end) = [];
[~,ct] = contour3(ax,delTs,Ls,(1+pB)*max(IAA(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);

plot3(ax,[0 5],[1 1],[1 1]*max(IAA(:)),"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax,delTs,p3u,ones(size(delTs))*max(IAA(:)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax,delTs,p3l,ones(size(delTs))*max(IAA(:)),"--","Color",[0 0 0],"LineWidth",1.5)

%surf(ax,delTs,Ls,IEE',"FaceColor",[0 0 0.8],"EdgeColor","none","FaceAlpha",0.5)
set(ax,"XLim",[0,max(delTs)-0.01],"YLim",[min(Ls)+0.1 max(Ls)-0.2],"FontSize",16)
xlabel(ax,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%% A:E correlations 
 
load("data_efficacy.mat");

idx = pB>10^-3;
IAE = (H(:,:,5) + H(:,:,10) - H(:,:,15))';
kernel = [1 2 1;2 4 2;1 2 1];kernel = kernel/sum(kernel(:));
IAE = conv2(padarray(IAE,[1 1],"replicate"),kernel,"valid");
IAE(not(idx)) = nan;

close all
fig = figure("Position",[0 0 650 600]);
ax = axes(fig,"Position",[0.13 0.13*6.5/6 0.73 0.73*6.5/6]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
surf(ax,delTs,Ls,IAE,"EdgeColor","none","FaceAlpha",1.0)
colormap(ax,myCmap)

cb = colorbar(ax,"Position",[sum(ax.Position([1 3]))+0.02 ax.Position(2) 0.025 ax.Position(4)]);
set(ax,"CLim",[0 max(IAE(:))*1.05])
cb.Label.String = "$I(E:A)$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 16; 
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.5;

cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95])*max(IAE(:)));
cts(end) = [];
[~,ct] = contour3(ax,delTs,Ls,(1+f1+f2)*max(IAE(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);

plot3(ax,[0 5],[1 1],[1 1]*max(IAE(:)),"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax,delTs,p3u,ones(size(delTs))*max(IAE(:)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax,delTs,p3l,ones(size(delTs))*max(IAE(:)),"--","Color",[0 0 0],"LineWidth",1.5)

%surf(ax,delTs,Ls,IEE',"FaceColor",[0 0 0.8],"EdgeColor","none","FaceAlpha",0.5)
set(ax,"XLim",[0,max(delTs)-0.01],"YLim",[min(Ls)+0.1 max(Ls)-0.2],"FontSize",16)
xlabel(ax,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%% f1:E correlations 
 
load("data_efficacy.mat");

idx = pB1>10^-3;
IA1E = (H(:,:,1) + H(:,:,5) - H(:,:,12))';
%kernel = [1 2 1;2 4 2;1 2 1];kernel = kernel/sum(kernel(:));
%IA1E = conv2(padarray(IAE,[1 1],"replicate"),kernel,"valid");
IA1E(not(idx)) = nan;

close all
fig = figure("Position",[0 0 650 600]);
ax = axes(fig,"Position",[0.13 0.13*6.5/6 0.73 0.73*6.5/6]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
surf(ax,delTs,Ls,IA1E,"EdgeColor","none","FaceAlpha",1.0)
colormap(ax,myCmap)

cb = colorbar(ax,"Position",[sum(ax.Position([1 3]))+0.02 ax.Position(2) 0.025 ax.Position(4)]);
set(ax,"CLim",[min(IA1E(:))*1.05 max(IA1E(:))*1.05])
cb.Label.String = "$I(E:f_B)$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 16; 
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.5;

cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95])*max(IA1E(:)));
cts(end) = [];
[~,ct] = contour3(ax,delTs,Ls,(1+pB1)*max(IA1E(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);

plot3(ax,[0 5],[1 1],[1 1]*max(IA1E(:)),"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax,delTs,p3u,ones(size(delTs))*max(IA1E(:)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax,delTs,p3l,ones(size(delTs))*max(IA1E(:)),"--","Color",[0 0 0],"LineWidth",1.5)

%surf(ax,delTs,Ls,IEE',"FaceColor",[0 0 0.8],"EdgeColor","none","FaceAlpha",0.5)
set(ax,"XLim",[0,max(delTs)-0.01],"YLim",[min(Ls)+0.1 max(Ls)-0.2],"FontSize",16)
xlabel(ax,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%% f2:E correlations 
 
load("data_efficacy.mat");

idx = pB2>10^-3;
IA2E = (H(:,:,2) + H(:,:,5) - H(:,:,11))';
%kernel = [1 2 1;2 4 2;1 2 1];kernel = kernel/sum(kernel(:));
%IAE = conv2(padarray(IAE,[1 1],"replicate"),kernel,"valid");
IA2E(not(idx)) = nan;

close all
fig = figure("Position",[0 0 650 600]);
ax = axes(fig,"Position",[0.13 0.13*6.5/6 0.73 0.73*6.5/6]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
surf(ax,delTs,Ls,IA2E,"EdgeColor","none","FaceAlpha",1.0)
colormap(ax,myCmap)

cb = colorbar(ax,"Position",[sum(ax.Position([1 3]))+0.02 ax.Position(2) 0.025 ax.Position(4)]);
set(ax,"CLim",[min(IA2E(:))*1.05 max(IA2E(:))*1.05])
cb.Label.String = "$I(E:f_W)$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 16; 
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.5;

cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95])*max(IA2E(:)));
cts(end) = [];
[~,ct] = contour3(ax,delTs,Ls,(1+pB2)*max(IA2E(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);

plot3(ax,[0 5],[1 1],[1 1]*max(IA2E(:)),"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax,delTs,p3u,ones(size(delTs))*max(IA2E(:)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax,delTs,p3l,ones(size(delTs))*max(IA2E(:)),"--","Color",[0 0 0],"LineWidth",1.5)

%surf(ax,delTs,Ls,IEE',"FaceColor",[0 0 0.8],"EdgeColor","none","FaceAlpha",0.5)
set(ax,"XLim",[0,max(delTs)-0.01],"YLim",[min(Ls)+0.1 max(Ls)-0.2],"FontSize",16)
xlabel(ax,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%% Syn(f2 f2:E) Redundancy between f1 and f2 in affecting E
 
load("data_efficacy.mat");

idx = (pB2>10^-8)&(pB1>10^-8);

Synergy = (H(:,:,15) + H(:,:,1) + H(:,:,2) + H(:,:,5) - H(:,:,10)-H(:,:,11)-H(:,:,12))';
kernel = [1 2 1;2 4 2;1 2 1];kernel = kernel/sum(kernel(:));
Synergy = conv2(padarray(Synergy,[1 1],"replicate"),kernel,"valid");
Synergy(not(idx)) = nan;

close all
fig = figure("Position",[0 0 650 600]);
ax = axes(fig,"Position",[0.13 0.13*6.5/6 0.73 0.73*6.5/6]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
surf(ax,delTs,Ls,Synergy,"EdgeColor","none","FaceAlpha",1.0)

r=abs(min(Synergy(:))/max(Synergy(:)));

cm = flipud([flipud(myCmap);myCmap(1:ceil(r*100),:)]);
colormap(ax,cm)

cb = colorbar(ax,"Position",[sum(ax.Position([1 3]))+0.02 ax.Position(2) 0.025 ax.Position(4)]);
set(ax,"CLim",[1.05*min(Synergy(:))*1.05 max(Synergy(:))*1.05])
cb.Label.String = "$I(A\rightarrow E)$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 16; 
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.5;

cts = sort((1+[logspace(-1.8,log10(0.5),12) 0.75 0.95])*max(Synergy(:)));
cts(end) = [];
[~,ct] = contour3(ax,delTs,Ls,(1+f1+f2)*max(Synergy(:)),cts,"LineWidth",0.8,"Color",[0 0 0]);
[~,ct2] = contour3(ax,delTs,Ls,max(Synergy(:))+Synergy,max(Synergy(:))+[-eps,0,eps],"LineWidth",1.5,"Color",[0 0 0]);
p3u = (theta(2)^0.25+0.5*delTs).^4/theta(2);
p3l = (theta(2)^0.25-0.5*delTs).^4/theta(2);

plot3(ax,[0 5],[1 1],[1 1]*max(IA2E(:)),"-","Color",[0 0 0],"LineWidth",1.5);
plot3(ax,delTs,p3u,ones(size(delTs))*max(IA2E(:)),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax,delTs,p3l,ones(size(delTs))*max(IA2E(:)),"--","Color",[0 0 0],"LineWidth",1.5)

%surf(ax,delTs,Ls,IEE',"FaceColor",[0 0 0.8],"EdgeColor","none","FaceAlpha",0.5)
set(ax,"XLim",[0,max(delTs)-0.01],"YLim",[min(Ls)+0.1 max(Ls)-0.2],"FontSize",16)
xlabel(ax,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%%

close all
fig = figure("Position",[0 0 400 600]);
ax = axes(fig,"Position",[0.16 0.13 0.8 0.8]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
idx = f1+f2>0.99*max(f1(:)+f2(:));
plot(ax,Synergy(idx),dIE(idx),'.',"Color",ResourceGreen,"MarkerSize",6)
set(ax,"XLim",[-1,3.5],"YLim",[-0.05 0.55],"FontSize",16)
xlabel(ax,"Synergy, $I(A\rightarrow E)$","FontSize",20,"Interpreter","latex")
ylabel(ax,"$\Delta I_E/I_E$","FontSize",20,"Interpreter","latex")

%%
close all
fig = figure("Position",[0 0 400 600]);
ax = axes(fig,"Position",[0.16 0.13 0.8 0.8]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
idx = (f1+f2)>(0*max(f1(:)+f2(:)));
loglog(ax,IAE(idx),((f1(idx)+f2(idx))/f),'.',"Color",ResourceGreen,"MarkerSize",6)
set(ax,"XLim",[0.1,10],"YLim",[0.001 1],"FontSize",12,"XScale","log","YScale","log");
xlabel(ax,"$I(A:E)$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Viability","FontSize",20,"Interpreter","latex")

%%
close all
fig = figure("Position",[0 0 400 600]);
ax = axes(fig,"Position",[0.16 0.13 0.8 0.8]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
idx = (f1+f2)>(0.95*max(f1(:)+f2(:)));
loglog(ax,IAE(idx),Efficacy(idx),'.',"Color",ResourceGreen,"MarkerSize",6)
set(ax,"XLim",[0.1,10],"YLim",[0.001 1],"FontSize",12,"XScale","log","YScale","log");
xlabel(ax,"$I(A:E)$","FontSize",20,"Interpreter","latex")
ylabel(ax,"Viability","FontSize",20,"Interpreter","latex")

%% Dynamics 

f       =   0.88;
Ag      =   0.3;
Ab      =   0.1;
Aw      =   0.6;
Topt    =   300;
Q       =   0.1;
gG      =   1;
gD      =   0.2;
tauS    =   3;
tauE    =   5;
dT      =   60;
N       = 500;
dt      = 0.1;
ll = 90; % starting T 1-128
ii = 55; % bandwidth 1-128

theta0 = [ ...
    f , ...             %    1:     f
    1-Ag , ...          %    2:     1 - A_G
    Ab-Ag, ...          %    3:     A_B - A_G
    Aw-Ag, ...          %    4:     A_W - A_G
    Q, ...              %    5:     Q
    gD/gG, ...          %    6:     gamma_D / gamma_G
    1/(gG*tauE), ...    %    7:     1 / gamma_G tau_E
    1/(gG*tauS), ...    %    8:     1 / gamma_G tau_S
    8*(Topt/dT)^4, ...  %    9:     8 (T_opt / Delta T)^4
    0.05, ...           %    10:    delta
    0 ...               %    11:    lambda
    ];

dTs = linspace(0,80,129);
dTs(1) = [];
lambda = linspace(-0.7,1.4,400);
Ls = (1+lambda);
Ts = Ls.^0.25;

theta = theta0;
theta(9) = 8*(Topt/dTs(ii))^4;

theta(11) = lambda(ll);


f1 = rand()*f; f2 = 1 - f1;
y = rand();
f1 = y*f1; f2 = y*f2;

close all
fig = figure("Position",[0 0 700 600]);
ax = axes(fig);
p1 = plot(ax,0,f1,"-","Color",ResourceGreen,"LineWidth",1.5);
hold(ax,"on");box(ax,"on");
p2 = plot(ax,0,f2,"-","Color",DartmouthGreen,"LineWidth",1.5);
p3 = plot(ax,0,Ts(ll),"-","Color",[0 0 1],"LineWidth",1.5);
p4 = plot(ax,0,Ls(ll),"-","Color",[1 0 0],"LineWidth",1.5);

set(ax,"XLim",[0,dt])
drawnow

x = [f1 f2 Ts(ll) Ls(ll)];
for tt = 1:30000
    x = updateExoDaisyWorld(x,dt,theta);
    p1.XData = [p1.XData p1.XData(end)+dt];
    p1.YData = [p1.YData x(1)];
    p2.XData = [p2.XData p2.XData(end)+dt];
    p2.YData = [p2.YData x(2)];
    p3.XData = [p3.XData p3.XData(end)+dt];
    p3.YData = [p3.YData x(3)];
    p4.XData = [p4.XData p4.XData(end)+dt];
    p4.YData = [p4.YData x(4)];
    set(ax,"XLim", [0 tt*dt])
    if mod(tt,10)==0
        drawnow
        pause(eps)
    end
end

