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

c1 = repmat(fliplr(linspace(0,1,40))',[1 3]);
c2 = [flipud(c1);c1;zeros(40,3)];
c3 = flipud(c2);
c1 = [c1;zeros(80,3)];
c4 = flipud(c1);

pink2   = 0.9*([222 215 217]-1)/255;
blue2   = 0.9*([181 195 210]-1)/255;
green2  = 0.9*([184 209 198]-1)/255;

% %myCmap = repmat([1 1 1],[100,1]).*cx+repmat(pink2,[100,1]).*cy +...
%     repmat(blue2,[100 1]).*cz;
% 
% myCmap = repmat([1 1 1],[100,1]).*cx+repmat(pink2,[100,1]).*cy +...
%     repmat(blue2,[100 1]).*cz;

myCmap = repmat([1 1 1],[120,1]).*c1 + ...
    repmat(pink2,[120,1]).*c2 +...
    repmat(blue2,[120 1]).*c3 + ...
    repmat(green2,[120 1]).*c4;

myCmap2 = repmat(green2,[120,1]).*c1 + ...
    repmat([1 1 1],[120,1]).*c2 +...
    repmat(pink2,[120 1]).*c3 + ...
    repmat(blue2,[120 1]).*c4;

myCmap3 = repmat(green2,[120,1]).*c1 + ...
    repmat(blue2,[120,1]).*c2 +...
    repmat(pink2,[120 1]).*c3 + ...
    repmat([1 1 1],[120 1]).*c4;

%% FIG3: Viability contours and Efficacy
theta = theta0;

load("data_temperaturesV.mat");
Efficacy = aT-1;
Viability = (f1+f2)/f;
mask = Viability<0.02;

kernel = [1 2 1]'; kernel = kernel/sum(kernel(:));
for kk = 1:1
    Viability = conv2(padarray(Viability,[1 0],"replicate"),kernel,"valid");
    Efficacy = conv2(padarray(Efficacy,[1 0],"replicate"),kernel,"valid");
end
Efficacy(mask) = nan;
Viability(mask) = nan;
close all 
fig = figure("Position",[0  0 600,600]);
ax1 = axes(fig,"Position",[0.11 0.15 0.75 0.76]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;

im1 = surf(ax1,delTs,Ls(1:(end)),Efficacy(1:(end),:),"EdgeColor","none","FaceAlpha",1);


cb = colorbar(ax1,"Position",[sum(ax1.Position([1 3]))+0.02 ax1.Position(2) 0.025 ax1.Position(4)],...
        "LineWidth", 1.5,"TickLabelInterpreter","latex","FontSize",12,"TickLabelsMode","auto");
cb.Label.String = "Efficacy";
cb.Label.FontSize = 18;
cb.Label.Interpreter = 'latex';
cb.Ruler.TickLabelRotation = 45;

cMax = max(Impact(:));
cMin = min(Impact(:));
cStart = ceil(cMax/abs(cMin)*100);
cmap = myCmap;%flipud(myCmap);
cmap = [flipud(cmap);cmap(1:cStart,:)];
colormap(ax1,cmap);

set(ax1,"CLim",[cMin cMax]);
cts = sort([1.042,(1+linspace(min(Viability(:)),max(Viability(:)),12))]*max(abs(Efficacy(:))));
cts(end) = [];

[~,ct] = contour3(ax1,delTs,Ls,(1+Viability)*max(abs(Efficacy(:))),cts,"LineWidth",1.5,"Color",[0 0 0]+0.1);     

plot3(ax1,delTs,(1+0.5*delTs).^4,ones(size(delTs))*max(abs(Impact(:))),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax1,delTs,(1-0.5*delTs).^4,ones(size(delTs))*max(abs(Impact(:))),"--","Color",[0 0 0],"LineWidth",1.5)
set(ax1,"YDir","normal","XLim",[0,max(delTs)-0.01],...
    "YLim",[min(Ls)+0.1 max(Ls)-0.2],"LineWidth",1.5,"FontSize",16)
xlabel(ax1,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax1,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%%
%% Testing FIg: Cooperation contours Drop in environmental correlation
theta = theta0;


load("data_entropyE.mat")
Iee0 = repmat(H(:,1) + H(:,2) - H(:,3),[1 128]);

load("data_temperaturesV.mat");
Viability = (f1+f2)/f;
Efficacy = aT-1;

%mask = (f1<10^-5)&(f2<10^-5);
mask = Viability < 0.01;
idx = mask;mask = 1.0*mask;
mask(idx) = nan;
mask(not(idx)) = inf;
for kk = 1:3
    mask = conv2(mask,ones(2),"same");
end
mask = isnan(mask);


load("data_entropy3.mat")
H = permute(H,[2 1 3]);

% Entropy H( Bandwidth, Luminosity, #) 
% 1:    H[a1]
% 2:    H[a2]
% 3:    H[e1]
% 4:    H[e2]
% 5:    H[e1e2]
% 6:    H[a2e2]
% 7:    H[a2e1]
% 8:    H[a1e2]
% 9:    H[a1e3]
% 10:   H[a1a2]
% 11:   H[a2e1e2]
% 12:   H[a1e1e2]
% 13:   H[a1a2e2]
% 14:   H[a1a2e1]
% 15:   H[a1a2e1e2]
% 16:   number of bins

% Interaction info I(a1:a2:E) = I(a1:e1e2) + I(a2:e1e2) - I(a1a2:e1e2)
%                             =  H(a1a2e1e2) - H(a1a2)  - H(a1e1e2) - H(a2e1e2) + H(a1) + H(a2) + H(e1e2)
%                              
%   sum(H(:,:,[]),3) - sum(H(:,:,[]),3);

Iee     = sum(H(:,:,[3 4]),3) - H(:,:,5);

dIE     = Iee - Iee0;
Coop    = -sum(H(:,:,[10 11 12]),3) + sum(H(:,:,[1 2 5 15]),3);

kernel1 = [1 2 1]'; kernel1 = kernel1/sum(kernel1(:));
kernel2 = [1 2 1]; kernel2 = kernel2/sum(kernel2(:));
for kk = 1:10
    Viability = conv2(padarray(Viability,[1 0],"replicate"),kernel1,"valid");
    Efficacy = conv2(padarray(Efficacy,[1 0],"replicate"),kernel1,"valid");
    
    dIE = conv2(padarray(dIE,[1 0],"replicate"),kernel1,"valid");
    Coop = conv2(padarray(Coop,[1 0],"replicate"),kernel1,"valid");
end
for kk = 1:1
    Viability = conv2(padarray(Viability,[0 1],"replicate"),kernel2,"valid");
    Efficacy = conv2(padarray(Efficacy,[0 1],"replicate"),kernel2,"valid");
    dIE = conv2(padarray(dIE,[0 1],"replicate"),kernel2,"valid");
    Coop = conv2(padarray(Coop,[0 1],"replicate"),kernel2,"valid");
end

Viability(mask) = nan;
Efficacy(mask) = nan;
dIE(mask) = nan;
Coop(mask) = nan;

close all 
fig = figure("Position",[0  0 600,600]);
ax1 = axes(fig,"Position",[0.11 0.15 0.75 0.76]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;

im1 = surf(ax1,delTs,Ls(1:(end)),Coop(1:(end),:),"EdgeColor","none","FaceAlpha",1);


cb = colorbar(ax1,"Position",[sum(ax1.Position([1 3]))+0.02 ax1.Position(2) 0.025 ax1.Position(4)],...
        "LineWidth", 1.5,"TickLabelInterpreter","latex","FontSize",12,"TickLabelsMode","auto");
cb.Label.String = "Efficacy";
cb.Label.FontSize = 18;
cb.Label.Interpreter = 'latex';
cb.Ruler.TickLabelRotation = 45;

cMax = max(Coop(:));
cMin = min(Coop(:));
cStart = ceil(abs(cMin)/abs(cMax)*100);
cmap = myCmap2;%flipud(myCmap);
%cmap = [flipud(cmap(1:cStart,:));cmap];
colormap(ax1,cmap);

set(ax1,"CLim",[cMin cMax]);
cts = sort((2+linspace(min(dIE(:)),max(dIE(:)),11))*max(abs(Coop(:))));
cts(end) = [];

[~,ct] = contour3(ax1,delTs,Ls,(2+dIE)*max(abs(Coop(:))),cts,"LineWidth",1.5,"Color",[0 0 0]+0.1);     

plot3(ax1,delTs,(1+0.5*delTs).^4,ones(size(delTs))*max(abs(Coop(:))),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax1,delTs,(1-0.5*delTs).^4,ones(size(delTs))*max(abs(Coop(:))),"--","Color",[0 0 0],"LineWidth",1.5)
set(ax1,"YDir","normal","XLim",[0,max(delTs)-0.01],...
    "YLim",[min(Ls)+0.1 max(Ls)-0.2],"LineWidth",1.5,"FontSize",16)
xlabel(ax1,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax1,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")

%% Testing FIg: Cooperation contours Drop in environmental correlation
theta = theta0;


load("data_entropyE.mat")
Iee0 = repmat(H(:,1) + H(:,2) - H(:,3),[1 128]);

load("data_temperaturesV.mat");
Viability = (f1+f2)/f;
Efficacy = aT-1;

%mask = (f1<10^-5)&(f2<10^-5);
mask = Viability < 0.0001;
idx = mask;mask = 1.0*mask;
mask(idx) = nan;
mask(not(idx)) = inf;
for kk = 1:3
    mask = conv2(mask,ones(2),"same");
end
mask = isnan(mask);


load("data_entropy3.mat")
H = permute(H,[2 1 3]);

% Entropy H( Bandwidth, Luminosity, #) 
% 1:    H[a1]
% 2:    H[a2]
% 3:    H[e1]
% 4:    H[e2]
% 5:    H[e1e2]
% 6:    H[a2e2]
% 7:    H[a2e1]
% 8:    H[a1e2]
% 9:    H[a1e3]
% 10:   H[a1a2]
% 11:   H[a2e1e2]
% 12:   H[a1e1e2]
% 13:   H[a1a2e2]
% 14:   H[a1a2e1]
% 15:   H[a1a2e1e2]
% 16:   number of bins

% Interaction info I(a1:a2:E) = I(a1:e1e2) + I(a2:e1e2) - I(a1a2:e1e2)
%                             =  H(a1a2e1e2) - H(a1a2)  - H(a1e1e2) - H(a2e1e2) + H(a1) + H(a2) + H(e1e2)
%                              
%   sum(H(:,:,[]),3) - sum(H(:,:,[]),3);

Iee     = sum(H(:,:,[3 4]),3) - H(:,:,5);

dIE     = Iee - Iee0;
Coop    = -sum(H(:,:,[10 11 12]),3) + sum(H(:,:,[1 2 5 15]),3);

kernel1 = [1 2 1]'; kernel1 = kernel1/sum(kernel1(:));
kernel2 = [1 2 1]; kernel2 = kernel2/sum(kernel2(:));
for kk = 1:7
    Viability = conv2(padarray(Viability,[1 0],"replicate"),kernel1,"valid");
    Efficacy = conv2(padarray(Efficacy,[1 0],"replicate"),kernel1,"valid");
    
    dIE = conv2(padarray(dIE,[1 0],"replicate"),kernel1,"valid");
    Coop = conv2(padarray(Coop,[1 0],"replicate"),kernel1,"valid");
end
for kk = 1:2
    Viability = conv2(padarray(Viability,[0 1],"replicate"),kernel2,"valid");
    Efficacy = conv2(padarray(Efficacy,[0 1],"replicate"),kernel2,"valid");
    dIE = conv2(padarray(dIE,[0 1],"replicate"),kernel2,"valid");
    Coop = conv2(padarray(Coop,[0 1],"replicate"),kernel2,"valid");
end

Viability(mask) = nan;
Efficacy(mask) = nan;
dIE(mask) = nan;
Coop(mask) = nan;

close all 
fig = figure("Position",[0  0 600,600]);
ax1 = axes(fig,"Position",[0.11 0.15 0.75 0.76]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;

im1 = surf(ax1,delTs,Ls(1:(end)),dIE(1:(end),:),"EdgeColor","none","FaceAlpha",1);


cb = colorbar(ax1,"Position",[sum(ax1.Position([1 3]))+0.02 ax1.Position(2) 0.025 ax1.Position(4)],...
        "LineWidth", 1.5,"TickLabelInterpreter","latex","FontSize",12,"TickLabelsMode","auto");
cb.Label.String = "$\Delta I_{\phi\rightarrow A}$";
cb.Label.FontSize = 18;
cb.Label.Interpreter = 'latex';
cb.Ruler.TickLabelRotation = 45;

cMax = max(dIE(:));
cMin = min(dIE(:));
cStart = ceil(abs(cMin)/abs(cMax)*100);
cmap = myCmap3;%flipud(myCmap);
%cmap = [flipud(cmap(1:cStart,:));cmap];
colormap(ax1,cmap);

set(ax1,"CLim",[cMin cMax]);
cts = sort((2+linspace(min(Coop(:)),max(Coop(:)),17))*max(abs(dIE(:))));
cts(end) = [];

idxP = (cts-2*max(abs(dIE(:))))>0;
[~,ct] = contour3(ax1,delTs,Ls,(2+Coop)*max(abs(dIE(:))),cts(idxP),"LineWidth",1.5,"Color",[0 0 0]+0.1);  
contour3(ax1,delTs,Ls,(2+Coop)*max(abs(dIE(:))),cts(not(idxP)),"LineWidth",1.5,"Color",[0 0 0]+0.1,"LineStyle","-.");  

plot3(ax1,delTs,(1+0.5*delTs).^4,ones(size(delTs))*max(abs(dIE(:))),"--","Color",[0 0 0],"LineWidth",1.5)
plot3(ax1,delTs,(1-0.5*delTs).^4,ones(size(delTs))*max(abs(dIE(:))),"--","Color",[0 0 0],"LineWidth",1.5)
set(ax1,"YDir","normal","XLim",[0,max(delTs)-0.01],...
    "YLim",[min(Ls)+0.1 max(Ls)-0.2],"LineWidth",1.5,"FontSize",16)
xlabel(ax1,"Bandwidth, $\Delta T/T_{opt}$","FontSize",20,"Interpreter","latex")
ylabel(ax1,"Luminosity, $L/L_{opt}$","FontSize",20,"Interpreter","latex")
