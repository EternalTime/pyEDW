%% I(E:E) drop due to Biome

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

myCmap = repmat([1 1 1],[120,1]).*c1 + ...
    repmat(pink2,[120,1]).*c2 +...
    repmat(blue2,[120 1]).*c3 + ...
    repmat(green2,[120 1]).*c4;

%%
load("data_entropyE.mat")
Iee0 = repmat(H(:,1) + H(:,2) - H(:,3),[1 128]);

load("data_temperaturesV.mat");
p1 = f1./(f1+f2);
p2 = f2./(f1+f2);

Viability = (f1+f2)/f;
Efficacy = aT-1;

mask = (f1<10^-5)&(f2<10^-5);
mask = Viability < 0.0001;
idx = mask;mask = 1.0*mask;
mask(idx) = nan;
mask(not(idx)) = inf;
for kk = 1:2
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

kernel = [1 2 1]'; kernel = kernel/sum(kernel(:));
for kk = 1:0
    Viability = conv2(padarray(Viability,[1 0],"replicate"),kernel,"valid");
    Efficacy = conv2(padarray(Efficacy,[1 0],"replicate"),kernel,"valid");
    
    dIE = conv2(padarray(dIE,[1 0],"replicate"),kernel,"valid");
    Coop = conv2(padarray(Coop,[1 0],"replicate"),kernel,"valid");
end

kernel = [1 2 1]; kernel = kernel/sum(kernel(:));
for kk = 1:0
    Viability = conv2(padarray(Viability,[0 1],"replicate"),kernel,"valid");
    Efficacy = conv2(padarray(Efficacy,[0 1],"replicate"),kernel,"valid");
    dIE = conv2(padarray(dIE,[0 1],"replicate"),kernel,"valid");
    Coop = conv2(padarray(Coop,[0 1],"replicate"),kernel,"valid");
end

Viability(mask) = nan;
Efficacy(mask) = nan;
dIE(mask) = nan;
Coop(mask) = nan;

%%

close all
plot3(Coop,dIE,Viability,'.')
%% Figure: cooperation vs correlation drop heat map

close all
fig = figure("Position",[0  0 600,600]);
ax1 = axes(fig,"Position",[0.12 0.12 0.8 0.8]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;

idx1 = Viability>0.75;
idx2 = Viability<0.3;

[Z1,X1,Y1] = histcounts2(Coop(idx1),dIE(idx1),linspace(-1.1,3.6,103),linspace(-0.8,0.2,101));
X1=conv(X1,[0.5 0.5],"valid");
Y1=conv(Y1,[0.5 0.5],"valid");
kernel = [0 1 1 1 0; 1 1 2 1 1; 1 2 5 2 1; 1 1 2 1 1; 0 1 1 1 0];
kernel = kernel/sum(kernel(:));
for kk = 1:1
    Z1 = conv2(Z1,kernel,"same");
end
%imagesc(ax1,"XData",X1,"YData",Y1,"CData",Z1')
colormap(myCmap)
contourf(ax1,X1,Y1,Z1',10,"EdgeColor",[0 0 0],"FaceAlpha",1,"LineWidth",1.2)
hold(ax1,"on")
set(ax1,"XLim",[-1.2,3.4],"YLim",[-0.77,0.18],"FontSize",16)
xlabel(ax1,"$C(a_1:a_2||E)$","Interpreter","latex","FontSize",20)
ylabel(ax1,"$\Delta I_{\phi \rightarrow A}$","Interpreter","latex","FontSize",20)

%%
[Z1,X1,Y1] = histcounts2(Coop(idx2),dIE(idx2),linspace(-1.1,3.6,103),linspace(-0.8,0.2,101));
X1=conv(X1,[0.5 0.5],"valid");
Y1=conv(Y1,[0.5 0.5],"valid");
kernel = [0 1 1 1 0; 1 1 2 1 1; 1 2 5 2 1; 1 1 2 1 1; 0 1 1 1 0];
kernel = kernel/sum(kernel(:));
for kk = 1:3
    Z1 = conv2(Z1,kernel,"same");
end
contour(ax1,X1,Y1,Z1',10,"EdgeColor",[0 0 0],"FaceAlpha",0.5,"LineStyle","--","LineWidth",1.2)
%plot(ax1,Coop(idx1),dIE(idx1),'.')

%% Figure: efficacy vs viability heat map

close all
fig = figure("Position",[0  0 600,600]);
ax1 = axes(fig,"Position",[0.12 0.12 0.8 0.8]);
hold(ax1,"on");box(ax1,"on");ax1.LineWidth = 1.5;

idx1 = IAE>-3;
idx2 = IAE<0;

[Z1,X1,Y1] = histcounts2(Viability(idx1),(Efficacy(idx1)),linspace(-0.02,1,203),linspace(-0.25,0.25,201));
X1=conv(X1,[0.5 0.5],"valid");
Y1=conv(Y1,[0.5 0.5],"valid");
kernel = [0 1 1 1 0; 1 1 2 1 1; 1 2 5 2 1; 1 1 2 1 1; 0 1 1 1 0];
kernel = kernel/sum(kernel(:));

Z1 = log(Z1+eps);

for kk = 1:10
    Z1 = conv2(Z1,kernel,"same");
end
%imagesc(ax1,"XData",X1,"YData",Y1,"CData",Z1')
colormap(myCmap)
contourf(ax1,X1,Y1,Z1',20,"EdgeColor",[0 0 0],"FaceAlpha",1,"LineWidth",1.2)
set(ax1,"XLim",[0,0.85],"YLim",[-0.22,0.22],"FontSize",16)
hold(ax1,"on")

%[Z1,X1,Y1] = histcounts2(Viability(idx2),Efficacy(idx2),linspace(-0.02,1,203),linspace(-0.25,0.25,201));
%X1=conv(X1,[0.5 0.5],"valid");
%Y1=conv(Y1,[0.5 0.5],"valid");
%kernel = [0 1 1 1 0; 1 1 2 1 1; 1 2 5 2 1; 1 1 2 1 1; 0 1 1 1 0];
%kernel = kernel/sum(kernel(:));
%for kk = 1:3
%    Z1 = conv2(Z1,kernel,"same");
%end
%contour(ax1,X1,Y1,Z1',10,"EdgeColor",[0 0 0],"FaceAlpha",0.5,"LineStyle","--","LineWidth",1.2)
%plot(ax1,Coop(idx1),dIE(idx1),'.')
%set(ax1,"XLim",[-1.2,3.7],"YLim",[-0.9,0.3],"FontSize",16)
xlabel(ax1,"$V$","Interpreter","latex","FontSize",20)
ylabel(ax1,"$\mathcal{E}$","Interpreter","latex","FontSize",20)
