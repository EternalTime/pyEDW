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

%% FIG1: Compilation of Daisy Worlds wiht different Bandwidths
limT0 = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899];

fnums = [8 30 100];

load("data_000.mat");
Lenv = data(:,:,4);
Tenv = data(:,:,3);


close all
dPlots = 0.01;
x0 = 0.06; y0 = 0.08;dy1 = 0.36;dy2 = 0.25;dy3 = 0.23;
dx1 = 0.2;dx2 = 0.28;dx3 = 0.39;
fig = figure("Position",[0 0 900 750],"Color",[1 1 1]);
ax11 = axes(fig,"Position",[x0 y0 dx1 dy1]);
hold(ax11,"on");box(ax11,"on");ax11.LineWidth = 1.5;
ax12 = axes(fig,"Position",[x0 y0+dy1+dPlots dx1 dy2]);
hold(ax12,"on");box(ax12,"on");ax12.LineWidth = 1.5;
ax13 = axes(fig,"Position",[x0 y0+dy1+dy2+2*dPlots dx1 dy3]);
hold(ax13,"on");box(ax13,"on");ax13.LineWidth = 1.5;

x0 = sum(ax11.Position([1 3]))+0.01*fig.Position(4)/fig.Position(3);
ax21 = axes(fig,"Position",[x0 y0 dx2 dy1]);
hold(ax21,"on");box(ax21,"on");ax21.LineWidth = 1.5;
ax22 = axes(fig,"Position",[x0 y0+dy1+dPlots dx2 dy2]);
hold(ax22,"on");box(ax22,"on");ax22.LineWidth = 1.5;
ax23 = axes(fig,"Position",[x0 y0+dy1+dy2+2*dPlots dx2 dy3]);
hold(ax23,"on");box(ax23,"on");ax23.LineWidth = 1.5;

x0 = sum(ax21.Position([1 3]))+0.01*fig.Position(4)/fig.Position(3);
ax31 = axes(fig,"Position",[x0 y0 dx3 dy1]);
hold(ax31,"on");box(ax31,"on");ax31.LineWidth = 1.5;
ax32 = axes(fig,"Position",[x0 y0+dy1+dPlots dx3 dy2]);
hold(ax32,"on");box(ax32,"on");ax32.LineWidth = 1.5;
ax33 = axes(fig,"Position",[x0 y0+dy1+dy2+2*dPlots dx3 dy3]);
hold(ax33,"on");box(ax33,"on");ax33.LineWidth = 1.5;

getA=  @(x,y) "ax"+num2str(x)+num2str(y);
for ff = 1:length(fnums)
    load("data_"+num2str(fnums(ff),"%03.0f")+".mat")
    L = data(:,:,4);
    T = data(:,:,3);
    f1 = data(:,:,1);
    f2 = data(:,:,2);
    A = 1-theta(2)+theta(3)*f1+theta(4)*f2;

    limT = limT0/(theta(9))^0.25;
    for ii = 1:length(limT)
        fill(eval(getA(ff,1)),[0.1 3 3 0.1],1+[-1 -1 1 1]*limT(ii),...
            ResourceGreen,"FaceColor","none","FaceAlpha",0.01,...
            "EdgeColor",[0 0 0],"LineWidth",0.1,"EdgeAlpha",0.3);
    end
    plot(eval(getA(ff,1)),[0.1 3],[1 1],"--","Color",[0 0 0],"LineWidth",1.0)
    plot(eval(getA(ff,1)),[1 1],[0.1 10],"--","Color",[0 0 0],"LineWidth",1.0)
    plot(eval(getA(ff,2)),[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)
    plot(eval(getA(ff,3)),[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)

    kernel = [1 2 2 2 1];kernel = kernel/sum(kernel);
    Lm = conv(mean(L,2),kernel,'valid');
    Tm = conv(mean(T,2),kernel,'valid');
    Am = conv(mean(A,2),kernel,"valid");
    f1m = conv(mean(f1,2),kernel,"valid");
    f2m = conv(mean(f2,2),kernel,"valid");
    ftot = conv(sum(f1 + f2,2)/500,kernel,"valid");
    
    scatter(eval(getA(ff,1)),L(:),T(:),"o","MarkerFaceColor",0.8*ResourceGreen,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
    plot(eval(getA(ff,1)),Lm,Tm,"-.","Color",ResourceGreen,"LineWidth",2.0)
    plot(eval(getA(ff,1)),mean(Lenv,2),mean(Tenv,2),"-","Color",[0 0 0],"LineWidth",2.0)

    scatter(eval(getA(ff,2)),L(:),A(:),"o","MarkerFaceColor",[0 0 0],"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
    plot(eval(getA(ff,2)),Lm,Am,"-.","Color",[0 0 0],"LineWidth",2.0)
    
    plot(eval(getA(ff,3)),Lm,ftot,"--","Color",[1 1 1]*0.2,"LineWidth",0.5)
    scatter(eval(getA(ff,3)),L(:),f1(:),"o","MarkerFaceColor",RochesterBlue,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
    scatter(eval(getA(ff,3)),L(:),f2(:),"o","MarkerFaceColor",ForagerPink,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
    plot(eval(getA(ff,3)),Lm,f1m,"-.","Color",RochesterBlue,"LineWidth",2.0)
    plot(eval(getA(ff,3)),Lm,f2m,"-.","Color",ForagerPink,"LineWidth",2.0)
    plot(eval(getA(ff,3)),[0.1 3],[1 1]*theta(1),"--","Color",[0 0 0],"LineWidth",0.5)
    
end
load("data_"+num2str(fnums(1),"%03.0f")+".mat","theta");
set(ax11,"XLim",[0.75 1.275],"YLim",[0.9 1.1],...
    "YTick",[1-limT0(18)/theta(9)^0.25 1 1+limT0(18)/theta(9)^0.25],...
    "YTickLabel",[],...["" "$T_{opt}$" ""],...
    "FontSize",20,...
    "XTick",[1 ], ...
    "XTickLabel",[ "$L_{opt}$" ])
set(ax12,"XLim",[0.75 1.275],"YLim",[0.1 0.6],...
    "XTick",[],...
    "YTick",[Ag],...
    "YTickLabel",[],..."$A_G$",...
    "FontSize",20)
set(ax13,"XLim",[0.75 1.275],"YLim",[0 1],...
    "XTick",[],...
    "YTick",[f],...
    "YTickLabel",[],..."$f$",...
    "FontSize",20)

load("data_"+num2str(fnums(2),"%03.0f")+".mat","theta");
set(ax21,"XLim",[0.65 1.525],"YLim",[0.88 1.12],...
    "YTick",[1-limT0(18)/theta(9)^0.25 1 1+limT0(18)/theta(9)^0.25],...
    "YTickLabel",[],...["" "$T_{opt}$" ""],...
    "FontSize",20,...
    "XTick",[1 ], ...
    "XTickLabel",[]);%[ "$L_{opt}$" ])
set(ax22,"XLim",[0.65 1.525],"YLim",[0.1 0.6],...
    "XTick",[],...
    "YTick",[Ag],...
    "YTickLabel",[],..."$A_G$",...
    "FontSize",20);
set(ax23,"XLim",[0.65 1.525],"YLim",[0 1],...
    "XTick",[],...
    "YTick",[f],...
    "YTickLabel",[],... "$f$",...
    "FontSize",20);

load("data_"+num2str(fnums(3),"%03.0f")+".mat","theta");
set(ax31,"XLim",[0.4 2.05],"YLim",[0.8 1.2],...
    "YTick",[1-limT0(18)/theta(9)^0.25 1 1+limT0(18)/theta(9)^0.25],...
    "YTickLabel",["" "$T_{opt}$" ""],...
    "FontSize",20,...
    "XTick",1 , ...
    "XTickLabel", "$L_{opt}$" ,...
    "YAxisLocation","right")
set(ax32,"XLim",[0.4 2.05],"YLim",[0.1 0.6],...
    "XTick",[],...
    "YTick",[Ag],...
    "YTickLabel","$A_G$",...
    "FontSize",20,...
    "YAxisLocation","right")
set(ax33,"XLim",[0.4 2.05],"YLim",[0 1],...
    "XTick",[],...
    "YTick",[f],...
    "YTickLabel","$f$",...
    "FontSize",20,...
    "YAxisLocation","right")

ax13.Title.String = "$\Delta T/T_{opt}="+num2str(delTs(fnums(1)),"%0.03f")+"$";
ax23.Title.String = "$\Delta T/T_{opt}="+num2str(delTs(fnums(2)),"%0.03f")+"$";
ax33.Title.String = "$\Delta T/T_{opt}="+num2str(delTs(fnums(3)),"%0.03f")+"$";

xlabel(ax21,"Stellar Luminosity","FontSize",25)
ylabel(ax11,"Planetary Temperature","FontSize",25)
ylabel(ax12,"Albedo","FontSize",25)
ylabel(ax13,"Daisy Fraction","FontSize",25)
%% Figure 1: Daisyworld effect
% conditioned over average luminosities 
% data_010

load("data_000.mat")


limT0 = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899]/theta0(9)^0.25;
load("data_008.mat")
L = data(:,:,4);
T = data(:,:,3);
f1 = data(:,:,1);
f2 = data(:,:,2);
A = 1-theta(2)+theta(3)*f1+theta(4)*f2;
limT = limT0*(theta0(9)/theta(9))^0.25;

close all
fig = figure("Position",[0 0 275 750],"Color",[1 1 1]);
ax = axes(fig,"Position",[0.13 0.1 0.82 0.36]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
ax2 = axes(fig,"Position",[0.13 0.75 0.82 0.23]);
hold(ax2,"on");box(ax2,"on");ax2.LineWidth = 1.5;
ax3 = axes(fig,"Position",[0.13 0.48 0.82 0.25]);
hold(ax3,"on");box(ax3,"on");ax3.LineWidth = 1.5;

% optimal temperature strip plot
for ii = 1:length(limT)
    fill(ax,[0.1 3 3 0.1],1+[-1 -1 1 1]*limT(ii),...
        ResourceGreen,"FaceColor",ResourceGreen,"FaceAlpha",0.01,...
        "EdgeColor",ResourceGreen,"LineWidth",0.1,"EdgeAlpha",0.3)
end
plot(ax,[0.1 3],[1 1],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax,[1 1],[0.1 10],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax2,[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax3,[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)

% Daisy effect on temperature plot

kernel = [1 2 2 2 1];kernel = kernel/sum(kernel);
Lm = conv(mean(L,2),kernel,'valid');
Tm = conv(mean(T,2),kernel,'valid');
Am = conv(mean(A,2),kernel,"valid");
f1m = conv(mean(f1,2),kernel,"valid");
f2m = conv(mean(f2,2),kernel,"valid");
ftot = conv(sum(f1 + f2,2)/500,kernel,"valid");

scatter(ax,L(:),T(:),"o","MarkerFaceColor",DartmouthGreen,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax,Lm,Tm,"-.","Color",ResourceGreen,"LineWidth",2.0)
plot(ax2,Lm,ftot,"--","Color",[1 1 1]*0.2,"LineWidth",0.5)
scatter(ax2,L(:),f1(:),"o","MarkerFaceColor",RochesterBlue,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
scatter(ax2,L(:),f2(:),"o","MarkerFaceColor",ForagerPink,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax2,Lm,f1m,"-.","Color",RochesterBlue,"LineWidth",2.0)
plot(ax2,Lm,f2m,"-.","Color",ForagerPink,"LineWidth",2.0)
plot(ax2,[0.1 3],[1 1]*theta(1),"--","Color",[0 0 0],"LineWidth",0.5)
scatter(ax3,L(:),A(:),"o","MarkerFaceColor",[0 0 0],"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax3,Lm,Am,"-.","Color",[0 0 0],"LineWidth",2.0)





annotation(fig,'textarrow', [0.9 0.9], [0.475 0.427]-0.11)
annotation(fig,'textarrow', [0.9 0.9], [0.235 0.282]-0.04,...
    "String","$\Delta T$","Interpreter","latex",...
    "FontSize",16)

set(ax,"XLim",[0.75 1.375],"YLim",[0.85 1.15],...
    "YTick",[1-limT(18) 1 1+limT(18)],...
    "YTickLabel",["" "$T_{opt}$" ""],...
    "FontSize",18,...
    "XTick",[1 ], ...
    "XTickLabel",[ "$L_{opt}$" ])
set(ax2,"XLim",[0.75 1.375],"YLim",[0 1],...
    "XTick",[],...
    "YTick",[f],...
    "YTickLabel","$f$",...
    "FontSize",18)
set(ax3,"XLim",[0.75 1.375],"YLim",[0.1 0.6],...
    "XTick",[],...
    "YTick",[Ag],...
    "YTickLabel","$A_G$",...
    "FontSize",18)

xlabel(ax,"Stellar Luminosity","Interpreter","latex","FontSize",20)
%ylabel(ax,"Planetary Temperature","Interpreter","latex","FontSize",20)
%ylabel(ax2,"Daisy Fractions","Interpreter","latex","FontSize",20)
%ylabel(ax3,"Planetary Albedo","Interpreter","latex","FontSize",20)

%% Figure 1: Daisyworld effect
% conditioned over average luminosities 

limT0 = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899]/theta0(9)^0.25;
load("data_010.mat")
L = data(:,:,4);
T = data(:,:,3);
f1 = data(:,:,1);
f2 = data(:,:,2);
A = 1-theta(2)+theta(3)*f1+theta(4)*f2;
limT = limT0*(theta0(9)/theta(9))^0.25;

close all
fig = figure("Position",[0 0 500 750],"Color",[1 1 1]);
ax = axes(fig,"Position",[0.13 0.1 0.82 0.36]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
ax2 = axes(fig,"Position",[0.13 0.75 0.82 0.23]);
hold(ax2,"on");box(ax2,"on");ax2.LineWidth = 1.5;
ax3 = axes(fig,"Position",[0.13 0.48 0.82 0.25]);
hold(ax3,"on");box(ax3,"on");ax3.LineWidth = 1.5;

% optimal temperature strip plot
for ii = 1:length(limT)
    fill(ax,[0.1 3 3 0.1],1+[-1 -1 1 1]*limT(ii),...
        ResourceGreen,"FaceColor",ResourceGreen,"FaceAlpha",0.01,...
        "EdgeColor",ResourceGreen,"LineWidth",0.1,"EdgeAlpha",0.3)
end
plot(ax,[0.1 3],[1 1],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax,[1 1],[0.1 10],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax2,[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax3,[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)

% Daisy effect on temperature plot

kernel = [1 2 2 2 1];kernel = kernel/sum(kernel);
Lm = conv(mean(L,2),kernel,'valid');
Tm = conv(mean(T,2),kernel,'valid');
Am = conv(mean(A,2),kernel,"valid");
f1m = conv(mean(f1,2),kernel,"valid");
f2m = conv(mean(f2,2),kernel,"valid");
ftot = conv(sum(f1 + f2,2)/500,kernel,"valid");

scatter(ax,L(:),T(:),"o","MarkerFaceColor",DartmouthGreen,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax,Lm,Tm,"-.","Color",ResourceGreen,"LineWidth",2.0)
plot(ax2,Lm,ftot,"--","Color",[1 1 1]*0.2,"LineWidth",0.5)
scatter(ax2,L(:),f1(:),"o","MarkerFaceColor",RochesterBlue,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
scatter(ax2,L(:),f2(:),"o","MarkerFaceColor",ForagerPink,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax2,Lm,f1m,"-.","Color",RochesterBlue,"LineWidth",2.0)
plot(ax2,Lm,f2m,"-.","Color",ForagerPink,"LineWidth",2.0)
plot(ax2,[0.1 3],[1 1]*theta(1),"--","Color",[0 0 0],"LineWidth",0.5)
scatter(ax3,L(:),A(:),"o","MarkerFaceColor",[0 0 0],"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax3,Lm,Am,"-.","Color",[0 0 0],"LineWidth",2.0)



% Agent Free temperature plot
load("data_000.mat")
L = data(:,:,4);
T = data(:,:,3);
plot(ax,mean(L,2),mean(T,2),"-","Color",[0 0 0],"LineWidth",2.0)


annotation(fig,'textarrow', [0.9 0.9], [0.475 0.427]-0.11)
annotation(fig,'textarrow', [0.9 0.9], [0.235 0.282]-0.04,...
    "String","$\Delta T$","Interpreter","latex",...
    "FontSize",16)

set(ax,"XLim",[0.75 1.375],"YLim",[0.85 1.15],...
    "YTick",[1-limT(18) 1 1+limT(18)],...
    "YTickLabel",["" "$T_{opt}$" ""],...
    "FontSize",18,...
    "XTick",[1 ], ...
    "XTickLabel",[ "$L_{opt}$" ])
set(ax2,"XLim",[0.75 1.375],"YLim",[0 1],...
    "XTick",[],...
    "YTick",[f],...
    "YTickLabel","$f$",...
    "FontSize",18)
set(ax3,"XLim",[0.75 1.375],"YLim",[0.1 0.6],...
    "XTick",[],...
    "YTick",[Ag],...
    "YTickLabel","$A_G$",...
    "FontSize",18)

xlabel(ax,"Stellar Luminosity","Interpreter","latex","FontSize",20)
ylabel(ax,"Planetary Temperature","Interpreter","latex","FontSize",20)
ylabel(ax2,"Daisy Fractions","Interpreter","latex","FontSize",20)
ylabel(ax3,"Planetary Albedo","Interpreter","latex","FontSize",20)
%% Figure 1: Daisyworld effect
% conditioned over average luminosities 

limT0 = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899]/theta0(9)^0.25;
load("data_010.mat")
L = data(:,:,4);
T = data(:,:,3);
f1 = data(:,:,1);
f2 = data(:,:,2);
A = 1-theta(2)+theta(3)*f1+theta(4)*f2;
limT = limT0*(theta0(9)/theta(9))^0.25;

close all
fig = figure("Position",[0 0 500 750],"Color",[1 1 1]);
ax = axes(fig,"Position",[0.13 0.1 0.82 0.36]);
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.5;
ax2 = axes(fig,"Position",[0.13 0.75 0.82 0.23]);
hold(ax2,"on");box(ax2,"on");ax2.LineWidth = 1.5;
ax3 = axes(fig,"Position",[0.13 0.48 0.82 0.25]);
hold(ax3,"on");box(ax3,"on");ax3.LineWidth = 1.5;

% optimal temperature strip plot
for ii = 1:length(limT)
    fill(ax,[0.1 3 3 0.1],1+[-1 -1 1 1]*limT(ii),...
        ResourceGreen,"FaceColor",ResourceGreen,"FaceAlpha",0.01,...
        "EdgeColor",ResourceGreen,"LineWidth",0.1,"EdgeAlpha",0.3)
end
plot(ax,[0.1 3],[1 1],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax,[1 1],[0.1 10],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax2,[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)
plot(ax3,[1 1],[0 10],"--","Color",[0 0 0],"LineWidth",1.0)

% Daisy effect on temperature plot

kernel = [1 2 2 2 1];kernel = kernel/sum(kernel);
Lm = conv(mean(L,2),kernel,'valid');
Tm = conv(mean(T,2),kernel,'valid');
Am = conv(mean(A,2),kernel,"valid");
f1m = conv(mean(f1,2),kernel,"valid");
f2m = conv(mean(f2,2),kernel,"valid");
ftot = conv(sum(f1 + f2,2)/500,kernel,"valid");

scatter(ax,L(:),T(:),"o","MarkerFaceColor",DartmouthGreen,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax,Lm,Tm,"-.","Color",ResourceGreen,"LineWidth",2.0)
plot(ax2,Lm,ftot,"--","Color",[1 1 1]*0.2,"LineWidth",0.5)
scatter(ax2,L(:),f1(:),"o","MarkerFaceColor",RochesterBlue,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
scatter(ax2,L(:),f2(:),"o","MarkerFaceColor",ForagerPink,"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax2,Lm,f1m,"-.","Color",RochesterBlue,"LineWidth",2.0)
plot(ax2,Lm,f2m,"-.","Color",ForagerPink,"LineWidth",2.0)
plot(ax2,[0.1 3],[1 1]*theta(1),"--","Color",[0 0 0],"LineWidth",0.5)
scatter(ax3,L(:),A(:),"o","MarkerFaceColor",[0 0 0],"MarkerEdgeColor","none","MarkerFaceAlpha",0.05,"SizeData",3)
plot(ax3,Lm,Am,"-.","Color",[0 0 0],"LineWidth",2.0)



% Agent Free temperature plot
load("data_000.mat")
L = data(:,:,4);
T = data(:,:,3);
plot(ax,mean(L,2),mean(T,2),"-","Color",[0 0 0],"LineWidth",2.0)


annotation(fig,'textarrow', [0.9 0.9], [0.475 0.427]-0.11)
annotation(fig,'textarrow', [0.9 0.9], [0.235 0.282]-0.04,...
    "String","$\Delta T$","Interpreter","latex",...
    "FontSize",16)

set(ax,"XLim",[0.5 1.75],"YLim",[0.85 1.15],...
    "YTick",[1-limT(18) 1 1+limT(18)],...
    "YTickLabel",["" "$T_{opt}$" ""],...
    "FontSize",18,...
    "XTick",[1 ], ...
    "XTickLabel",[ "$L_{opt}$" ])
set(ax2,"XLim",[0.5 1.75],"YLim",[0 1],...
    "XTick",[],...
    "YTick",[f],...
    "YTickLabel","$f$",...
    "FontSize",18)
set(ax3,"XLim",[0.5 1.75],"YLim",[0.1 0.6],...
    "XTick",[],...
    "YTick",[Ag],...
    "YTickLabel","$A_G$",...
    "FontSize",18)

xlabel(ax,"Stellar Luminosity","Interpreter","latex","FontSize",20)
ylabel(ax,"Planetary Temperature","Interpreter","latex","FontSize",20)
ylabel(ax2,"Daisy Fractions","Interpreter","latex","FontSize",20)
ylabel(ax3,"Planetary Albedo","Interpreter","latex","FontSize",20)