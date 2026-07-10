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
dT      =   2;

dt = 0.1;

theta = [ ...
            f , ...             %    1:     f
            1-Ag , ...          %    2:     1 - A_G
            Ab-Ag, ...          %    3:     A_B - A_G
            Aw-Ag, ...          %    4:     A_W - A_G
            Q, ...              %    5:     Q
            gD/gG, ...          %    6:     gamma_D / gamma_G
            1/(gG*tauE), ...    %    7:     1 / gamma_G tau_E
            1/(gG*tauS), ...    %    8:     1 / gamma_G tau_S
            8*(Topt/dT)^4, ...  %    9:     8 (T_opt / Delta T)^4
            0, ...              %    10:    delta
            0 ...               %    11:    lambda
        ];

dTs = linspace(0,80,129);
%% Plot Viability


maxIdx = 128;
close all
fig = figure("Position",[0 0 800 500]);
ax = axes(fig);
load("data_001.mat","Ls","Ts");
%Ls = sum(dataE(:,:,2),2);
%Ts = sum(dataE(:,:,1),2);
%plot(ax,Ls,Ts,"-","Color",[0 0 0],"LineWidth",2.0);
hold(ax,"on");box(ax,"on");
myCol = @(n) [1 0 0]+(n-1)/(maxIdx-1)*[-1 0 1];
X = zeros([1,maxIdx]);
Y = zeros([1,maxIdx]);
for ii = 1:maxIdx
    load("data_"+num2str(ii,"%03.0f")+".mat","data")
    L = squeeze(sum(data(:,:,4),2)/500)';
    T = squeeze(sum(data(:,:,3),2)/500)';
    eff = abs(T./Ts-1);
    plot(ax,dTs(ii),sum(eff)*mean(diff(L)),".","Color",myCol(ii),"MarkerSize",15)
    X(ii) = dTs(ii);
    Y(ii) = sum(eff)*mean(diff(L));
end
set(ax,"XScale","log","YScale","log")


%% Plot Viability

maxIdx = 128;
close all
fig = figure("Position",[0 0 800 500]);
ax = axes(fig);
load("data_001.mat","Ls","Ts");
%Ls = sum(dataE(:,:,2),2);
%Ts = sum(dataE(:,:,1),2);
%plot(ax,Ls,Ts,"-","Color",[0 0 0],"LineWidth",2.0);
hold(ax,"on");box(ax,"on");
myCol = @(n) [1 0 0]+(n-1)/(maxIdx-1)*[-1 0 1];
for ii = 1:maxIdx
    load("data_"+num2str(ii,"%03.0f")+".mat","data")
    L = squeeze(sum(data(:,:,4),2)/500)';
    T = squeeze(sum(data(:,:,3),2)/500)';
    plot(ax,L,abs(T./Ts-1),"Color",myCol(ii),"LineWidth",1.5)
end
set(ax,"XLim",[0.1 2.6],"YLim",[0.0 0.1])
%% Plot Efficacy

maxIdx = 128;
close all
fig = figure("Position",[0 0 800 500]);
ax = axes(fig);
hold(ax,"on");box(ax,"on");
myCol = @(n) [1 0 0]+(n-1)/(maxIdx-1)*[-1 0 1];
for ii = 1:maxIdx
    load("data_"+num2str(ii,"%03.0f")+".mat","data")
    L = squeeze(sum(data(:,:,4),2)/500)';
    T = squeeze(sum(data(:,:,3),2)/500)';
    plot(ax,L,T,".","Color",myCol(ii),"MarkerSize",15)
end
set(ax,"XLim",[0.2 2.4],"YLim",[0.0 3])
    

%% Get Environment Only Data

lambda = linspace(-0.7,1.4,400);
Ls = (1+lambda);
Ts = Ls.^0.25;

N = 500;
theta(10) = 0.05;

data = zeros([length(lambda),N,4]);
for ll = 1:length(Ls)
    tic
    x0 = [0 0 Ts(ll) Ls(ll)];
    theta(11) = lambda(ll);
    parfor nn = 1:N
        x = x0;
        for tt = 1:200
            x = updateExoDaisyWorld(x,dt,theta);
        end
        data(ll,nn,:) = x;
    end
    save("data_000.mat","data","theta","lambda","Ls","Ts")
    drawnow;pause(eps);
    disp("Lambda Run " + num2str(ll) +"/"+num2str(length(Ls))+...
        " complete in " + num2str(toc/60,"%.2f")+"min")
end


%% Test Environment and Daisies

close all
fig = figure("Position",[0 0 600 600]);
ax = axes(fig,"Position",[0.13 0.15 0.8 0.5]);
ax2 = axes(fig,"Position",[0.13 0.68 0.8 0.27]);

lambda = linspace(-0.7,1.4,200);
Ls = (1+lambda);
Ts = Ls.^0.25;

plot(ax,Ls,Ts,"--","Color",[0 0 0],"Linewidth", 2)
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.0;ax.FontSize = 16;
set(ax,"XLim",[min(Ls),max(Ls)])
limT = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899]/theta(9)^0.25;
for nn = 1:length(limT)
    fill(ax,[min(Ls) max(Ls) max(Ls) min(Ls)],1+[-1 -1 1 1]*limT(nn),...
        [0 0.4 0],"EdgeColor",[0 0 0],"FaceColor",[0 0.4 0],...
        "FaceAlpha",0.015,"LineWidth",0.1,"EdgeAlpha",0.3);
end
plot(ax2,[min(Ls) max(Ls)],[f f],"--","Color",[0 0 0],"LineWidth",0.5)
hold(ax2,"on");box(ax2,"on");ax2.LineWidth = 1.0;ax2.FontSize = 16;
set(ax2,"XLim",[min(Ls),max(Ls)],"XTickLabel",[],"YLim",[0 1])

plt = plot(ax,Ls(1),Ts(1),".","Color",[0 0 0],"MarkerSize",4);
pltB = plot(ax2,Ls(1),0,".","Color",newBlue,"MarkerSize",4);
pltW = plot(ax2,Ls(1),0,".","Color",ForagerPink,"MarkerSize",4);
pltD = plot(ax2,Ls(1),0,"-","Color",[0 0 0],"LineWidth",0.5);

N = 10;
theta(10) = 0.05;
for ll = 1:length(Ls)
    f0 = rand()*theta(1)/2;
    x0 = [f0 f0 Ts(ll) Ls(ll)];
    theta(11) = lambda(ll);
    for nn = 1:N
        x = x0;
        for tt = 1:1000
            x = updateExoDaisyWorld(x,dt,theta);
        end
        plt.XData = [plt.XData x(4)];
        plt.YData = [plt.YData x(3)];
        pltB.XData = plt.XData;
        pltW.XData = plt.XData;
        pltB.YData = [pltB.YData x(1)];
        pltW.YData = [pltW.YData x(2)];
    end
    drawnow;pause(eps);
end
for ii = 2:length(Ls)
    idx = (plt.XData>Ls(ii-1))&(plt.XData<=Ls(ii));
    pltD.XData = [pltD.XData 0.5*(Ls(ii)+Ls(ii-1))];
    pltD.YData = [pltD.YData max(pltB.YData(idx)+pltW.YData(idx))];
end
xlabel(ax,"Luminosity, $L/L_{opt}$","Interpreter","latex","FontSize",20)
ylabel(ax,"Temperature, $T/T_{opt}$","Interpreter","latex","FontSize",20)
ylabel(ax2,"Daisy Fractions, $f_\alpha$","Interpreter","latex","FontSize",20)


%% Test Environment

close all
fig = figure("Position",[0 0 600 460]);
ax = axes(fig,"Position",[0.13 0.15 0.8 0.8]);

lambda = linspace(-0.7,0.9,100);
Ls = (1+lambda);
Ts = Ls.^0.25;

plot(ax,Ls,Ts,"--","Color",[0 0 0],"Linewidth", 2)
hold(ax,"on");box(ax,"on");ax.LineWidth = 1.0;ax.FontSize = 16;
set(ax,"XLim",[min(Ls),max(Ls)])
limT = [1.31561, 1.23184, 1.17361, 1.12634, 1.08509, 1.0475, 1.01223, ...
    0.978382, 0.9453, 0.912444, 0.879317, 0.845412, 0.810148, 0.772802, ...
    0.732366, 0.6873, 0.63493, 0.569731, 0.475899]/theta(9)^0.25;
for nn = 1:length(limT)
    fill(ax,[min(Ls) max(Ls) max(Ls) min(Ls)],1+[-1 -1 1 1]*limT(nn),...
        [0 0.4 0],"EdgeColor",[0 0 0],"FaceColor",[0 0.4 0],...
        "FaceAlpha",0.015,"LineWidth",0.1,"EdgeAlpha",0.3);
end

plt = plot(ax,Ls(1),Ts(1),".","Color",[0 0 0],"MarkerSize",3);
N = 20;
theta(10) = 0.02;
for ll = 1:length(Ls)
    x0 = [0 0 Ts(ll) Ls(ll)];
    theta(11) = lambda(ll);
    for nn = 1:N
        x = x0;
        for tt = 1:300
            x = updateExoDaisyWorld(x,dt,theta);
        end
        plt.XData = [plt.XData x(4)];
        plt.YData = [plt.YData x(3)];
    end
    drawnow;pause(eps);
end
xlabel(ax,"Luminosity, $L/L_{opt}$","Interpreter","latex","FontSize",20)
ylabel(ax,"Temperature, $T/T_{opt}$","Interpreter","latex","FontSize",20)




