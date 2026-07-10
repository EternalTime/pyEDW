load("data_040.mat")

Nex = 170;
f1 = squeeze(data(Nex,:,1));
f2 = squeeze(data(Nex,:,2));
L = squeeze(data(Nex,:,4));
T = squeeze(data(Nex,:,3));
A = 1-theta(2)+theta(3)*f1+theta(4)*f2;

minF = 0.9;maxF = 1.1;
%f1min = minF*min(f1);f1max = maxF*max(f1);
df1 =5*std(f1);f1min = mean(f1)-df1;f1max = mean(f1)+df1;
%f2min = minF*min(f2);f2max = maxF*max(f2);
df2 =5*std(f2);f2min = mean(f2)-df2;f2max = mean(f2)+df2;
%Tmin = minF*min(T);Tmax = maxF*max(T);
dT1 =5*std(T);Tmin = mean(T)-dT1;Tmax = mean(T)+dT1;
%Lmin = minF*min(L);Lmax = maxF*max(L);
dL1 =5*std(L);Lmin = mean(L)-dL1;Lmax = mean(L)+dL1;

Nc = 9;
Nb = 41;
close all
fig = figure("Position",[0 0 1400 700]);

ax14 = axes(fig,"Position",[0.1*0.5 0.1 0.21*0.5 0.21],"Box","on")
ax13 = axes(fig,"Position",[0.1*0.5 0.32 0.21*0.5 0.21],"Box","on")
ax12 = axes(fig,"Position",[0.1*0.5 0.54 0.21*0.5 0.21],"Box","on")
ax11 = axes(fig,"Position",[0.1*0.5 0.76 0.21*0.5 0.21],"Box","on")

ax24 = axes(fig,"Position",[0.32*0.5 0.1 0.21*0.5 0.21],"Box","on")
ax23 = axes(fig,"Position",[0.32*0.5 0.32 0.21*0.5 0.21],"Box","on")
ax22 = axes(fig,"Position",[0.32*0.5 0.54 0.21*0.5 0.21],"Box","on")

ax34 = axes(fig,"Position",[0.54*0.5 0.1 0.21*0.5 0.21],"Box","on")
ax33 = axes(fig,"Position",[0.54*0.5 0.32 0.21*0.5 0.21],"Box","on")

ax44 = axes(fig,"Position",[0.76*0.5 0.1 0.21*0.5 0.21],"Box","on")



[y,x] = histcounts(f1,linspace(f1min,f1max,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt11 = fill(ax11,[x fliplr(x)],[y zeros(size(y))-1],RochesterBlue,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",RochesterBlue);
set(ax11,"LineWidth",1.5,"XTickLabel",[],"XLim",[f1min f1max],"YLim",[0,1.1*max(y)],"FontSize",16)


[y,x] = histcounts(f2,linspace(f2min,f2max,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt22 = fill(ax22,[x fliplr(x)],[y -ones(size(y))],ForagerPink,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",ForagerPink);
set(ax22,"LineWidth",1.5,"XTickLabel",[],"YTickLabel",[],"XLim",[f2min f2max],"YLim",[0 1.1*max(y)],"FontSize",16)

[y,x] = histcounts(T,linspace(Tmin,Tmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt33 = fill(ax33,[x fliplr(x)],[y -ones(size(y))],ResourceGreen,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",ResourceGreen);
set(ax33,"LineWidth",1.5,"XTickLabel",[],"YTickLabel",[],"XLim",[Tmin Tmax],"YLim",[0,1.1*max(y)],"FontSize",16)

[y,x] = histcounts(L,linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt44 = fill(ax44,[x fliplr(x)],[y -ones(size(y))],DartmouthGreen,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",DartmouthGreen);
set(ax44,"LineWidth",1.5,"YTickLabel",[],"XLim",[Lmin,Lmax],"YLim",[0,1.1*max(y)],"FontSize",16)

[z,x,y] = histcounts2(f1,f2,linspace(f1min,f1max,Nb),linspace(f2min,f2max,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
kernel = [1 2 1;2 3 2;1 2 1];kernel = kernel/sum(kernel(:));
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(ax12,x,y,z,Nc,"Color",[0 0 0],"LineWidth",1.5);
set(ax12,"LineWidth",1.5,"XTickLabel",[],"XLim",[f1min f1max],"YLim",[f2min f2max],"FontSize",16)

[z,x,y] = histcounts2(f1,T,linspace(f1min,f1max,Nb),linspace(Tmin,Tmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(ax13,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(ax13,"LineWidth",1.5,"XTickLabel",[],"XLim",[f1min f1max],"YLim",[Tmin Tmax],"FontSize",16)

[z,x,y] = histcounts2(f1,L,linspace(f1min,f1max,Nb),linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(ax14,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(ax14,"LineWidth",1.5,"XLim",[f1min f1max],"YLim",[Lmin,Lmax],"FontSize",16)

[z,x,y] = histcounts2(f2,T,linspace(f2min,f2max,Nb),linspace(Tmin,Tmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(ax23,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(ax23,"LineWidth",1.5,"XTickLabel",[],"YTickLabel",[],"XLim",[f2min f2max],"YLim",[Tmin Tmax],"FontSize",16)

[z,x,y] = histcounts2(f2,L,linspace(f2min,f2max,Nb),linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(ax24,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(ax24,"LineWidth",1.5,"YTickLabel",[],"XLim",[f2min f2max],"YLim",[Lmin,Lmax],"FontSize",16)

[z,x,y] = histcounts2(T,L,linspace(Tmin,Tmax,Nb),linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(ax34,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(ax34,"LineWidth",1.5,"YTickLabel",[],"XLim",[Tmin Tmax],"YLim",[Lmin Lmax],"FontSize",16)


xlabel(ax14,"$f_B$","Interpreter","latex","FontSize",20)
ylabel(ax11,"$f_B$","Interpreter","latex","FontSize",20)
xlabel(ax24,"$f_W$","Interpreter","latex","FontSize",20)
ylabel(ax12,"$f_W$","Interpreter","latex","FontSize",20)
xlabel(ax34,"$T$","Interpreter","latex","FontSize",20)
ylabel(ax13,"$T$","Interpreter","latex","FontSize",20)
xlabel(ax44,"$L$","Interpreter","latex","FontSize",20)
ylabel(ax14,"$L$","Interpreter","latex","FontSize",20)


Nd = 80;
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

Nc = 9;
Nb = 81;


axF14 = axes(fig,"Position",[0.5+0.1*0.5 0.1 0.21*0.5 0.21],"Box","on");
axF13 = axes(fig,"Position",[0.5+0.1*0.5 0.32 0.21*0.5 0.21],"Box","on");
axF12 = axes(fig,"Position",[0.5+0.1*0.5 0.54 0.21*0.5 0.21],"Box","on");
axF11 = axes(fig,"Position",[0.5+0.1*0.5 0.76 0.21*0.5 0.21],"Box","on");

axF24 = axes(fig,"Position",[0.5+0.32*0.5 0.1 0.21*0.5 0.21],"Box","on");
axF23 = axes(fig,"Position",[0.5+0.32*0.5 0.32 0.21*0.5 0.21],"Box","on");
axF22 = axes(fig,"Position",[0.5+0.32*0.5 0.54 0.21*0.5 0.21],"Box","on");

axF34 = axes(fig,"Position",[0.5+0.54*0.5 0.1 0.21*0.5 0.21],"Box","on");
axF33 = axes(fig,"Position",[0.5+0.54*0.5 0.32 0.21*0.5 0.21],"Box","on");

axF44 = axes(fig,"Position",[0.5+0.76*0.5 0.1 0.21*0.5 0.21],"Box","on");

hold(axF13,"on");
fill(axF13,[0 1 1 0],[1 1 1 1]+[-1 -1 1 1]*0.4*delTs(Nd),...
    DartmouthGreen, "FaceColor",DartmouthGreen,"FaceAlpha",0.1,"EdgeColor","none")
hold(axF23,"on");
fill(axF23,[0 1 1 0],[1 1 1 1]+[-1 -1 1 1]*0.4*delTs(Nd),...
    DartmouthGreen, "FaceColor",DartmouthGreen,"FaceAlpha",0.1,"EdgeColor","none")
hold(axF34,"on");
fill(axF34,[1 1 1 1]+[-1 -1 1 1]*0.4*delTs(Nd),[0 10 10 0], ...
    DartmouthGreen, "FaceColor",DartmouthGreen,"FaceAlpha",0.1,"EdgeColor","none")

[y,x] = histcounts(f1,linspace(f1min,f1max,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt11 = fill(axF11,[x fliplr(x)],[y zeros(size(y))-1],RochesterBlue,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",RochesterBlue);
set(axF11,"LineWidth",1.5,"XTickLabel",[],"XLim",[f1min f1max],"YLim",[0,1.1*max(y)],"FontSize",16)


[y,x] = histcounts(f2,linspace(f2min,f2max,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt22 = fill(axF22,[x fliplr(x)],[y -ones(size(y))],ForagerPink,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",ForagerPink);
set(axF22,"LineWidth",1.5,"XTickLabel",[],"YTickLabel",[],"XLim",[f2min f2max],"YLim",[0 1.1*max(y)],"FontSize",16)

[y,x] = histcounts(T,linspace(Tmin,Tmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
plt33 = fill(axF33,[x fliplr(x)],[y -ones(size(y))],ResourceGreen,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",ResourceGreen);
set(axF33,"LineWidth",1.5,"XTickLabel",[],"YTickLabel",[],"XLim",[Tmin Tmax],"YLim",[0,1.1*max(y)],"FontSize",16)

[y,x] = histcounts(L,linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[1 2 2 1]/6,"same");
fill(axF44,[x fliplr(x)],[y -ones(size(y))],DartmouthGreen,"FaceAlpha",0.3,"LineWidth",2.0,"EdgeColor",DartmouthGreen);
set(axF44,"LineWidth",1.5,"YTickLabel",[],"XLim",[Lmin,Lmax],"YLim",[0,1.1*max(y)],"FontSize",16)

[z,x,y] = histcounts2(f1,f2,linspace(f1min,f1max,Nb),linspace(f2min,f2max,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
kernel = [1 2 1;2 3 2;1 2 1];kernel = kernel/sum(kernel(:));
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(axF12,x,y,z,Nc,"Color",[0 0 0],"LineWidth",2.0);
set(axF12,"LineWidth",1.5,"XTickLabel",[],"XLim",[f1min f1max],"YLim",[f2min f2max],"FontSize",16)

[z,x,y] = histcounts2(f1,T,linspace(f1min,f1max,Nb),linspace(Tmin,Tmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(axF13,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(axF13,"LineWidth",1.5,"XTickLabel",[],"XLim",[f1min f1max],"YLim",[Tmin Tmax],"FontSize",16)

[z,x,y] = histcounts2(f1,L,linspace(f1min,f1max,Nb),linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(axF14,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(axF14,"LineWidth",1.5,"XLim",[f1min f1max],"YLim",[Lmin,Lmax],"FontSize",16)

[z,x,y] = histcounts2(f2,T,linspace(f2min,f2max,Nb),linspace(Tmin,Tmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(axF23,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(axF23,"LineWidth",1.5,"XTickLabel",[],"YTickLabel",[],"XLim",[f2min f2max],"YLim",[Tmin Tmax],"FontSize",16)

[z,x,y] = histcounts2(f2,L,linspace(f2min,f2max,Nb),linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
contour(axF24,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(axF24,"LineWidth",1.5,"YTickLabel",[],"XLim",[f2min f2max],"YLim",[Lmin,Lmax],"FontSize",16)

[z,x,y] = histcounts2(T,L,linspace(Tmin,Tmax,Nb),linspace(Lmin,Lmax,Nb),"Normalization","Probability");
x = conv(x,[0.5 0.5],"valid");
y = conv(y,[0.5 0.5],"valid");
for ii = 1:2
z = conv2(z,kernel,"same");
end
[~,plt34] = contour(axF34,x,y,z',Nc,"Color",[0 0 0],"LineWidth",1.5);
set(axF34,"LineWidth",1.5,"YTickLabel",[],"XLim",[Tmin Tmax],"YLim",[Lmin Lmax],"FontSize",16)



xlabel(axF14,"$f_B$","Interpreter","latex","FontSize",20)
ylabel(axF11,"$f_B$","Interpreter","latex","FontSize",20)
xlabel(axF24,"$f_W$","Interpreter","latex","FontSize",20)
ylabel(axF12,"$f_W$","Interpreter","latex","FontSize",20)
xlabel(axF34,"$T$","Interpreter","latex","FontSize",20)
ylabel(axF13,"$T$","Interpreter","latex","FontSize",20)
xlabel(axF44,"$L$","Interpreter","latex","FontSize",20)
ylabel(axF14,"$L$","Interpreter","latex","FontSize",20)
