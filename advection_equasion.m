clear variables;
%% settings advection equasion

% boundry conditions
u0 = @(u) (1).*(u<=-1/2)  + (0) .* (u>-1/2 & u <+1/2)  +  (1).*(u>1/2); %piecwise defined function
u0_sin= @(u) sin((pi)*u);

%naive
%Lax_Friedrichs
%Lax_Wendroff

method = 'Lax_Wendroff';
boundry = u0;
N =200;
t_end = 0.2;

%flux function
advection_eq = @(u) 2.*u;

%% approximate advection equasion


[u,all,distance_t] = finiteVolume(N,t_end,boundry,advection_eq,method);


x = linspace(-1,1,N);
t = zeros(1);

% creating unequal timegrid
for i=2:size(distance_t,1)
    t(i) = t(i-1) + distance_t(i);
end

%% exact advection equasion
a = -1;
b = 1;

uEx =@(x,t) boundry(x-2*t);


time = linspace(0,t_end,1000);
axis = linspace(-1,1,3*N);

for i=1:size(axis,2)
    for j=1:size(time,2)
        tmp(i,j) = uEx(axis(i),time(j));
    end
end

%% plot
figure;
sgtitle('advection equasion F(u) =2u');
FaceColor = 'interp';
LineWidth = 7;
colormap('parula');
%imagesc(x,y,all);

%approximate Solution
subplot(1,2,1)
pcolor(x,t,all');
title('approximate solution');

colorbar;

% exact solution
subplot(1,2,2)

plot2 = pcolor(axis,time,tmp');
title('exact solution');
set(plot2, 'EdgeColor', 'none');
colorbar;

