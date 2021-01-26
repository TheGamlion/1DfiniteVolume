clear variables;
%% settings burgers equasion
u0 = @(u) (1).*(u<=-1/2)  + (0) .* (u>-1/2 & u <+1/2)  +  (1).*(u>1/2); %piecwise defined function
u0_sin= @(u) sin((pi)*u);

% method options
%naive
%Lax_Friedrichs
%Lax_Wendroff

method = 'Lax_Friedrichs';
boundry = u0;
N =60;
t_end =0.5;

burgers_eq = @(u) 1/2.*u.^2;

%% approximate burgers equasion

[u,all,distance_t] = finiteVolume(N,t_end,boundry,burgers_eq,method);


x = linspace(-1,1,N);
t = zeros(1);

% creating unequal timegrid
for i=2:size(distance_t,1)
    t(i) = t(i-1) + distance_t(i);
end

%% exact burgers equasion
a = -1;
b = 1;

[xmin,u0min] = fminbnd(@(x) boundry(x),a,b);
[xmax,u0max] = fminbnd(@(x) -boundry(x),a,b);
u0max = -u0max;
uEx =@(x,t) fminbnd(@(u) (u - boundry(x-u*t)).^2,u0min,u0max);


time = linspace(0,t_end,100);
axis = linspace(-1,1,4*N);

for i=1:size(axis,2)
    for j=1:size(time,2)
        tmp(i,j) = uEx(axis(i),time(j));
    end
end

%% plot
figure;
sgtitle('burgers equasion F(u) =1/2u^u');

LineWidth = 0.01;
colormap('parula');
%imagesc(x,y,all);

%approximate Solution
subplot(1,2,1)
plot1 = pcolor(x,t,all');
plot1.LineWidth = 0.1;
% plot1.FaceColor = 'interp';
% set(plot1, 'EdgeColor', 'none');
title('approximate solution');

colorbar;

%exact solution
subplot(1,2,2)

plot2 = pcolor(axis,time,tmp');
%plot2.FaceColor = 'interp';
title('exact solution');
set(plot2, 'EdgeColor', 'none');
title('exact solution');
colorbar;
