%% norm
clear variables;
%get data for diffrent delta t

advection_eq = @(u) 2.*u;

u0 = @(u) (1).*(u<=-1/2)  + (0) .* (u>-1/2 & u <+1/2)  +  (1).*(u>1/2); %piecwise defined function
u0_sin= @(u) sin((pi)*u);

% method options
%naive
%Lax_Friedrichs
%Lax_Wendroff

boundry = u0_sin;
t_end =0.1;
a = -1;
b = 1;

[xmin,u0min] = fminbnd(@(x) boundry(x),a,b);
[xmax,u0max] = fminbnd(@(x) -boundry(x),a,b);
u0max = -u0max;
uEx =@(x,t) fminbnd(@(u) (u - boundry(x-u*t)).^2,u0min,u0max);



for i=1:3
    [u,all,distance_t] = finiteVolume(10^i,t_end,boundry,advection_eq,'naive');
    
    x = linspace(-1,1,10^i);
    t = zeros(1);
    
    % creating unequal timegrid
    for k=2:size(distance_t,1)
        t(k) = t(k-1) + distance_t(k);
    end
    
    tmp = zeros(size(x,2),size(t,2));
    
    for k=1:size(x,2)
        for j=1:size(t,2)
            tmp(k,j) = uEx(x(k),t(j));
        end
    end
    residum = 0;
    
%     for k=1:size(all,2)
%         tmp1 = sum(all(:,k)) * (x(2)-x(1));
%         tmp2 = integral(@(x) uEx(x,t(k)),a,b);
%         residum(k) = abs(tmp1-tmp2);
%     end
%     l1(i) = norm(residum,1);
%     
    tmp1 = sum(abs(all(:,end))) * (x(2)-x(1));
    tmp2 = integral(@(x) uEx(x,t(end)),a,b);
    l1_naive(i) = norm(tmp1-tmp2,1);
    
end

for i=1:3
    [u,all,distance_t] = finiteVolume(10^i,t_end,boundry,advection_eq,'Lax_Friedrichs');
    
    x = linspace(-1,1,10^i);
    t = zeros(1);
    
    % creating unequal timegrid
    for k=2:size(distance_t,1)
        t(k) = t(k-1) + distance_t(k);
    end
    
    
    tmp = zeros(size(x,2),size(t,2));
    
    for k=1:size(x,2)
        for j=1:size(t,2)
            tmp(k,j) = uEx(x(k),t(j));
        end
    end
    residum = 0;
     
    tmp1 = sum(abs(all(:,end))) * (x(2)-x(1));
    tmp2 = integral(@(x) uEx(x,t(end)),a,b);
    l1_lf(i) = norm(tmp1-tmp2,1);
    
end

for i=1:3
    [u,all,distance_t] = finiteVolume(10^i,t_end,boundry,advection_eq,'Lax_Wendroff');
    
    x = linspace(-1,1,10^i);
    t = zeros(1);
    
    % creating unequal timegrid
    for k=2:size(distance_t,1)
        t(k) = t(k-1) + distance_t(k);
    end
    
    tmp = zeros(size(x,2),size(t,2));
    
    for k=1:size(x,2)
        for j=1:size(t,2)
            tmp(k,j) = uEx(x(k),t(j));
        end
    end
    residum = 0;
     
    tmp1 = sum(abs(all(:,end))) * (x(2)-x(1));
    tmp2 = integral(@(x) uEx(x,t(end)),a,b);
    l1_lw(i) = norm(tmp1-tmp2,1);
    
end




figure
loglog([0,10,100],l1_naive,[0,10,100],l1_lf,[0,10,100],l1_lw);
legend({'naive','Lax-Friedrich','Lax-Wendroff'},'Location','southwest')
title('l1 error');
xlabel('discretisation steps')
ylabel('error')