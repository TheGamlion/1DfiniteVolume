function [u,u_total,distance_t] = finiteVolume(N,t_end,u0,flux_function,num_flux)

%allocate memory
delta_x = 2/N;


%fill inital condition
u_prev = zeros(N,1);
for i=1:N
    u_prev(i) = 1/delta_x*integral(u0,-1+(delta_x*(i-1)),-1+(delta_x*i));
end

u_total = zeros(N,1);
u_total(:,1)= u_prev;
i = 2;

distance_t = zeros(2,1);

u = zeros(N,1);

t = 0;
%update u:
while t<t_end
    %stepsize needs to be eval. every step to meet cfl condition
    der_max = max_der(flux_function,u_prev,N,delta_x);
    %disp(der_max)
    delta_t = delta_x/abs(der_max);
    for x=1:N
        %special cases due to boundry conditions
        if(x == 1)
            if(strcmp(num_flux,'naive'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_naive(flux_function,u_prev(end),u_prev(x),delta_x,delta_t)-flux_naive(flux_function,u_prev(x),u_prev(x+1),delta_x,delta_t));
            end
            if(strcmp(num_flux,'Lax_Friedrichs'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_Lax_Friedrichs(flux_function,u_prev(end),u_prev(x),delta_x,delta_t)-flux_Lax_Friedrichs(flux_function,u_prev(x),u_prev(x+1),delta_x,delta_t));
            end
            if(strcmp(num_flux,'Lax_Wendroff'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_Lax_Wendroff(flux_function,u_prev(end),u_prev(x),delta_x,delta_t)-flux_Lax_Wendroff(flux_function,u_prev(x),u_prev(x+1),delta_x,delta_t));
            end
        elseif(x == N)
            if(strcmp(num_flux,'naive'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_naive(flux_function,u_prev(x-1),u_prev(x),delta_x,delta_t)-flux_naive(flux_function,u_prev(x),u_prev(1),delta_x,delta_t));
            end
            if(strcmp(num_flux,'Lax_Friedrichs'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_Lax_Friedrichs(flux_function,u_prev(x-1),u_prev(x),delta_x,delta_t)-flux_Lax_Friedrichs(flux_function,u_prev(x),u_prev(1),delta_x,delta_t));
            end
            if(strcmp(num_flux,'Lax_Wendroff'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_Lax_Wendroff(flux_function,u_prev(x-1),u_prev(x),delta_x,delta_t)-flux_Lax_Wendroff(flux_function,u_prev(x),u_prev(1),delta_x,delta_t));
            end
        else
            %main calculation part
            if(strcmp(num_flux,'naive'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_naive(flux_function,u_prev(x-1),u_prev(x),delta_x,delta_t)-flux_naive(flux_function,u_prev(x),u_prev(x+1),delta_x,delta_t));
            end
            if(strcmp(num_flux,'Lax_Friedrichs'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_Lax_Friedrichs(flux_function,u_prev(x-1),u_prev(x),delta_x,delta_t)-flux_Lax_Friedrichs(flux_function,u_prev(x),u_prev(x+1),delta_x,delta_t));
            end
            if(strcmp(num_flux,'Lax_Wendroff'))
                u(x) = u_prev(x) + delta_t/delta_x*(flux_Lax_Wendroff(flux_function,u_prev(x-1),u_prev(x),delta_x,delta_t)-flux_Lax_Wendroff(flux_function,u_prev(x),u_prev(x+1),delta_x,delta_t));
            end
        end
    end
    u_total(:,i)= u;
    distance_t(i,1) = delta_t;
    i = i+1;
    u_prev = u;
    t = t + delta_t;
    disp(t)
end

end

function [der_max] =max_der(flux_function,u_prev,N,delta_x)

der_max = abs(flux_function(u_prev(1)) -flux_function(u_prev(end)))/delta_x;

for i=2:N
    if(der_max < (flux_function(u_prev(i)) -flux_function(u_prev(i-1)))/delta_x)
        der_max = abs(flux_function(u_prev(i)) -flux_function(u_prev(i-1)))/delta_x;
    end
end
end

function F = flux_naive(flux_function,u_l,u_r,~,~)
F = 1/2*(flux_function(u_l)+flux_function(u_r));
end
function F = flux_Lax_Friedrichs(flux_function,u_l,u_r,delta_x,delta_t)
F = 1/2*(flux_function(u_l)+flux_function(u_r)) - delta_x/(2*delta_t)*(u_r-u_l);
end
function F = flux_Lax_Wendroff(flux_function,u_l,u_r,delta_x,delta_t)
if((u_r-u_l) ~=0)
    F = 1/2*(flux_function(u_l)+flux_function(u_r)) - delta_t/(2*delta_x)*(u_r-u_l)*((flux_function(u_l)-flux_function(u_r))/(u_l-u_r))^2;
else
    F = 1/2*(flux_function(u_l)+flux_function(u_r));
end
end