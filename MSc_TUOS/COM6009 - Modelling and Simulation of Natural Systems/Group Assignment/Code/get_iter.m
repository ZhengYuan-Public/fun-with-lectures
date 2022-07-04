function [health_cows_local_days,infected_cows_local_days,in_rate_days,t_days] = get_iter(size_farm,health_cows_local,infected_cows_local,R_infected,full_rate_infected,full_rate_heal,Time_record_get_infected,D_infected,D_move,get_quarantine,quarantine_T,in_rate,T,t,T_can_infect)

%assume the size of the farm is 1000*1000

% population of the cows
%Popu=150;

%assume the size of the farm is 1000*1000
%size_farm=[-1000,1000]; %size of farm;
%leng_size=size_farm(1,2);

%number of infected cows in beginings
%num_infected=5;

%we assume the cows will move 10m in 1 hours and we use sin and cos to
%represent the random angle they move.
%D_move=10;

%we assume the cows will get affected when they are closed in 10 meter.
%D_infected=10;

%we assume the rate of getting infected =0.5 when they are closed
%R_infected=0.05;

%We assume that the time it takes when it can be examined is 87 days
%T_has_symptom

%See the result on 30 refer as days 24 refer as hours 
%T

% the time choosed to quarantine
%quarantine_T

leng_size=size_farm(1,2);
Angle=0:(2*pi/1000):2*pi;
    
while t<T
    
    
    for x=sort(1:1:length(health_cows_local),2,"descend")
        for y=sort(1:1:length(infected_cows_local),2,"descend")
            if x==length(health_cows_local)
                break
            end
            dis=sqrt((health_cows_local(x,1)-infected_cows_local(y,1))^2+(health_cows_local(x,2)-infected_cows_local(y,2))^2); %the distance of 2 cows
            if dis<=D_infected && binornd(1,R_infected)==1 &&Time_record_get_infected(y,1)==1&&(Time_record_get_infected(y,1)==1||(t-Time_record_get_infected(y,1))>=T_can_infect)
               %conditions that make cows infected
               infected_cows_local=[infected_cows_local;health_cows_local(x,:)];
               health_cows_local(x,:)=[];
               full_rate_infected=[full_rate_infected;full_rate_heal(x,:)];
               full_rate_heal(x,:)=[];
               Time_record_get_infected=[Time_record_get_infected;t];
            end
        end
    end
    
        
    for i=sort(1:1:length(health_cows_local),2,"descend")
        angle_rand=randperm(1000);
        x_cow=health_cows_local(i,1);
        y_cow=health_cows_local(i,2); % conditions that will make make the cows eat
          if (x_cow>=-leng_size*0.1 && x_cow<=-leng_size*0.1+leng_size*0.25) && (y_cow>=-leng_size*0.1 && y_cow<=-leng_size*0.1+leng_size*0.25)&&full_rate_heal(i,1)<=0.5
            full_rate_heal(i,1)=1;
        elseif (x_cow>=-leng_size*0.75&&x_cow<=-leng_size*0.75+leng_size*0.25)&&(y_cow>=leng_size*0.5 && y_cow<=leng_size*0.5+leng_size*0.25)&&full_rate_heal(i,1)<=0.5
            full_rate_heal(i,1)=1;
        elseif (x_cow>=leng_size*0.5&&x_cow<=leng_size*0.5+leng_size*0.25)&&(y_cow>=leng_size*0.5 && y_cow<=leng_size*0.5+leng_size*0.25)&&full_rate_heal(i,1)<=0.5
            full_rate_heal(i,1)=1;
        elseif (x_cow>=-leng_size*0.75&&x_cow<=-leng_size*0.75+leng_size*0.25)&&(y_cow>=-leng_size*0.75&&y_cow<=-leng_size*0.75+leng_size*0.25)&&full_rate_heal(i,1)<=0.5            
            full_rate_heal(i,1)=1;
        elseif (x_cow>=leng_size*0.5&&x_cow<=leng_size*0.5+leng_size*0.25)&&(y_cow>=-leng_size*0.75&&y_cow<=-leng_size*0.75+leng_size*0.25)&&full_rate_heal(i,1)<=0.5
            full_rate_heal(i,1)=1;
          else
            full_rate_heal(i,1)=full_rate_heal(i,1)-0.1;
          end
        %the movement of each cows which is random
        if binornd(1,4/24)==0 
            %the cows sleep 4 hours a day and the time is random, if the cow sleep it will not move.
        health_cows_local(i,1)=health_cows_local(i,1)+D_move*cos(Angle(1,angle_rand(1,1)));
        health_cows_local(i,2)=health_cows_local(i,2)+D_move*sin(Angle(1,angle_rand(1,1)));
        end
    end
    
    
    
    for j=1:1:length(infected_cows_local)% conditions that will make make the cows eat
        angle_rand=randperm(1000);
        x_cow=infected_cows_local(i,1);
        y_cow=infected_cows_local(i,2);
          if (x_cow>=-leng_size*0.1 && x_cow<=-leng_size*0.1+leng_size*0.25) && (y_cow>=-leng_size*0.1 && y_cow<=-leng_size*0.1+leng_size*0.25)&&full_rate_infected(i,1)<=0.5
            full_rate_infected(i,1)=1;
        elseif (x_cow>=-leng_size*0.75&&x_cow<=-leng_size*0.75+leng_size*0.25)&&(y_cow>=leng_size*0.5 && y_cow<=leng_size*0.5+leng_size*0.25)&&full_rate_infected(i,1)<=0.5
            full_rate_infected(i,1)=1;
        elseif (x_cow>=leng_size*0.5&&x_cow<=leng_size*0.5+leng_size*0.25)&&(y_cow>=leng_size*0.5 && y_cow<=leng_size*0.5+leng_size*0.25)&&full_rate_infected(i,1)<=0.5
            full_rate_infected(i,1)=1;
        elseif (x_cow>=-leng_size*0.75&&x_cow<=-leng_size*0.75+leng_size*0.25)&&(y_cow>=-leng_size*0.75&&y_cow<=-leng_size*0.75+leng_size*0.25)&&full_rate_infected(i,1)<=0.5            
            full_rate_infected(i,1)=1;
        elseif (x_cow>=leng_size*0.5&&x_cow<=leng_size*0.5+leng_size*0.25)&&(y_cow>=-leng_size*0.75&&y_cow<=-leng_size*0.75+leng_size*0.25)&&full_rate_infected(i,1)<=0.5
            full_rate_infected(i,1)=1;
          else
            full_rate_infected(i,1)=full_rate_infected(i,1)-0.1;
          end
        if binornd(1,4/24)==0 
            %the cows sleep 4 hours a day and the time is random, if the cow sleep it will not move.
        infected_cows_local(j,1)=infected_cows_local(j,1)+D_move*cos(Angle(1,angle_rand(1,1)));
        infected_cows_local(j,2)=infected_cows_local(j,2)+D_move*sin(Angle(1,angle_rand(1,1)));
        end
    end
    
    %control methods
    %the condition that will quarantine cows
    x=randperm(360);
    T_has_symptom=x(1,1)*24;
    if t==quarantine_T&& get_quarantine==true % the quarantine process
        for i=sort(1:1:length(Time_record_get_infected),2,"descend")
            if Time_record_get_infected(i,1)+T_has_symptom<quarantine_T
               infected_cows_local(i,:)=[]; 
            end
            if length(infected_cows_local)==1
                break
            end
        end
    end
    
    if rem(t,24)==0 % record the change of infected rate through days
    in_rate(1,t/24+1)=length(infected_cows_local(:,1))/(length(infected_cows_local(:,1))+length(health_cows_local(:,1)));
    end
 
    
    
    t=t+1;




end
health_cows_local_days=health_cows_local;
infected_cows_local_days=infected_cows_local;
in_rate_days=in_rate;
t_days=t;

end
