%% This function generates the coordinates of nodes

function [cor_user] = coordinates(K,radius_BS,min_dis_user)

Mokh_user=0;
t=0;
while 6>0 
    Dis_user_BS=min_dis_user+(radius_BS-min_dis_user)*unifrnd(0,1);
    angle_user = 2*pi*unifrnd(0,1);
    crd_user=Dis_user_BS.*( cos(angle_user) + i*sin(angle_user) );
    t=t+1;
    Mokh_user=[Mokh_user,crd_user];
    if (t==K)
        break
    end
end
Mokh_user(1)=[]; %Coordinated of users in the cell
cor_user=Mokh_user;  
end

