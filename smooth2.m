function S=smooth2(H, smooth2)

%smooth= number of bins to be averaged.
% this gives you colum. 

dt=length(H);
%S = zeros(1,dt);
S = zeros(dt,1);
switch  smooth2
    case 0
        S=H;
    
    case 3	
        S(2:dt-2)=(H(1:dt-3)+H(2:dt-2)+H(3:dt-1))/3;
        %S(1:3)=H(2:3); S(dt-1:dt)=H(dt-1:dt);
        %S(1)=sum(H(2:4))/3; S(2)=sum(H(3:6))/4;
        %S(dt-1)=sum(H(dt-3:dt))/4; S(dt)=sum(H(dt-2:dt-1))/2;
         S(2)=sum((H(3)+S(3:4)))/3; S(1)=sum(S(2:4))/3; 
         S(dt-1)=sum(H(dt-1)+S(dt-2:dt-1))/3; S(dt)=sum(S(dt-2:dt))/3;
       
    case 5
       S(3:dt-2)=(H(1:dt-4)+H(2:dt-3)+H(3:dt-2)+H(4:dt-1)+H(5:dt))/5; 
       S(1)=H(1); S(2)=sum(H(1:3))/3; 
       S(dt-1)=sum(H(dt-2:dt))/3; S(dt)=H(dt);
    case 9	
       S(5:dt-4)=(H(1:dt-8)+H(2:dt-7)+H(3:dt-6)+H(4:dt-5)+H(5:dt-4)+H(6:dt-3)+H(7:dt-2)+H(8:dt-1)+H(9:dt))/9;
       S(1)=H(1); S(2)=sum(H(1:3))/3; S(3)=sum(H(1:5))/5; S(4)=sum(H(1:7))/7; 
       S(dt-3)=sum(H(dt-6:dt))/7; S(dt-2)=sum(H(dt-4:dt))/5; S(dt-1)=sum(H(dt-2:dt))/3; S(dt)=H(dt);
       
          
    case 11	
       S(6:dt-5)= (H(1:dt-10)+H(2:dt-9)+H(3:dt-8)+H(4:dt-7)+H(5:dt-6)+H(6:dt-5)+H(7:dt-4)+H(8:dt-3)+H(9:dt-2)+H(10:dt-1)+H(11:dt))/9;
       S(1)=H(1); S(2)=sum(H(1:3))/3; S(3)=sum(H(1:5))/5; S(4)=sum(H(1:7))/7; S(5)=sum(H(1:9))/9;
       S(dt-4)=sum(H(dt-8:dt))/9; S(dt-3)=sum(H(dt-6:dt))/7; S(dt-2)=sum(H(dt-4:dt))/5; S(dt-1)=sum(H(dt-2:dt))/3; S(dt)=H(dt);
       

       
       
    otherwise 
		 error('****select from 3,5,or 9 ****');
end

