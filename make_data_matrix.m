%NAME:  ultra_delay_base;
%Calculate baseline firing and reward firing, then compare by t-test;

close all
clear all

%PARAMETERS;
%For reward analysis
align = 25; %all_trials(:,25) = reward delivery;
pre = 0.1; %start XX sec after align;
post = 0.5; %end XX sec after align;

base = 1; %all_trials(:,1) = light onset; 24 = well entry
duration = 0.4; %duration for baseline;
base_pre = -0.5;

%LOAD FILES
for cell_count = 1:1; 											
											
if cell_count ==	1	load	DA_unit1	;	cell_1=	unit	; end	%	


%CLEAR VARIABLES;
clear all_trials all_trials_valid;
clear lo_left_A lo_left_B lo_right_A lo_right_B sh_left_A sh_left_B sh_right_A sh_right_B;
clear FR_1st

if length(light_on) > length(light_off);  
    light_off(length(light_off)+1) = light_on(length(light_on))+30;
end


%USE Lights on (LO) and Light off (LF) for start and stop of each trilas.
cell=cell_1;
all_trials=cat(2,light_on,light_off); 

%THIS GETS RID OF ZEROS IN MATRIX TO AVOID CONFUSION 
uall_trials(:,3:30) = zeros(1:length(all_trials),1);
all_trials(:,3:30) = -999;

%COLUMN 3 =FIND TIME OF WHEN BROKE BEAM IN ODOR PORT; 
b=1;
for a = 1:length(all_trials),
    if b > length(odor_port_in(:,1))
        break
    end
    if odor_port_in(b,1) < all_trials(a,2) &  odor_port_in(b,1) > all_trials(a,1)
        all_trials(a,3) = odor_port_in(b,1);
        b=b+1;
    else
        all_trials(a,3) = -999;
    end
end

%Column 4 = odor 2 (Odor for Left well);
b=1;
for a = 1:length(all_trials),
    if b > length(odor_left(:,1))
        break
    end
    if odor_left(b,1) < all_trials(a,2) &  odor_left(b,1) > all_trials(a,1)
        all_trials(a,4) =  odor_left(b,1);
        b=b+1;
    else
        all_trials(a,4) = -999;
    end
end

%Column 5 = odor 12 (Odor for Right well); 
b=1;
for a = 1:length(all_trials),
    if b > length(odor_right(:,1))
        break
    end
    if odor_right(b,1) < all_trials(a,2) &  odor_right(b,1) > all_trials(a,1)
        all_trials(a,5) = odor_right(b,1);
        b=b+1;
    else
        all_trials(a,5) = -999;
    end
end

%Column 6 = odor 13 (Choice Odor); 
b=1;
for a = 1:length(all_trials),
    if b > length(odor_choice(:,1))
        break
    end
    if odor_choice(b,1) < all_trials(a,2) &  odor_choice(b,1) > all_trials(a,1),
        all_trials(a,6) = odor_choice(b,1);
        b=b+1;
    else
        all_trials(a,6) = -999;
    end
end

%COLUMN 7 =FIND TIME OF WHEN LEFT ODOR PORT; 
b=1;
for a = 1:length(all_trials),
    if b > length(odor_port_out(:,1))
        break
    end
    if odor_port_out(b,1) < all_trials(a,2) &  odor_port_out(b,1) > all_trials(a,1)
        all_trials(a,7) = odor_port_out(b,1);
        b=b+1;
    else
        all_trials(a,7) = -999;
    end
end

%Column 8 =FIND TIME OF WHEN BROKE BEAM IN LEFT WELL; 
b=1;
for a = 1:length(all_trials)
    if b > length(well_left_in(:,1))
        break
    end
    if well_left_in(b,1) < all_trials(a,2) &  well_left_in(b,1) > all_trials(a,1)
        all_trials(a,8) = well_left_in(b,1);
        b=b+1;
    else
        all_trials(a,8) = -999;
    end
end

%Column 9 =FIND TIME OF WHEN BROKE BEAM IN RIGHT WELL; 
b=1;
for a = 1:length(all_trials),
    if b > length(well_right_in(:,1))
        break
    end
    if well_right_in(b,1) < all_trials(a,2) &  well_right_in(b,1) > all_trials(a,1)
        all_trials(a,9) = well_right_in(b,1);
        b=b+1;
    else
        all_trials(a,9) = -999;
    end
end

%Column 10 = Reward Delivery on LEFT well, 
%Column 12 = Ultra delay on LEFT well, 
for a=1:length(all_trials),
    for b = 1:length(reward_left);
        if reward_left(b) < all_trials(a,8)+4,
            all_trials(a,10) = reward_left(b);
        else if reward_left(b) > all_trials(a,10) & reward_left(b) < all_trials(a,2) & all_trials(a,10) > -999,
                all_trials(a,12) = reward_left(b);
            end
        end
    end
end
clear a b;

%Column 11 = Reward Delivery on Right well, 
%Column 13 = Ultra delay on RIGHT well, 
for a=1:length(all_trials),
    for b = 1:length(reward_right);
        if reward_right(b) < all_trials(a,9)+4,
            all_trials(a,11) = reward_right(b);
        else if reward_right(b) > all_trials(a,11) & reward_right(b) < all_trials(a,2) & all_trials(a,11) > -999,
                all_trials(a,13) = reward_right(b);
            end
        end
    end
end
clear a b;

%Column 14 = FIND TIME OF WHEN LEFT THE LEFT WELL;
b=1;
for a = 1:length(all_trials),
    if b > length(well_left_out(:,1))
        break
    end
    if well_left_out(b,1) < all_trials(a,2) &  well_left_out(b,1) > all_trials(a,1)
        all_trials(a,14) = well_left_out(b,1);
        b=b+1;
    else
        all_trials(a,14) = -999;
    end
end

%COLUMN 15 =FIND TIME OF WHEN LEFT THE RIGHT WELL; 
b=1;
for a = 1:length(all_trials),
    if b > length(well_right_out(:,1))
        break
    end
    if well_right_out(b,1) < all_trials(a,2) &  well_right_out(b,1) > all_trials(a,1)
        all_trials(a,15) = well_right_out(b,1);
        b=b+1;
    else
        all_trials(a,15) = -999;
    end
end

%COLUMN 16 = LENGTH OF DELAY;
for a = 1:length(all_trials(:,1));
    if all_trials(a,10) > -999;
        all_trials(a,16) = all_trials(a,10) - all_trials(a,8);
    else if all_trials(a,11) > -999;
            all_trials(a,16) = all_trials(a,11) - all_trials(a,9);
        else all_trials(a,16) = -999;
        end
    end
end

%COLUMN 17 = DELIVERY OF FLUID A;
b=1;
for a = 1:length(all_trials),
    if b > length(fluid_A(:,1))
        break
    end
    if fluid_A(b,1) < all_trials(a,2) &  fluid_A(b,1) > all_trials(a,1)
        all_trials(a,17) = fluid_A(b,1);
        b=b+1;
    else
        all_trials(a,17) = -999;
    end
end
        
%COLUMN 18 = DELIVERY OF FLUID B;
b=1;
for a = 1:length(all_trials),
    if b > length(fluid_B(:,1))
        break
    end
    if fluid_B(b,1) < all_trials(a,2) &  fluid_B(b,1) > all_trials(a,1)
        all_trials(a,18) = fluid_B(b,1);
        b=b+1;
    else
        all_trials(a,18) = -999;
    end
end

%COLUMN 19 = FIND LONG DELAY ON LEFT ;
b=1;
for a = 1:length(all_trials),
    if b > length(long_delay_left(:,1))
        break
    end
    if a ~= length(all_trials) & long_delay_left(b,1) <= all_trials(a+1,1) &  long_delay_left(b,1) >= all_trials(a,2)
        all_trials(a,19) = long_delay_left(b,1);
        b=b+1;
        %FOR THE LAST TRIAL THERE WILL BE NO UPCOMING LIGHT ON
        else if  a == length(all_trials)   &  long_delay_left(b,1) > all_trials(a,2)
                all_trials(a,19) = long_delay_left(b,1);
                b=b+1;
            else
                all_trials(a,19) = -999;
            end
    end
end

%Column 20 =FIND LONG DELAY ON RIGHT 
b=1;
for a = 1:length(all_trials),
    if b > length(long_delay_right(:,1))
        break
    end
    if a ~= length(all_trials) & long_delay_right(b,1) <= all_trials(a+1,1) &  long_delay_right(b,1) >= all_trials(a,2)
        all_trials(a,20) = long_delay_right(b,1);
        b=b+1;
        %FOR THE LAST TRIAL THERE WILL BE NO UPCOMING LIGHT ON
        else if  a == length(all_trials)   &  long_delay_right(b,1) > all_trials(a,2)
                all_trials(a,20) = long_delay_right(b,1);
                b=b+1;
            else
                all_trials(a,20) = -999;
            end
    end
end

%REMOVE INVALID TRIALS AND MAKE ALL_TRIALS_VALID;
valid_trials = find(all_trials(:,10) > -999 | all_trials(:,11) > -999);
all_trials_valid = all_trials(valid_trials,:);

clear valid_trials;

%Find EACH CONDITION;
%FIND long_left_A;
lo_left_A = find(all_trials_valid(:,10) > -999 & (all_trials_valid(:,10) - all_trials_valid(:,8) > 0.75) & all_trials_valid(:,17) > -999);

%FIND long_left_B;
lo_left_B = find(all_trials_valid(:,10) > -999 & (all_trials_valid(:,10) - all_trials_valid(:,8) > 0.75) & all_trials_valid(:,18) > -999);

%FIND long_right_A;
lo_right_A = find(all_trials_valid(:,11) > -999 & (all_trials_valid(:,11) - all_trials_valid(:,9) > 0.75) & all_trials_valid(:,17) > -999);

%FIND long_right_B;
lo_right_B = find(all_trials_valid(:,11) > -999 & (all_trials_valid(:,11) - all_trials_valid(:,9) > 0.75) & all_trials_valid(:,18) > -999);

%FIND short_left_A;
sh_left_A = find(all_trials_valid(:,10) > -999 & (all_trials_valid(:,10) - all_trials_valid(:,8) < 0.75) & all_trials_valid(:,17) > -999);

%FIND short_left_B;
sh_left_B = find(all_trials_valid(:,10) > -999 & (all_trials_valid(:,10) - all_trials_valid(:,8) < 0.75) & all_trials_valid(:,18) > -999);

%FIND short_right_A;
sh_right_A = find(all_trials_valid(:,11) > -999 & (all_trials_valid(:,11) - all_trials_valid(:,9) < 0.75) & all_trials_valid(:,17) > -999);

%FIND short_right_B;
sh_right_B = find(all_trials_valid(:,11) > -999 & (all_trials_valid(:,11) - all_trials_valid(:,9) < 0.75) & all_trials_valid(:,18) > -999);

%CASE 1; 1st block started sh_left_B - lo_right_A;
if isempty(sh_right_A) == 1;
    %Find start of block2;
    if lo_left_A(1,1) < sh_right_B(1,1);
        bk2_start = lo_left_A(1,1);
    else bk2_start = sh_right_B(1,1);
    end
    
    %Find start of block3;
    if sh_left_A(1,1) < lo_right_B(1,1);
        bk3_start = sh_left_A(1,1);
    else bk3_start = lo_right_B(1,1);
    end
    
    %Find start of block4;
    if sh_left_A(length(sh_left_A),1) > lo_right_B(length(lo_right_B),1) ;
        bk4_start = sh_left_A(length(sh_left_A),1)+1;
    else bk4_start = lo_right_B(length(lo_right_B),1)+1;
    end
    
    %Find_start of block5;
    if lo_left_A(length(lo_left_A),1) > sh_right_B(length(sh_right_B),1) ;
        bk5_start = lo_left_A(length(lo_left_A),1)+1;
    else bk5_start = sh_right_B(length(sh_right_B),1)+1;
    end       
end

%CASE 2; 1st block started sh_left_A - lo_right_B;
if isempty(sh_right_B) == 1;
    %Find start of block2;
    if lo_left_B(1,1) < sh_right_A(1,1);
        bk2_start = lo_left_B(1,1);
    else bk2_start = sh_right_A(1,1);
    end
    
    %Find start of block3;
    if sh_left_B(1,1) < lo_right_A(1,1);
        bk3_start = sh_left_B(1,1);
    else bk3_start = lo_right_A(1,1);
    end
    
    %Find start of block4;
    if sh_left_B(length(sh_left_B),1) > lo_right_A(length(lo_right_A),1) ;
        bk4_start = sh_left_B(length(sh_left_B),1)+1;
    else bk4_start = lo_right_A(length(lo_right_A),1)+1;
    end
    
    %Find_start of block5;
    if lo_left_B(length(lo_left_B),1) > sh_right_A(length(sh_right_A),1) ;
        bk5_start = lo_left_B(length(lo_left_B),1)+1;
    else bk5_start = sh_right_A(length(sh_right_A),1)+1;
    end       
end    
    
%CASE 3; 1st block started lo_left_A - sh_right_B;
if isempty(sh_left_A) == 1;
    %Find start of block2;
    if sh_left_B(1,1) < lo_right_A(1,1);
        bk2_start = sh_left_B(1,1);
    else bk2_start = lo_right_A(1,1);
    end
    
    %Find start of block3;
    if lo_left_B(1,1) < sh_right_A(1,1);
        bk3_start = lo_left_B(1,1);
    else bk3_start = sh_right_A(1,1);
    end
    
    %Find start of block4;
    if lo_left_B(length(lo_left_B),1) > sh_right_A(length(sh_right_A),1) ;
        bk4_start = lo_left_B(length(lo_left_B),1)+1;
    else bk4_start = sh_right_A(length(sh_right_A),1)+1;
    end
    
    %Find_start of block5;
    if sh_left_B(length(sh_left_B),1) > lo_right_A(length(lo_right_A),1) ;
        bk5_start = sh_left_B(length(sh_left_B),1)+1;
    else bk5_start = lo_right_A(length(lo_right_A),1)+1;
    end       
end

%CASE 4; 1st block started lo_left_B - sh_right_A;
if isempty(sh_left_B) == 1;
    %Find start of block2;
    if sh_left_A(1,1) < lo_right_B(1,1);
        bk2_start = sh_left_A(1,1);
    else bk2_start = lo_right_B(1,1);
    end
    
    %Find start of block3;
    if lo_left_A(1,1) < sh_right_B(1,1);
        bk3_start = lo_left_A(1,1);
    else bk3_start = sh_right_B(1,1);
    end
    
    %Find start of block4;
    if lo_left_A(length(lo_left_A),1) > sh_right_B(length(sh_right_B),1) ;
        bk4_start = lo_left_A(length(lo_left_A),1)+1;
    else bk4_start = sh_right_B(length(sh_right_B),1)+1;
    end
    
    %Find_start of block5;
    if sh_left_A(length(sh_left_A),1) > lo_right_B(length(lo_right_B),1) ;
        bk5_start = sh_left_A(length(sh_left_A),1)+1;
    else bk5_start = lo_right_B(length(lo_right_B),1)+1;
    end       
end

%COLUMN 21 = BLOCK ID;
all_trials_valid(1:bk2_start-1,21) = 1;
all_trials_valid(bk2_start:bk3_start-1,21) = 2;
all_trials_valid(bk3_start:bk4_start-1,21) = 3;
all_trials_valid(bk4_start:bk5_start-1,21) = 4;
all_trials_valid(bk5_start:length(all_trials_valid(:,1)),21) = 5;

%COLUMN 22 = CONDITION ID;
%bk1_sh = 1; bk1_lo = 2;
%bk2_sh = 3; bk2_lo = 4;
%bk3_sh = 5; bk3_lo = 6;
%bk4_sh = 7; bk4_lo = 8;
%bk5_sh = 9; bk5_lo = 10;
for a = 1:length(all_trials_valid(:,1));
    if all_trials_valid(a,21) == 1 & all_trials_valid(a,16) < 0.75;
        all_trials_valid(a,22) = 1;
    else if all_trials_valid(a,21) == 1 & all_trials_valid(a,16) > 0.75;
            all_trials_valid(a,22) = 2;
        else if all_trials_valid(a,21) == 2 & all_trials_valid(a,16) < 0.75;
                all_trials_valid(a,22) = 3;
            else if all_trials_valid(a,21) == 2 & all_trials_valid(a,16) > 0.75;
                    all_trials_valid(a,22) = 4;
                else if all_trials_valid(a,21) == 3 & all_trials_valid(a,16) < 0.75;
                        all_trials_valid(a,22) = 5;
                    else if all_trials_valid(a,21) == 3 & all_trials_valid(a,16) > 0.75;
                            all_trials_valid(a,22) = 6;
                        else if all_trials_valid(a,21) == 4 & all_trials_valid(a,16) < 0.75;
                                all_trials_valid(a,22) = 7;
                            else if all_trials_valid(a,21) == 4 & all_trials_valid(a,16) > 0.75;
                                    all_trials_valid(a,22) = 8;
                                else if all_trials_valid(a,21) == 5 & all_trials_valid(a,16) < 0.75;
                                        all_trials_valid(a,22) = 9;
                                    else if all_trials_valid(a,21) == 5 & all_trials_valid(a,16) > 0.75;
                                            all_trials_valid(a,22) = 10;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%Column 23 =TIME OF ODOR ONSET for all odors
for a = 1:length(all_trials_valid),
    if all_trials_valid(a,4) > -999
        all_trials_valid(a,23) = all_trials_valid(a,4); %LEFT ODOR
    end
    if all_trials_valid(a,5) > -999
        all_trials_valid(a,23) = all_trials_valid(a,5); %RIGHT ODOR
    end
    if all_trials_valid(a,6) > -999
        all_trials_valid(a,23) = all_trials_valid(a,6); %CHOICE ODOR
    end
end

%Column 24 =WELL ENTRY right and left
for a = 1:length(all_trials_valid),
    if all_trials_valid(a,8) > -999 
        all_trials_valid(a,24) = all_trials_valid(a,8);
    end
    if all_trials_valid(a,9) > -999
        all_trials_valid(a,24) = all_trials_valid(a,9);
    end
end

%Column 25 = 1st reward in both wells;
for a = 1:length(all_trials_valid),
    if all_trials_valid(a,10) > -999 
        all_trials_valid(a,25) = all_trials_valid(a,10);
    end
    if all_trials_valid(a,11) > -999
        all_trials_valid(a,25) = all_trials_valid(a,11);
    end
end

%Column 26 = Ultra delayed reward in both wells;
for a = 1:length(all_trials_valid),
    if all_trials_valid(a,12) > -999 
        all_trials_valid(a,26) = all_trials_valid(a,12);
    end
    if all_trials_valid(a,13) > -999
        all_trials_valid(a,26) = all_trials_valid(a,13);
    end
end

end