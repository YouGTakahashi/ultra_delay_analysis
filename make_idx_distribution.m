close all
clear all

%PARAMETERS;
%For reward analysis
align = 25; %all_trials(:,25) = reward delivery;
pre = 0.1; %start XX sec after align;
post = 0.5; %end XX sec after align;

base = 26; %all_trials(:,1) = light onset; 24 = well entry
duration = 1.0; %duration for baseline;
base_pre = 0.5;

normalize = 0; %normalized by baseline subtraction;

%for PPE and NPE;
First_T = 5;
Last_T = 5;

%for long omission;
First_T2 = 1;

%LOAD FILES
for cell_count = 1:120; 								
								
if cell_count ==	1	load	DA_unit1	;	cell_1=	unit	; end	%
if cell_count ==	2	load	DA_unit2	;	cell_1=	unit	; end	%
if cell_count ==	3	load	DA_unit3	;	cell_1=	unit	; end	%
if cell_count ==	4	load	DA_unit4	;	cell_1=	unit	; end	%
if cell_count ==	5	load	DA_unit5	;	cell_1=	unit	; end	%
if cell_count ==	6	load	DA_unit6	;	cell_1=	unit	; end	%
if cell_count ==	7	load	DA_unit7	;	cell_1=	unit	; end	%
if cell_count ==	8	load	DA_unit8	;	cell_1=	unit	; end	%
if cell_count ==	9	load	DA_unit9	;	cell_1=	unit	; end	%
if cell_count ==	10	load	DA_unit10	;	cell_1=	unit	; end	%
if cell_count ==	11	load	DA_unit11	;	cell_1=	unit	; end	%
if cell_count ==	12	load	DA_unit12	;	cell_1=	unit	; end	%
if cell_count ==	13	load	DA_unit13	;	cell_1=	unit	; end	%
if cell_count ==	14	load	DA_unit14	;	cell_1=	unit	; end	%
if cell_count ==	15	load	DA_unit15	;	cell_1=	unit	; end	%
if cell_count ==	16	load	DA_unit16	;	cell_1=	unit	; end	%
if cell_count ==	17	load	DA_unit17	;	cell_1=	unit	; end	%
if cell_count ==	18	load	DA_unit18	;	cell_1=	unit	; end	%
if cell_count ==	19	load	DA_unit19	;	cell_1=	unit	; end	%
if cell_count ==	20	load	DA_unit20	;	cell_1=	unit	; end	%
if cell_count ==	21	load	DA_unit21	;	cell_1=	unit	; end	%
if cell_count ==	22	load	DA_unit22	;	cell_1=	unit	; end	%
if cell_count ==	23	load	DA_unit23	;	cell_1=	unit	; end	%
if cell_count ==	24	load	DA_unit24	;	cell_1=	unit	; end	%
if cell_count ==	25	load	DA_unit25	;	cell_1=	unit	; end	%
if cell_count ==	26	load	DA_unit26	;	cell_1=	unit	; end	%
if cell_count ==	27	load	DA_unit27	;	cell_1=	unit	; end	%
if cell_count ==	28	load	DA_unit28	;	cell_1=	unit	; end	%
if cell_count ==	29	load	DA_unit29	;	cell_1=	unit	; end	%
if cell_count ==	30	load	DA_unit30	;	cell_1=	unit	; end	%
if cell_count ==	31	load	DA_unit31	;	cell_1=	unit	; end	%
if cell_count ==	32	load	DA_unit32	;	cell_1=	unit	; end	%
if cell_count ==	33	load	DA_unit33	;	cell_1=	unit	; end	%
if cell_count ==	34	load	DA_unit34	;	cell_1=	unit	; end	%
if cell_count ==	35	load	DA_unit35	;	cell_1=	unit	; end	%
if cell_count ==	36	load	DA_unit36	;	cell_1=	unit	; end	%
if cell_count ==	37	load	DA_unit37	;	cell_1=	unit	; end	%
if cell_count ==	38	load	DA_unit38	;	cell_1=	unit	; end	%
if cell_count ==	39	load	DA_unit39	;	cell_1=	unit	; end	%
if cell_count ==	40	load	DA_unit40	;	cell_1=	unit	; end	%
if cell_count ==	41	load	DA_unit41	;	cell_1=	unit	; end	%
if cell_count ==	42	load	DA_unit42	;	cell_1=	unit	; end	%
if cell_count ==	43	load	DA_unit43	;	cell_1=	unit	; end	%
if cell_count ==	44	load	DA_unit44	;	cell_1=	unit	; end	%
if cell_count ==	45	load	DA_unit45	;	cell_1=	unit	; end	%
if cell_count ==	46	load	DA_unit46	;	cell_1=	unit	; end	%
if cell_count ==	47	load	DA_unit47	;	cell_1=	unit	; end	%
if cell_count ==	48	load	DA_unit48	;	cell_1=	unit	; end	%
if cell_count ==	49	load	DA_unit49	;	cell_1=	unit	; end	%
if cell_count ==	50	load	DA_unit50	;	cell_1=	unit	; end	%
if cell_count ==	51	load	DA_unit51	;	cell_1=	unit	; end	%
if cell_count ==	52	load	DA_unit52	;	cell_1=	unit	; end	%
if cell_count ==	53	load	DA_unit53	;	cell_1=	unit	; end	%
if cell_count ==	54	load	DA_unit54	;	cell_1=	unit	; end	%
if cell_count ==	55	load	DA_unit55	;	cell_1=	unit	; end	%
if cell_count ==	56	load	DA_unit56	;	cell_1=	unit	; end	%
if cell_count ==	57	load	DA_unit57	;	cell_1=	unit	; end	%
if cell_count ==	58	load	DA_unit58	;	cell_1=	unit	; end	%
if cell_count ==	59	load	DA_unit59	;	cell_1=	unit	; end	%
if cell_count ==	60	load	DA_unit60	;	cell_1=	unit	; end	%
if cell_count ==	61	load	DA_unit61	;	cell_1=	unit	; end	%
if cell_count ==	62	load	DA_unit62	;	cell_1=	unit	; end	%
if cell_count ==	63	load	DA_unit63	;	cell_1=	unit	; end	%
if cell_count ==	64	load	DA_unit64	;	cell_1=	unit	; end	%
if cell_count ==	65	load	DA_unit65	;	cell_1=	unit	; end	%
if cell_count ==	66	load	DA_unit66	;	cell_1=	unit	; end	%
if cell_count ==	67	load	DA_unit67	;	cell_1=	unit	; end	%
if cell_count ==	68	load	DA_unit68	;	cell_1=	unit	; end	%
if cell_count ==	69	load	DA_unit69	;	cell_1=	unit	; end	%
if cell_count ==	70	load	DA_unit70	;	cell_1=	unit	; end	%
if cell_count ==	71	load	DA_unit71	;	cell_1=	unit	; end	%
if cell_count ==	72	load	DA_unit72	;	cell_1=	unit	; end	%
if cell_count ==	73	load	DA_unit73	;	cell_1=	unit	; end	%
if cell_count ==	74	load	DA_unit74	;	cell_1=	unit	; end	%
if cell_count ==	75	load	DA_unit75	;	cell_1=	unit	; end	%
if cell_count ==	76	load	DA_unit76	;	cell_1=	unit	; end	%
if cell_count ==	77	load	DA_unit77	;	cell_1=	unit	; end	%
if cell_count ==	78	load	DA_unit78	;	cell_1=	unit	; end	%
if cell_count ==	79	load	DA_unit79	;	cell_1=	unit	; end	%
if cell_count ==	80	load	DA_unit80	;	cell_1=	unit	; end	%
if cell_count ==	81	load	DA_unit81	;	cell_1=	unit	; end	%
if cell_count ==	82	load	DA_unit82	;	cell_1=	unit	; end	%
if cell_count ==	83	load	DA_unit83	;	cell_1=	unit	; end	%
if cell_count ==	84	load	DA_unit84	;	cell_1=	unit	; end	%
if cell_count ==	85	load	DA_unit85	;	cell_1=	unit	; end	%
if cell_count ==	86	load	DA_unit86	;	cell_1=	unit	; end	%
if cell_count ==	87	load	DA_unit87	;	cell_1=	unit	; end	%
if cell_count ==	88	load	DA_unit88	;	cell_1=	unit	; end	%
if cell_count ==	89	load	DA_unit89	;	cell_1=	unit	; end	%
if cell_count ==	90	load	DA_unit90	;	cell_1=	unit	; end	%
if cell_count ==	91	load	DA_unit91	;	cell_1=	unit	; end	%
if cell_count ==	92	load	DA_unit92	;	cell_1=	unit	; end	%
if cell_count ==	93	load	DA_unit93	;	cell_1=	unit	; end	%
if cell_count ==	94	load	DA_unit94	;	cell_1=	unit	; end	%
if cell_count ==	95	load	DA_unit95	;	cell_1=	unit	; end	%
if cell_count ==	96	load	DA_unit96	;	cell_1=	unit	; end	%
if cell_count ==	97	load	DA_unit97	;	cell_1=	unit	; end	%
if cell_count ==	98	load	DA_unit98	;	cell_1=	unit	; end	%
if cell_count ==	99	load	DA_unit99	;	cell_1=	unit	; end	%
if cell_count ==	100	load	DA_unit100	;	cell_1=	unit	; end	%
if cell_count ==	101	load	DA_unit101	;	cell_1=	unit	; end	%
if cell_count ==	102	load	DA_unit102	;	cell_1=	unit	; end	%
if cell_count ==	103	load	DA_unit103	;	cell_1=	unit	; end	%
if cell_count ==	104	load	DA_unit104	;	cell_1=	unit	; end	%
if cell_count ==	105	load	DA_unit105	;	cell_1=	unit	; end	%
if cell_count ==	106	load	DA_unit106	;	cell_1=	unit	; end	%
if cell_count ==	107	load	DA_unit107	;	cell_1=	unit	; end	%
if cell_count ==	108	load	DA_unit108	;	cell_1=	unit	; end	%
if cell_count ==	109	load	DA_unit109	;	cell_1=	unit	; end	%
if cell_count ==	110	load	DA_unit110	;	cell_1=	unit	; end	%
if cell_count ==	111	load	DA_unit111	;	cell_1=	unit	; end	%
if cell_count ==	112	load	DA_unit112	;	cell_1=	unit	; end	%
if cell_count ==	113	load	DA_unit113	;	cell_1=	unit	; end	%
if cell_count ==	114	load	DA_unit114	;	cell_1=	unit	; end	%
if cell_count ==	115	load	DA_unit115	;	cell_1=	unit	; end	%
if cell_count ==	116	load	DA_unit116	;	cell_1=	unit	; end	%
if cell_count ==	117	load	DA_unit117	;	cell_1=	unit	; end	%
if cell_count ==	118	load	DA_unit118	;	cell_1=	unit	; end	%
if cell_count ==	119	load	DA_unit119	;	cell_1=	unit	; end	%
if cell_count ==	120	load	DA_unit120	;	cell_1=	unit	; end	%

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

%CALCULATE FIRING RATE DURING REWARD DELIVERY, OMISSION AND BASELINE;
%FR_1st 
%column 1 = reward response on 1st drop, 
%column 2 = short omission;
%column 3 = long omission;
%column 4 = ultra delay reward;
%column 5 = baseline;
for a = 1:length(all_trials_valid),
   if all_trials_valid(a,25)>0;
       FR_1st(a,1) = (length(find(cell(:,1) <=(all_trials_valid(a,align)+post) & cell(:,1) >=(all_trials_valid(a,align)+pre)))...
           /((all_trials_valid(a,align)+post)-(all_trials_valid(a,align)+pre)));
       FR_1st(a,2) = (length(find(cell(:,1) <=(all_trials_valid(a,24)+0.5+post) & cell(:,1) >=(all_trials_valid(a,24)+0.5+pre)))...
           /((all_trials_valid(a,24)+0.5+post)-(all_trials_valid(a,24)+0.5+pre)));
       FR_1st(a,3) = (length(find(cell(:,1) <=(all_trials_valid(a,24)+3.0+post) & cell(:,1) >=(all_trials_valid(a,24)+3.0+pre)))...
           /((all_trials_valid(a,24)+3.0+post)-(all_trials_valid(a,24)+3.0+pre))); %time of long reward delivery or omission;
       FR_1st(a,4) = (length(find(cell(:,1) <=(all_trials_valid(a,26)+post) & cell(:,1) >=(all_trials_valid(a,26)+pre)))...
           /((all_trials_valid(a,26)+post)-(all_trials_valid(a,26)+pre)));
       FR_1st(a,5) = length(find(cell(:,1) <=(all_trials_valid(a,base)+base_pre+duration) & cell(:,1) >=(all_trials_valid(a,base)+base_pre)))/duration;
   else
       FR_1st(a,1:5) = -999;
   end
end

%Normalization;
if normalize == 1;
    FR_1st_n(:,1) = FR_1st(:,1) - FR_1st(:,5);
    FR_1st_n(:,2) = FR_1st(:,2) - FR_1st(:,5);
    FR_1st_n(:,3) = FR_1st(:,3) - FR_1st(:,5);
    FR_1st_n(:,4) = FR_1st(:,4) - FR_1st(:,5);
end

if normalize == 0;
    FR_1st_n = FR_1st;
end

%EXTRACT EACH CONDITION;
bk1_sh = find(all_trials_valid(:,22) == 1 & all_trials_valid(:,26) > 0);
bk1_lo = find(all_trials_valid(:,22) == 2 & all_trials_valid(:,26) > 0);
bk2_sh = find(all_trials_valid(:,22) == 3 & all_trials_valid(:,26) > 0);
bk2_lo = find(all_trials_valid(:,22) == 4 & all_trials_valid(:,26) > 0);
bk3_sh = find(all_trials_valid(:,22) == 5 & all_trials_valid(:,26) > 0);
bk3_lo = find(all_trials_valid(:,22) == 6 & all_trials_valid(:,26) > 0);
bk4_sh = find(all_trials_valid(:,22) == 7 & all_trials_valid(:,26) > 0);
bk4_lo = find(all_trials_valid(:,22) == 8 & all_trials_valid(:,26) > 0);
bk5_sh = find(all_trials_valid(:,22) == 9 & all_trials_valid(:,26) > 0);
bk5_lo = find(all_trials_valid(:,22) == 10 & all_trials_valid(:,26) > 0);

%EXTRACT FIRING RATE IN EACH CONDITION;
FR_bk1_sh = FR_1st_n(bk1_sh,:);
FR_bk1_lo = FR_1st_n(bk1_lo,:);
FR_bk2_sh = FR_1st_n(bk2_sh,:);
FR_bk2_lo = FR_1st_n(bk2_lo,:);
FR_bk3_sh = FR_1st_n(bk3_sh,:);
FR_bk3_lo = FR_1st_n(bk3_lo,:);
FR_bk4_sh = FR_1st_n(bk4_sh,:);
FR_bk4_lo = FR_1st_n(bk4_lo,:);
FR_bk5_sh = FR_1st_n(bk5_sh,:);
FR_bk5_lo = FR_1st_n(bk5_lo,:);

%EXTRACE FIRST AND LAST SOME TRIALS, DEFINED BY First_T and Last_T;
%FIRST TRIALS;
FR_bk1_sh_F = FR_bk1_sh(1:First_T,:);
FR_bk1_lo_F = FR_bk1_lo(1:First_T,:);
FR_bk2_sh_F = FR_bk2_sh(1:First_T,:);
FR_bk2_lo_F = FR_bk2_lo(1:First_T,:);
FR_bk3_sh_F = FR_bk3_sh(1:First_T,:);
FR_bk3_lo_F = FR_bk3_lo(1:First_T,:);
FR_bk4_sh_F = FR_bk4_sh(1:First_T,:);
FR_bk4_lo_F = FR_bk4_lo(1:First_T,:);
FR_bk5_sh_F = FR_bk5_sh(1:First_T,:);
FR_bk5_lo_F = FR_bk5_lo(1:First_T,:);

%LAST TRIALS;
FR_bk1_sh_L = FR_bk1_sh(length(bk1_sh)-Last_T +1:length(bk1_sh),:);
FR_bk1_lo_L = FR_bk1_lo(length(bk1_lo)-Last_T +1:length(bk1_lo),:);
FR_bk2_sh_L = FR_bk2_sh(length(bk2_sh)-Last_T +1:length(bk2_sh),:);
FR_bk2_lo_L = FR_bk2_lo(length(bk2_lo)-Last_T +1:length(bk2_lo),:);
FR_bk3_sh_L = FR_bk3_sh(length(bk3_sh)-Last_T +1:length(bk3_sh),:);
FR_bk3_lo_L = FR_bk3_lo(length(bk3_lo)-Last_T +1:length(bk3_lo),:);
FR_bk4_sh_L = FR_bk4_sh(length(bk4_sh)-Last_T +1:length(bk4_sh),:);
FR_bk4_lo_L = FR_bk4_lo(length(bk4_lo)-Last_T +1:length(bk4_lo),:);
FR_bk5_sh_L = FR_bk5_sh(length(bk5_sh)-Last_T +1:length(bk5_sh),:);
FR_bk5_lo_L = FR_bk5_lo(length(bk5_lo)-Last_T +1:length(bk5_lo),:);

%COMBINE ALL CELLS;
cat_bk1_sh(cell_count,:) = cat(1,cat(2,mean(FR_bk1_sh_F(:,1)),mean(FR_bk1_sh_L(:,1)),...
    mean(FR_bk1_sh_F(:,2)),mean(FR_bk1_sh_L(:,2)),mean(FR_bk1_sh_F(:,3)),mean(FR_bk1_sh_L(:,3)),...
    mean(FR_bk1_sh_F(:,4)),mean(FR_bk1_sh_L(:,4))));

cat_bk2_sh(cell_count,:) = cat(1,cat(2,mean(FR_bk2_sh_F(:,1)),mean(FR_bk2_sh_L(:,1)),...
    mean(FR_bk2_sh_F(:,2)),mean(FR_bk2_sh_L(:,2)),mean(FR_bk2_sh_F(1:First_T2,3)),mean(FR_bk2_sh_L(1:First_T2,3)),...
    mean(FR_bk2_sh_F(:,4)),mean(FR_bk2_sh_L(:,4))));

cat_bk3_sh(cell_count,:) = cat(1,cat(2,mean(FR_bk3_sh_F(:,1)),mean(FR_bk3_sh_L(:,1)),...
    mean(FR_bk3_sh_F(:,2)),mean(FR_bk3_sh_L(:,2)),mean(FR_bk3_sh_F(1:First_T2,3)),mean(FR_bk3_sh_L(1:First_T2,3)),...
    mean(FR_bk3_sh_F(:,4)),mean(FR_bk3_sh_L(:,4))));

cat_bk4_sh(cell_count,:) = cat(1,cat(2,mean(FR_bk4_sh_F(:,1)),mean(FR_bk4_sh_L(:,1)),...
    mean(FR_bk4_sh_F(:,2)),mean(FR_bk4_sh_L(:,2)),mean(FR_bk4_sh_F(1:First_T2,3)),mean(FR_bk4_sh_L(1:First_T2,3)),...
    mean(FR_bk4_sh_F(:,4)),mean(FR_bk4_sh_L(:,4))));

cat_bk5_sh(cell_count,:) = cat(1,cat(2,mean(FR_bk5_sh_F(:,1)),mean(FR_bk5_sh_L(:,1)),...
    mean(FR_bk5_sh_F(:,2)),mean(FR_bk5_sh_L(:,2)),mean(FR_bk5_sh_F(1:First_T2,3)),mean(FR_bk5_sh_L(1:First_T2,3)),...
    mean(FR_bk5_sh_F(:,4)),mean(FR_bk5_sh_L(:,4))));

cat_bk1_lo(cell_count,:) = cat(1,cat(2,mean(FR_bk1_lo_F(:,1)),mean(FR_bk1_lo_L(:,1)),...
    mean(FR_bk1_lo_F(:,2)),mean(FR_bk1_lo_L(:,2)),mean(FR_bk1_lo_F(:,3)),mean(FR_bk1_lo_L(:,3)),...
    mean(FR_bk1_lo_F(:,4)),mean(FR_bk1_lo_L(:,4))));

cat_bk2_lo(cell_count,:) = cat(1,cat(2,mean(FR_bk2_lo_F(:,1)),mean(FR_bk2_lo_L(:,1)),...
    mean(FR_bk2_lo_F(:,2)),mean(FR_bk2_lo_L(:,2)),mean(FR_bk2_lo_F(:,3)),mean(FR_bk2_lo_L(:,3)),...
    mean(FR_bk2_lo_F(:,4)),mean(FR_bk2_lo_L(:,4))));

cat_bk3_lo(cell_count,:) = cat(1,cat(2,mean(FR_bk3_lo_F(:,1)),mean(FR_bk3_lo_L(:,1)),...
    mean(FR_bk3_lo_F(:,2)),mean(FR_bk3_lo_L(:,2)),mean(FR_bk3_lo_F(:,3)),mean(FR_bk3_lo_L(:,3)),...
    mean(FR_bk3_lo_F(:,4)),mean(FR_bk3_lo_L(:,4))));

cat_bk4_lo(cell_count,:) = cat(1,cat(2,mean(FR_bk4_lo_F(:,1)),mean(FR_bk4_lo_L(:,1)),...
    mean(FR_bk4_lo_F(:,2)),mean(FR_bk4_lo_L(:,2)),mean(FR_bk4_lo_F(:,3)),mean(FR_bk4_lo_L(:,3)),...
    mean(FR_bk4_lo_F(:,4)),mean(FR_bk4_lo_L(:,4))));

cat_bk5_lo(cell_count,:) = cat(1,cat(2,mean(FR_bk5_lo_F(:,1)),mean(FR_bk5_lo_L(:,1)),...
    mean(FR_bk5_lo_F(:,2)),mean(FR_bk5_lo_L(:,2)),mean(FR_bk5_lo_F(:,3)),mean(FR_bk5_lo_L(:,3)),...
    mean(FR_bk5_lo_F(:,4)),mean(FR_bk5_lo_L(:,4))));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Combine blocks;
cat_cis_sh = (cat_bk3_sh + cat_bk4_sh)/2;
cat_cis_lo = (cat_bk3_lo + cat_bk4_lo)/2;

cat_tra_sh = (cat_bk2_sh + cat_bk5_sh)/2;
cat_tra_lo = (cat_bk2_lo + cat_bk5_lo)/2;

cat_sh = (cat_cis_sh + cat_tra_sh)/2;
cat_lo = (cat_cis_lo + cat_tra_lo)/2;

%MAKE DISTRIBUTION PLOTS AND CALCULATE INDEX;
idx_cis_PPE = cat_cis_sh(:,1) - cat_cis_sh(:,2);
idx_tra_PPE = cat_tra_sh(:,1) - cat_tra_sh(:,2);
idx_cis_NPE = cat_cis_lo(:,3) - cat_cis_lo(:,4);
idx_tra_NPE = cat_tra_lo(:,3) - cat_tra_lo(:,4);
idx_cis_OmitL = cat_cis_sh(:,5) - cat_cis_sh(:,6);
idx_tra_OmitL = cat_tra_sh(:,5) - cat_tra_sh(:,6);

idx_PPE = cat_sh(:,1) - cat_sh(:,2);
idx_NPE = cat_lo(:,3) - cat_lo(:,4);

cat_PE_array = cat(2,idx_cis_PPE,idx_cis_NPE,idx_cis_OmitL,idx_tra_PPE,idx_tra_NPE,idx_tra_OmitL);

cat_idx(:,1) = idx_cis_PPE;
cat_idx(:,2) = idx_cis_NPE;
cat_idx(:,3) = idx_cis_OmitL;
cat_idx(:,4) = idx_tra_PPE;
cat_idx(:,5) = idx_tra_NPE;
cat_idx(:,6) = idx_tra_OmitL;

x = -10:1:10;
idx_cis_PPE_score = histc(idx_cis_PPE,x);
idx_tra_PPE_score = histc(idx_tra_PPE,x);
idx_cis_NPE_score = histc(idx_cis_NPE,x);
idx_tra_NPE_score = histc(idx_tra_NPE,x);
idx_cis_OmitL_score = histc(idx_cis_OmitL,x);
idx_tra_OmitL_score = histc(idx_tra_OmitL,x);

idx_PPE_score = histc(idx_PPE,x);
idx_NPE_score = histc(idx_NPE,x);


x2 = x+0.5;
%MAKE FIGURES;
figure %figure1 cisPPE;
bar(x2,idx_cis_PPE_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig1 cis PPE');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure2 transPPE;
bar(x2,idx_tra_PPE_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig2 trans PPE');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure3 cisNPE;
bar(x2,idx_cis_NPE_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig3 cis NPE');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure4 transNPE;
bar(x2,idx_tra_NPE_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig4 trans NPE');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure5 cis Long Omission;
bar(x2,idx_cis_OmitL_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig5 cis LONG OMISSION');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure6 trans Long Omission;
bar(x2,idx_tra_OmitL_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig6 trans LONG OMISSION');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure7 PPE combined;
bar(x2,idx_PPE_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig7 PPE combined');
%axis square;
xlim([-10 10]);

%MAKE FIGURES;
figure %figure8 NPE combined;
bar(x2,idx_NPE_score(:,1),'BarWidth',1,'FaceColor',[0.5 0.5 0.5]);hold on;
bar(0,60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',0.1); hold on
box off;
title('Fig8 NPE combined');
%axis square;
xlim([-10 10]);