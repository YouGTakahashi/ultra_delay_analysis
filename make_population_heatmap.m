close all
clear all

%PARAMETERS;
%For reward analysis
align = 24; %25 = reward delivery; 24 = well entry;
pre = -0.5; %start XX sec after align;
post = 6.5; %end XX sec after align;
bin_size = 0.1;

base = 24; %all_trials(:,1) = light onset; 24 = well entry
duration = 0.4; %duration for baseline;
base_pre = 0.1;

First_T = 10;
Last_T = 10;

normalize = 0; %normalized by baseline subtraction;

filtering = 1;
windowSize = 2;

mag = 7.5;

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
clear FR_base spike_idx spike_idx_n;

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

%CALCULATE FIRING RATE IN EACH BIN;
bin_centers = pre: bin_size: post;
bin_centers = bin_centers - bin_size/2;

for a = 1:length(all_trials_valid);
    spike_times_idx = find(cell(:,1) <= (all_trials_valid(a,align)+post)  & cell(:,1) >= (all_trials_valid(a,align)+pre)); %Find Y-axix of Spikes in Cell during Reward Delivery 
    spike_times = cell(spike_times_idx); %Find Times of Spikes based on Y-axis,
    spike_times_normalized = spike_times - all_trials_valid(a,align); %Normalized by Reward Delivery,
    if filtering == 0
        spike_idx(a,:) = smooth2(hist(spike_times_normalized,bin_centers),0)/bin_size; %Smoothing
    else if filtering == 1
            spike_idx(a,:) = filter(ones(1,windowSize)/windowSize,1,(smooth2(hist(spike_times_normalized,bin_centers),0)/bin_size)); %Smoothing
        end
    end
end

if normalize == 1;
    for a = 1:length(all_trials_valid(:,1));
        FR_base(a,1) = length(find(cell(:,1) <=(all_trials_valid(a,base)+base_pre+duration) & cell(:,1) >=(all_trials_valid(a,base)+base_pre)))/duration;
        spike_idx_n(a,:) = spike_idx(a,:)- FR_base(a,1);
    end
end

if normalize == 0;
    spike_idx_n = spike_idx; %change raw data name for consistency;
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
%block 1;
FR_bk1_sh = spike_idx_n(bk1_sh,:); %firing on short reward;
FR_bk1_lo = spike_idx_n(bk1_lo,:); %firing on long reward;

%block 2;
FR_bk2_sh = spike_idx_n(bk2_sh,:); %firing on short reward; 
FR_bk2_lo = spike_idx_n(bk2_lo,:); %firing on long reward;

%block 3;
FR_bk3_sh = spike_idx_n(bk3_sh,:); %firing on short reward; 
FR_bk3_lo = spike_idx_n(bk3_lo,:); %firing on long reward;

%block 4;
FR_bk4_sh = spike_idx_n(bk4_sh,:); %firing on short reward;
FR_bk4_lo = spike_idx_n(bk4_lo,:); %iring on long reward;

%block 5;
FR_bk5_sh = spike_idx_n(bk5_sh,:); %firing on short reward;
FR_bk5_lo = spike_idx_n(bk5_lo,:); %iring on long reward;

%EXTRACT FIRST AND LAST SOME TRIALS, DEFINED BY First_T and Last_T;
%First trials;
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

%last trials;
FR_bk1_sh_L = FR_bk1_sh(length(FR_bk1_sh(:,1))-Last_T+1:length(FR_bk1_sh(:,1)),:);
FR_bk1_lo_L = FR_bk1_lo(length(FR_bk1_lo(:,1))-Last_T+1:length(FR_bk1_lo(:,1)),:);

FR_bk2_sh_L = FR_bk2_sh(length(FR_bk2_sh(:,1))-Last_T+1:length(FR_bk2_sh(:,1)),:);
FR_bk2_lo_L = FR_bk2_lo(length(FR_bk2_lo(:,1))-Last_T+1:length(FR_bk2_lo(:,1)),:);

FR_bk3_sh_L = FR_bk3_sh(length(FR_bk3_sh(:,1))-Last_T+1:length(FR_bk3_sh(:,1)),:);
FR_bk3_lo_L = FR_bk3_lo(length(FR_bk3_lo(:,1))-Last_T+1:length(FR_bk3_lo(:,1)),:);

FR_bk4_sh_L = FR_bk4_sh(length(FR_bk4_sh(:,1))-Last_T+1:length(FR_bk4_sh(:,1)),:);
FR_bk4_lo_L = FR_bk4_lo(length(FR_bk4_lo(:,1))-Last_T+1:length(FR_bk4_lo(:,1)),:);

FR_bk5_sh_L = FR_bk5_sh(length(FR_bk5_sh(:,1))-Last_T+1:length(FR_bk5_sh(:,1)),:);
FR_bk5_lo_L = FR_bk5_lo(length(FR_bk5_lo(:,1))-Last_T+1:length(FR_bk5_lo(:,1)),:);

%COMBINE ALL CELLS;
cat_bk1_sh_F(:,:,cell_count) = cat(3,FR_bk1_sh_F);
cat_bk1_sh_L(:,:,cell_count) = cat(3,FR_bk1_sh_L);
cat_bk1_lo_F(:,:,cell_count) = cat(3,FR_bk1_lo_F);
cat_bk1_lo_L(:,:,cell_count) = cat(3,FR_bk1_lo_L);

cat_bk2_sh_F(:,:,cell_count) = cat(3,FR_bk2_sh_F);
cat_bk2_sh_L(:,:,cell_count) = cat(3,FR_bk2_sh_L);
cat_bk2_lo_F(:,:,cell_count) = cat(3,FR_bk2_lo_F);
cat_bk2_lo_L(:,:,cell_count) = cat(3,FR_bk2_lo_L);

cat_bk3_sh_F(:,:,cell_count) = cat(3,FR_bk3_sh_F);
cat_bk3_sh_L(:,:,cell_count) = cat(3,FR_bk3_sh_L);
cat_bk3_lo_F(:,:,cell_count) = cat(3,FR_bk3_lo_F);
cat_bk3_lo_L(:,:,cell_count) = cat(3,FR_bk3_lo_L);

cat_bk4_sh_F(:,:,cell_count) = cat(3,FR_bk4_sh_F);
cat_bk4_sh_L(:,:,cell_count) = cat(3,FR_bk4_sh_L);
cat_bk4_lo_F(:,:,cell_count) = cat(3,FR_bk4_lo_F);
cat_bk4_lo_L(:,:,cell_count) = cat(3,FR_bk4_lo_L);

cat_bk5_sh_F(:,:,cell_count) = cat(3,FR_bk5_sh_F);
cat_bk5_sh_L(:,:,cell_count) = cat(3,FR_bk5_sh_L);
cat_bk5_lo_F(:,:,cell_count) = cat(3,FR_bk5_lo_F);
cat_bk5_lo_L(:,:,cell_count) = cat(3,FR_bk5_lo_L);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Average all cells;
avg_bk1_sh_F = mean(cat_bk1_sh_F,3);
avg_bk1_sh_L = mean(cat_bk1_sh_L,3);
avg_bk1_lo_F = mean(cat_bk1_lo_F,3);
avg_bk1_lo_L = mean(cat_bk1_lo_L,3);

avg_bk2_sh_F = mean(cat_bk2_sh_F,3);
avg_bk2_sh_L = mean(cat_bk2_sh_L,3);
avg_bk2_lo_F = mean(cat_bk2_lo_F,3);
avg_bk2_lo_L = mean(cat_bk2_lo_L,3);

avg_bk3_sh_F = mean(cat_bk3_sh_F,3);
avg_bk3_sh_L = mean(cat_bk3_sh_L,3);
avg_bk3_lo_F = mean(cat_bk3_lo_F,3);
avg_bk3_lo_L = mean(cat_bk3_lo_L,3);

avg_bk4_sh_F = mean(cat_bk4_sh_F,3);
avg_bk4_sh_L = mean(cat_bk4_sh_L,3);
avg_bk4_lo_F = mean(cat_bk4_lo_F,3);
avg_bk4_lo_L = mean(cat_bk4_lo_L,3);

avg_bk5_sh_F = mean(cat_bk5_sh_F,3);
avg_bk5_sh_L = mean(cat_bk5_sh_L,3);
avg_bk5_lo_F = mean(cat_bk5_lo_F,3);
avg_bk5_lo_L = mean(cat_bk5_lo_L,3);

%Combine blocks;
avg_cis_sh_F = (avg_bk3_sh_F + avg_bk4_sh_F)/2;
avg_cis_sh_L = (avg_bk3_sh_L + avg_bk4_sh_L)/2;
avg_cis_lo_F = (avg_bk3_lo_F + avg_bk4_lo_F)/2;
avg_cis_lo_L = (avg_bk3_lo_L + avg_bk4_lo_L)/2;

avg_tra_sh_F = (avg_bk2_sh_F + avg_bk5_sh_F)/2;
avg_tra_sh_L = (avg_bk2_sh_L + avg_bk5_sh_L)/2;
avg_tra_lo_F = (avg_bk2_lo_F + avg_bk5_lo_F)/2;
avg_tra_lo_L = (avg_bk2_lo_L + avg_bk5_lo_L)/2;

%Combine early and late blocks;
all_bk1_sh = cat(1,avg_bk1_sh_F,avg_bk1_sh_L);
all_bk1_lo = cat(1,avg_bk1_lo_F,avg_bk1_lo_L);
all_bk2_sh = cat(1,avg_bk2_sh_F,avg_bk2_sh_L);
all_bk2_lo = cat(1,avg_bk2_lo_F,avg_bk2_lo_L);
all_bk3_sh = cat(1,avg_bk3_sh_F,avg_bk3_sh_L);
all_bk3_lo = cat(1,avg_bk3_lo_F,avg_bk3_lo_L);
all_bk4_sh = cat(1,avg_bk4_sh_F,avg_bk4_sh_L);
all_bk4_lo = cat(1,avg_bk4_lo_F,avg_bk4_lo_L);
all_bk5_sh = cat(1,avg_bk5_sh_F,avg_bk5_sh_L);
all_bk5_lo = cat(1,avg_bk5_lo_F,avg_bk5_lo_L);

all_cis_sh = cat(1,avg_cis_sh_F,avg_cis_sh_L);
all_cis_lo = cat(1,avg_cis_lo_F,avg_cis_lo_L);
all_tra_sh = cat(1,avg_tra_sh_F,avg_tra_sh_L);
all_tra_lo = cat(1,avg_tra_lo_F,avg_tra_lo_L);


x1 = [12 12]; y1 = [0 First_T+Last_T+1];
x2 = [37 37]; y2 = [0 First_T+Last_T+1];
x3 = [57 57]; y3 = [0 First_T+Last_T+1];
x4 = 1:length(bin_centers);
y4(1:length(bin_centers)) = First_T+0.5; 


%MAKE FIGURES;
%figure 1; block1;
figure;
subplot(1,2,1);
image(all_bk1_sh*mag); hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig1. block1 short');
subplot(1,2,2);
image(all_bk1_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig1. block1 long');

%figure 2; block2;
figure;
subplot(1,2,1);
image(all_bk2_sh*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig2. block2 short');
subplot(1,2,2);
image(all_bk2_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig2. block2 long');

%figure 3; block3;
figure;
subplot(1,2,1);
image(all_bk3_sh*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig3. block3 short');
subplot(1,2,2);
image(all_bk3_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig3. block3 long');

%figure 4; block4;
figure;
subplot(1,2,1);
image(all_bk4_sh*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig4. block4 short');
subplot(1,2,2);
image(all_bk4_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig4. block4 long');

%figure 5; block5;
figure;
subplot(1,2,1);
image(all_bk5_sh*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig5. block5 short');
subplot(1,2,2);
image(all_bk5_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig5. block5 long');

%figure 6; delay only;
figure;
subplot(1,2,1);
image(all_cis_sh*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig6. delay only short');
subplot(1,2,2);
image(all_cis_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig6. delay only long');

%figure 7; trans;
figure;
subplot(1,2,1);
image(all_tra_sh*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig7. delay-flavor short');
subplot(1,2,2);
image(all_tra_lo*mag);hold on;
plot(x1,y1,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x2,y2,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x3,y3,'LineStyle',':','Color',[0 0 0],'LineWidth',1);hold on;
plot(x4,y4,'LineStyle','-','Color',[0 0 0],'LineWidth',1);hold on;
title('Fig7. delay_flavors long');