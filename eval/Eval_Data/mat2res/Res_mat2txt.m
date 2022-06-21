for res = 11
Restype = ['GMMCP_CvMHATB_gt']; %,num2str(res)];
% Restype = ['Res_',num2str(res)];

Trackertype = 'CvMHT_baseline2.0';

% filepath = ['F:\Jerome\',Trackertype,'\Data\Tracking_Res\'];
filepath = 'F:/Jerome/Evaluation/Result_all/';

addpath('F:\Jerome\CvMHT_baseline2.0');


seq = configSeqs_benchmark;

for seq_i = 1 : length(seq)
    
scene_name = seq{seq_i}.name; % 'V1-S_square-G_3';
num_hor_seq =  seq{seq_i}.num_hor; 

views = {'t','h1','h2', 'h3','h4'};

    for view_i = 1 : num_hor_seq + 1

    view = views{view_i};

%     res_mat = cell2mat(struct2cell(load([filepath,Restype,'\track_res_',scene_name,'_',view,'.mat'])));
    res_mat = cell2mat(struct2cell(load([filepath,Restype,'\',scene_name,'_',view,'.mat'])));
    
    if ~exist(['../All_Res_txt/',Restype],'dir')
        mkdir(['../All_Res_txt/',Restype]);
    end

    fid = fopen(['../All_Res_txt/',Restype,'/',scene_name,'_',view,'.txt'],'w');

        for i = 1 : length(res_mat)

        clear res_mat_i
        res_mat_i(1:4) =  res_mat(i,1:4);
        res_i_frm = res_mat_i(1);
        if (seq_i == 26 && res_i_frm >= 1501) || (seq_i == 29 && res_i_frm >= 1301) || (seq_i == 30 && res_i_frm >= 1201)
            continue
        end
        res_mat_i(5) =  res_mat(i,5) - res_mat(i,3);
        res_mat_i(6) =  res_mat(i,6) - res_mat(i,4);

        fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',[res_mat_i,-1,-1,-1,-1]);

        end

    fclose(fid);
    
    end

end

end