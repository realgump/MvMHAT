clear;
addpath('/home/jerome/Documents/CVMOT/CVMOT_Tracker6.0');
addpath('./unti/');

seqs = configSeqs;
% seqs = seqs(8);

allSequences = cell(1,2*length(seqs));

for seq_i = 1 :  length(seqs)

top_i =2 *(seq_i)-1;
allSequences{top_i} =  strcat(seqs{seq_i}.name,'_top');

hor_i =2 *(seq_i);
allSequences{hor_i} =  strcat(seqs{seq_i}.name,'_hor');

end

Sequences = allSequences(2:2:end);

Res_path = './Eval_Data/All_Res_txt/Res_app_11/';
GT_path = './Eval_Data/GT_txt/';

tracker = {'DMAN','GMMCP','MDP','Ours'};

compareMat = zeros(length(Sequences), 600);

for ind = 1:1:length(Sequences)
    
    sequenceName = char(Sequences(ind));

    resFilename = [Res_path, sequenceName, '.txt'];
    gtFilename = [GT_path, sequenceName, '.txt'];
    res = dlmread(resFilename);
    gt = dlmread(gtFilename);

    res = res(:,1:6);
    gt = gt(:,1:6);

    res_sort = sortrows(res,1);
    gt_sort = sortrows(gt,1);
    gt_first_id = [];
    res_first_id = [];
    for frm_i = 1:600
        res_i = res_sort(res_sort(:,1) == frm_i,:);
        gt_i = gt_sort(gt_sort(:,1) == frm_i,:);
        if frm_i == 1
%             %%%%%%%%%%%%top
%             thr_dis = 20;
%             gt_i_id = gt_i(:,2);
%             res_i_id = res_i(:,2);
%             top_match = top_match_box(res_i,gt_i,thr_dis);     
%             [top_gt,top_det] = find(top_match == 1);
% %             top_ids = containers.Map(res_i_id(top_det),gt_i_id(top_gt));
% %             top_ids_convert = containers.Map(gt_i_id(top_gt),res_i_id(top_det));
%             gt_first_id = gt_i_id(top_gt)';
%             res_first_id = res_i_id(top_det)';
            %%%%%%%%%%%%hor
            thr_IOU = 0.5;
            gt_i_id = gt_i(:,2);
            res_i_id = res_i(:,2);            
            hor_match = hor_match_box(res_i,gt_i,thr_IOU);
            [hor_gt,hor_det] = find(hor_match == 1);
%             hor_ids = containers.Map(res_i_id(hor_det),gt_i_id(hor_gt));
%             hor_ids_convert = containers.Map(gt_i_id(hor_gt),res_i_id(hor_det));
            gt_first_id = gt_i_id(hor_gt)';
            res_first_id = res_i_id(hor_det)';
        end
        gt_i_id = gt_i(:,2)';
        res_i_id = res_i(:,2)';
        gt_true_num = length(intersect(gt_i_id,gt_first_id));
        res_true_num = length(intersect(res_i_id,res_first_id));
        compareMat(ind, frm_i) = min(res_true_num, gt_true_num) / gt_true_num;

    end    

end
compareMat = mean(compareMat);
a = 1;
        
