clear;
addpath('F:\Jerome\CvMHT_baseline2.0\');

addpath('./unti/');
seqs = configSeqs_benchmark;
views = {'t','h1','h2', 'h3','h4'};

Res_path = 'F:/Jerome/Evaluation/Result_all/GMMCP_CvMHATB/';
GT_path = './Eval_Data/GT_txt/';
newRes_path = './Eval_Data/All_Res_txt/GMMCP/';


for seq_i = 1:length(seqs)
    
    scene_name = seqs{seq_i}.name; % 'V1-S_square-G_3';
    num_hor_seq =  seqs{seq_i}.num_hor; 
    maxframe = seqs{seq_i}.endFrame;
    all_view_num_seq = num_hor_seq + 1;
    
    res_sort = cell(1,all_view_num_seq);
    gt_sort = cell(1,all_view_num_seq);
 
    for view_i = 1:all_view_num_seq
  
        sequenceName = strcat(scene_name,'_',views{view_i});
     
        resFilename = [Res_path, sequenceName, '.txt'];
        gtFilename = [GT_path, sequenceName, '.txt'];
        res = dlmread(resFilename);
        gt = dlmread(gtFilename);
        
        res = res(:,1:6);
        gt = gt(:,1:6);
        
        res_sort{view_i} = sortrows(res,1);
        gt_sort{view_i} = sortrows(gt,1);
    
    end
   
        
    for first_view = 1 %: all_view_num_seq - 1
        for last_view = first_view + 1 : all_view_num_seq
            top_res = res_sort{first_view};
            hor_res = res_sort{last_view};
            top_gt = gt_sort{first_view};
            hor_gt = gt_sort{last_view};

            fid_top = fopen([newRes_path,scene_name,'_',views{first_view},'.txt'],'w');
            fid_hor = fopen([newRes_path,scene_name,'_',views{last_view},'.txt'],'w');
            
            top_id_init = -1;
            hor_id_init = -1;
            
            for frm_i = 1 : maxframe
                if (seq_i == 26 && frm_i >= 1501) || (seq_i == 29 && frm_i >= 1301) || (seq_i == 30 && frm_i >= 1201)
                    break
                end
                
                top_res_i = top_res(top_res(:,1) == frm_i,:);
                top_gt_i = top_gt(top_gt(:,1) == frm_i,:);
                hor_res_i = hor_res(hor_res(:,1) == frm_i,:);
                hor_gt_i = hor_gt(hor_gt(:,1) == frm_i,:);

                if frm_i == 1
                    thr_dis = 20;
                    thr_IOU = 0.5;
                    top_id_gt = top_gt_i(:,2);
                    hor_id_gt = hor_gt_i(:,2);
                    top_id_res = top_res_i(:,2);
                    hor_id_res = hor_res_i(:,2);
                    top_match = top_match_box(top_res_i,top_gt_i,thr_dis);     
                    [top_gt,top_det] = find(top_match == 1);
                    top_ids = containers.Map(top_id_res(top_det),top_id_gt(top_gt));
                    top_ids_convert = containers.Map(top_id_gt(top_gt),top_id_res(top_det));
                    hor_match = hor_match_box(hor_res_i,hor_gt_i,thr_IOU);
                    [hor_gt,hor_det] = find(hor_match == 1);
                    hor_ids = containers.Map(hor_id_res(hor_det),hor_id_gt(hor_gt));
                    hor_ids_convert = containers.Map(hor_id_gt(hor_gt),hor_id_res(hor_det));
                end

                top_id_res = top_res_i(:,2);
                hor_id_res = hor_res_i(:,2);
                [top_id_res_num,~] = size(top_id_res);
                new_top_id_res = [];
                new_hor_id_res = [];
                for top_i = 1:top_id_res_num
                    if top_ids.isKey(top_id_res(top_i)) ~= 1 && top_ids_convert.isKey(top_id_res(top_i)) ~= 1
                        new_top_id_res(top_i) = top_id_res(top_i);
                    elseif top_ids.isKey(top_id_res(top_i)) ~= 1 && top_ids_convert.isKey(top_id_res(top_i)) == 1
                        new_top_id_res(top_i) = top_id_init;   
                        top_id_init = top_id_init - 1;
                    else
                        new_top_id_res(top_i) = top_ids(top_id_res(top_i));
                    end
                end

                [hor_id_res_num,~] = size(hor_id_res);
                for hor_i = 1:hor_id_res_num
                    if hor_ids.isKey(hor_id_res(hor_i)) ~= 1 && hor_ids_convert.isKey(hor_id_res(hor_i)) ~= 1
                        new_hor_id_res(hor_i) = hor_id_res(hor_i);
                    elseif hor_ids.isKey(hor_id_res(hor_i)) ~= 1 && hor_ids_convert.isKey(hor_id_res(hor_i)) == 1    
                        new_hor_id_res(hor_i) = hor_id_init;
                        hor_id_init = hor_id_init - 1;
                    else
                        new_hor_id_res(hor_i) = hor_ids(hor_id_res(hor_i));
                    end
                end
                top_res_i(:,2) = new_top_id_res;
                hor_res_i(:,2) = new_hor_id_res;
                for top_i = 1:top_id_res_num 
                    new_top_res_i = top_res_i(top_i,:);
                    fprintf(fid_top,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',[new_top_res_i,-1,-1,-1,-1]);
                end
                for hor_i = 1:hor_id_res_num
                    new_hor_res_i = hor_res_i(hor_i,:);
                    fprintf(fid_hor,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',[new_hor_res_i,-1,-1,-1,-1]);
                end
        
            end
            fclose(fid_top);
            fclose(fid_hor);
        end
    end

end

