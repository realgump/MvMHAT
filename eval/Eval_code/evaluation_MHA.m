function eval_res = evaluation_MHA(seqs,Res_path,GT_path)

addpath('./unti/');
views = {'1','2', '3','4'};
for seq_i = 1:length(seqs)
    scene_name = seqs{seq_i}.name; % 'V1-S_square-G_3';
    num_hor_seq =  seqs{seq_i}.num_hor; 
    start_frame = seqs{seq_i}.startFrame;
    maxframe = seqs{seq_i}.endFrame;
    numframe = maxframe - start_frame + 1;
    all_view_num_seq = num_hor_seq;
%     frame_nums = frames(2) - frames(1);
    
    res_sort = cell(1,all_view_num_seq);
    gt_sort = cell(1,all_view_num_seq);

    for view_i = 1:all_view_num_seq
  
        sequenceName = strcat(scene_name,'_',views{view_i});
     
        resFilename = [Res_path, sequenceName, '.txt'];
        gtFilename = [GT_path, sequenceName, '.txt'];
        res = dlmread(resFilename);
        gt = dlmread(gtFilename);
        
%         mask1 = gt(:, 1) > frames(1);
%         mask2 = gt(:, 1) < frames(2);
%         mask = mask1 & mask2;
%         gt = gt(mask(:,1), :);
        
        res = res(:,1:6);
        gt = gt(:,1:6);
        
        res_sort{view_i} = sortrows(res,1);
        gt_sort{view_i} = sortrows(gt,1);
    
    end
    
    all_bbx = zeros(numframe,1);
    precision = zeros(numframe,1);
    recall = zeros(numframe,1);
    %%% new_m, new_fp, new_mme
    new_missed = zeros(numframe,1);
    new_missed_avgview = zeros(numframe,1);
    new_fp = zeros(numframe,1);
    new_mme = zeros(numframe,1);
    
      
    %% Cn2 cal multi-view metrix (top----view1, hor-----view2)
    for frm_i_all = start_frame : maxframe
      frm_i =  frm_i_all - start_frame + 1;
%     for frm_i = frames(1) : frames(2)
       %% pair view Cn2
%         if frm_i == 10
%            a = 1 
%         end
        all_TP = 0;
        all_recall_denominator = 0;
        all_prec_denominator = 0;
        bbx_tmp = 0;
              
        for first_view = 1 : all_view_num_seq
            for last_view = first_view + 1 : all_view_num_seq
                top_res = res_sort{first_view};
                hor_res = res_sort{last_view};
                top_gt = gt_sort{first_view};
                hor_gt = gt_sort{last_view};

                top_res_i = top_res(top_res(:,1) == frm_i_all,:);
                top_gt_i = top_gt(top_gt(:,1) == frm_i_all,:);
                hor_res_i = hor_res(hor_res(:,1) == frm_i_all,:);
                hor_gt_i = hor_gt(hor_gt(:,1) == frm_i_all,:);

                thr_dis = 20;
                thr_IOU = 0.5;
                
                %%gt null, res not null
                if size(hor_gt_i,1) == 0
%                     new_fp(frm_i) = new_fp(frm_i) + size(hor_res_i, 1);
                    all_prec_denominator = all_prec_denominator + size(hor_res_i, 1);
                    continue;
                end
                if size(top_gt_i,1) == 0
%                     new_fp(frm_i) = new_fp(frm_i) + size(top_res_i, 1);
                    all_prec_denominator = all_prec_denominator + size(top_res_i, 1);
                    continue;
                end
                             
              %% view1 is top or view1 is hor
%                 if first_view == 1
%                     top_match = top_match_box(top_res_i,top_gt_i,thr_dis);
%                 else
%                     top_match = hor_match_box(top_res_i,top_gt_i,thr_IOU);
%                 end

                top_match = hor_match_box(top_res_i,top_gt_i,thr_IOU);
                top_id_gt = top_gt_i(:,2);
                top_id_res = top_res_i(:,2);

                det2gt_top = zeros(size(top_id_gt,1),1);
                [gt,det] = find(top_match == 1);
                if ~isempty(gt) && ~isempty(det)
                    try
                        top_ids_gt = containers.Map(top_id_gt(gt),top_id_res(det));
                    catch err
                    end
                else
                    top_ids_gt = [];
                end
        %         top_ids_gt = containers.Map(top_id_gt(gt),top_id_res(det));
                det2gt_top(gt) =  top_id_res(det);  
%                 det2gt_top(det) =  top_id_gt(gt);%%%%label ID
                lab_top_end = det2gt_top'; 
                               
                hor_match = hor_match_box(hor_res_i,hor_gt_i,thr_IOU);
                hor_id_gt = hor_gt_i(:,2);
                hor_id_res = hor_res_i(:,2);

                det2gt = zeros(size(hor_id_gt,1),1);
                [gt,det] = find(hor_match == 1);

                if ~isempty(gt) && ~isempty(det)
                    try
                        hor_ids_gt = containers.Map(hor_id_gt(gt),hor_id_res(det));
                    catch err
                    end
                else
                    hor_ids_gt = [];
                end

                det2gt(gt) =  hor_id_res(det);   
%                 det2gt(det) =  hor_id_gt(gt); 
                lab_ego_end = det2gt';  

                %%%%%%%% --------- new_fp-------
                diss_res_top = setdiff(top_id_res,det2gt_top);
                if ~isempty(diss_res_top)
                    diss_top_num = size(diss_res_top,1);
                    for diss_i = 1:diss_top_num
                        diss_label_i = diss_res_top(diss_i);
                        diss_index = find(hor_id_res == diss_label_i, 1);
                        if ~isempty(diss_index)
                            new_fp(frm_i) = new_fp(frm_i) + 1;
                        end
                    end
                end

                diss_res_hor = setdiff(hor_id_res,det2gt);
                if ~isempty(diss_res_hor)
                    diss_hor_num = size(diss_res_hor,1);
                    for diss_i = 1:diss_hor_num
                        diss_label_i = diss_res_hor(diss_i);
                        diss_index = find(top_id_res == diss_label_i, 1);
                        if ~isempty(diss_index)
                            new_fp(frm_i) = new_fp(frm_i) + 1;
                        end
                    end
                end

                [gt_match_id,top_id_gt_trans,hor_id_gt_trans] = transform_top_hor(top_id_gt,hor_id_gt);  % transform the unmatched id into 0
                [res_match_id,top_id_res_trans,hor_id_res_trans] = transform_top_hor(lab_top_end,lab_ego_end);

                %%%%%%-------new_metric-----
                new_top_id_gt_trans = top_id_gt_trans;
                new_hor_id_gt_trans = hor_id_gt_trans;
                gt_match_num = size(gt_match_id,1);
                res_match_num = size(res_match_id,2);
                bbx = size(top_id_gt,1) + size(hor_id_gt,1);
%                 all_bbx = all_bbx + bbx;
                bbx_tmp = bbx_tmp + bbx;
                %%%%%%%% --------- new_missed-------
%                 new_top_id_gt_trans(find(new_top_id_gt_trans==0)) = [];
%                 new_hor_id_gt_trans(find(new_hor_id_gt_trans==0)) = [];
%                 [~,top_gt2res_num] = size(new_top_id_gt_trans);
%                 for top_i = 1:top_gt2res_num
%                     %%%%% miss detection
%                     if ~isempty(top_ids_gt)
%                         if top_ids_gt.isKey(new_top_id_gt_trans(top_i)) ~= 1
%                             new_missed(frm_i) = new_missed(frm_i) + 1;
%                         else
%                             res_top_label = top_ids_gt(new_top_id_gt_trans(top_i));
%                             if find(hor_id_res == res_top_label) == 0
%                                 new_mme(frm_i) = new_mme(frm_i) + 1;
% %                                 new_missed(frm_i) = new_missed(frm_i) + 1;
%                             end
%                         end
%                     else
%                         new_missed(frm_i) = new_missed(frm_i) + size(new_top_id_gt_trans,2); 
%                     end
%                 end
%                 [~,hor_gt2res_num] = size(new_hor_id_gt_trans);
%                 for hor_i = 1:hor_gt2res_num
%                     %%%%% miss detection
%                     if ~isempty(hor_ids_gt)
%                         if hor_ids_gt.isKey(new_hor_id_gt_trans(hor_i)) ~= 1
%                             new_missed(frm_i) = new_missed(frm_i) + 1;
%                         else
%                             res_hor_label = hor_ids_gt(new_hor_id_gt_trans(hor_i));
%                             %%%%% miss match            
%                             if find(top_id_res == res_hor_label) == 0
%                                 new_mme(frm_i) = new_mme(frm_i) + 1;
% %                                 new_missed(frm_i) = new_missed(frm_i) + 1;
%                             end
%                         end  
%                     else
%                         new_missed(frm_i) = new_missed(frm_i) + size(new_hor_id_gt_trans,2);
%                     end
%                 end
                %%%%%%

                % ----top--------find the match number in GT & lab --------------------
                match_idx = intersect(find(top_id_gt_trans ~= 0) , find(top_id_res_trans ~= 0));
                match_top_GT = top_id_gt_trans(match_idx);
                match_top_final = top_id_res_trans(match_idx);

                % ----ego--------find the match number in GT & lab --------------------
                match_idx_ego = intersect(find(hor_id_gt_trans ~= 0), find(hor_id_res_trans ~= 0));   
                match_ego_GT = hor_id_gt_trans(match_idx_ego);
                match_ego_final = hor_id_res_trans(match_idx_ego);

                % ------------change the num and order of ego and top as the intersect of them-----
                match_inter_GT = intersect(match_ego_GT,match_top_GT);   %return match gt label
                [~,pos1]= ismember(match_inter_GT,match_top_GT);
                [~,pos2]= ismember(match_inter_GT,match_ego_GT);   
                match_top_GT = match_top_GT(pos1);
                match_top_final = match_top_final(pos1);
                match_ego_GT = match_ego_GT(pos2);
                match_ego_final = match_ego_final(pos2);

                diff_match_top = match_top_GT - match_top_final;
                diff_match_ego = match_ego_GT - match_ego_final;

                match_err = diff_match_top - diff_match_ego;        
                right = length(match_err(match_err==0));

                new_mme(frm_i) = new_mme(frm_i) + res_match_num - right;
                new_missed(frm_i) = new_missed(frm_i) + (gt_match_num - res_match_num) * 2;
                TP = right;
                
                all_TP = all_TP + TP;
                all_recall_denominator = all_recall_denominator + length(top_id_gt_trans(top_id_gt_trans > 0));
                all_prec_denominator = all_prec_denominator + length(top_id_res_trans(top_id_res_trans > 0));
%                 recall(frm_i) = TP/length(top_id_gt_trans(top_id_gt_trans > 0));
%                 if length(top_id_res_trans(top_id_res_trans > 0)) == 0
%                     precision(frm_i) = 0;
%                 else
%                     precision(frm_i) = TP/length(top_id_res_trans(top_id_res_trans > 0));
%                 end

                %accuracy = 2 * TP/(length(lab_top_GT) + length(lab_ego_GT));
                %eval_res.acc = accuracy;

        %       top_id_res_match = top_id_res(top_det_match);
        %       top_id_gt_match = top_id_gt(top_gt_match);
            end
        end
        recall(frm_i) = all_TP/all_recall_denominator;
        precision(frm_i) = all_TP/all_prec_denominator;
        new_missed(frm_i) = new_missed(frm_i); %/ max(all_view_num_seq * (all_view_num_seq - 1) / 2, 1); %/ bbx;% / all_bbx;
        new_missed_avgview(frm_i) = new_missed(frm_i) / (all_view_num_seq); %/ max(all_view_num_seq * (all_view_num_seq - 1) / 2, 1); %/ bbx;% / all_bbx;
        new_fp(frm_i) = new_fp(frm_i); %/ max(all_view_num_seq * (all_view_num_seq - 1) / 2, 1); %/ bbx;% / all_bbx;
        new_mme(frm_i) = new_mme(frm_i); %/ max(all_view_num_seq * (all_view_num_seq - 1) / 2, 1) ;% / bbx;% / all_bbx;
        all_bbx(frm_i) = all_bbx(frm_i) + bbx_tmp; %/ max(all_view_num_seq * (all_view_num_seq - 1) / 2, 1);
 
    end
    
        eval_res.pre{seq_i} = precision;
        eval_res.rec{seq_i} = recall;
        eval_res.new_missed{seq_i} = new_missed;
        eval_res.new_fp{seq_i} = new_fp;
        eval_res.new_mme{seq_i} = new_mme;
        eval_res.all_bbx{seq_i} = all_bbx;
        eval_res.new_missed_avgview{seq_i} = new_missed_avgview;
end

prec_avg = zeros(1,length(eval_res.pre));
reca_avg = zeros(1,length(eval_res.pre));
missed_avgview = zeros(1,length(eval_res.pre));
new_missed_avg = zeros(1,length(eval_res.pre));
new_fp_avg = zeros(1,length(eval_res.pre));
new_mme_avg = zeros(1,length(eval_res.pre));
new_bbox_avg = zeros(1,length(eval_res.pre));

prec_1 = zeros(1,length(eval_res.pre)); 
reca_1 = zeros(1,length(eval_res.pre));
% CVMOTA = zeros(1,length(eval_res.pre));

for i = 1 : length(eval_res.pre)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    prec = eval_res.pre{i};
    reca = eval_res.rec{i};
    a = (prec == 1);
    num = sum(a(:));
    num_frame = size(prec,1);
    prec_1(i) = num / num_frame;
    b = (reca == 1);
    num_rec = sum(b(:));
    reca_1(i) = num_rec / num_frame;
    new_missed = eval_res.new_missed{i};
    new_fp = eval_res.new_fp{i};
    new_mme = eval_res.new_mme{i};
    new_bbox = eval_res.all_bbx{i};
    prec((isnan(prec)==1)) = 0;
    reca((isnan(reca)==1)) = 0;
    new_missed((isnan(new_missed)==1)) = 0;
    new_fp((isnan(new_fp)==1)) = 0;
    new_mme((isnan(new_mme)==1)) = 0;
    new_bbox((isnan(new_bbox)==1)) = 0;
    prec_avg(i) = mean(prec);
    reca_avg(i) = mean(reca);
%     new_missed_avg(i) = mean(new_missed);
%     new_fp_avg(i) = mean(new_fp);
%     new_mme_avg(i) = mean(new_mme);
    new_missed_avg(i) = sum(new_missed);
    new_fp_avg(i) = sum(new_fp);
    new_mme_avg(i) = sum(new_mme);
    new_bbox_avg(i) = sum(new_bbox);
%     CVMOTA(i) = 1 - (2*(new_mme_avg) + (new_fp_avg) + new_missed_avg)/sum(new_bbox);
    missed_avgview(i) = sum(eval_res.new_missed_avgview{i});
end
% CVMOTA = mean(CVMOTA);
CVMOTA = 1 - (2*sum(new_mme_avg) + sum(new_fp_avg) + sum(new_missed_avg))/sum(new_bbox_avg);
CVIDP = mean(prec_avg);
CVIDR = mean(reca_avg);
CVIDF1 = (2*CVIDP*CVIDR)/(CVIDP + CVIDR);
disp(['MVIDP----',num2str(CVIDP)])
disp(['MVIDR----',num2str(CVIDR)])
disp(['MVIDF1----',num2str(CVIDF1)])
% disp(['MVFP----',num2str(sum(new_fp_avg))])
% disp(['MVMME----',num2str(sum(new_mme_avg))])
disp(['MVMISSED----',num2str(sum(new_missed_avg))])
% disp(['ALLBOX----',num2str(all_bbx)])
disp(['MVMA----',num2str(CVMOTA)])
disp([num2str(CVIDP), ' ', num2str(CVIDR), ' ', num2str(CVIDF1), ' ', num2str(sum(new_missed_avg)), ' ', num2str(CVMOTA)])
% mean(prec_avg)
% mean(reca_avg)
% disp(CVMOTA)
% disp(CVIDF1)
% mean(prec_1)
% mean(reca_1)


