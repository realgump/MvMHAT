function eval_res = evaluation_MHA(allSequences,Res_path,GT_path)

addpath('./unti/');
views = {'top','hor'};
seq_i = 1;
all_bbx = 0;
for ind = 1:2:length(allSequences)
    
    res_sort = cell(1,2);
    gt_sort = cell(1,2);
 
    for view_i = 1:2
        idx = ind + view_i - 1;   
        sequenceName = char(allSequences(idx));
     
        resFilename = [Res_path, sequenceName, '.txt'];
        gtFilename = [GT_path, sequenceName, '.txt'];
        res = dlmread(resFilename);
        gt = dlmread(gtFilename);
        
        res = res(:,1:6);
        gt = gt(:,1:6);
        
        res_sort{view_i} = sortrows(res,1);
        gt_sort{view_i} = sortrows(gt,1);
    
    end
    
    top_res = res_sort{1};
    hor_res = res_sort{2};
    top_gt = gt_sort{1};
    hor_gt = gt_sort{2};
    
    maxframe = max([top_res(end,1),hor_res(end,1),top_gt(end,1),hor_gt(end,1)]);
    precision = zeros(maxframe,1);
    recall = zeros(maxframe,1);
    %%% new_m, new_fp, new_mme
    new_missed = zeros(maxframe,1);
    new_fp = zeros(maxframe,1);
    new_mme = zeros(maxframe,1);
    bbx_num = zeros(maxframe,1);
    
    for frm_i = 1 : maxframe
        top_res_i = top_res(top_res(:,1) == frm_i,:);
        top_gt_i = top_gt(top_gt(:,1) == frm_i,:);
        hor_res_i = hor_res(hor_res(:,1) == frm_i,:);
        hor_gt_i = hor_gt(hor_gt(:,1) == frm_i,:);
        
        thr_dis = 20;
        thr_IOU = 0.5;
        top_match = top_match_box(top_res_i,top_gt_i,thr_dis);
        top_id_gt = top_gt_i(:,2);
        top_id_res = top_res_i(:,2);
                         
        det2gt_top = zeros(size(top_id_gt,1),1);
        [gt,det] = find(top_match == 1);
        top_ids_gt = containers.Map(top_id_gt(gt),top_id_res(det));
        det2gt_top(gt) =  top_id_res(det);   
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
        bbx = size(top_id_gt,1) + size(hor_id_gt,1);
        bbx_num(frm_i) = bbx;
        all_bbx = all_bbx + bbx;
        %%%%%%%% --------- new_missed-------
        new_top_id_gt_trans(find(new_top_id_gt_trans==0)) = [];
        new_hor_id_gt_trans(find(new_hor_id_gt_trans==0)) = [];
        [~,top_gt2res_num] = size(new_top_id_gt_trans);
        for top_i = 1:top_gt2res_num
            %%%%% miss detection
            if top_ids_gt.isKey(new_top_id_gt_trans(top_i)) ~= 1
                new_missed(frm_i) = new_missed(frm_i) + 1;
            else
                res_top_label = top_ids_gt(new_top_id_gt_trans(top_i));
                if find(hor_id_res == res_top_label) == 0
                    new_missed(frm_i) = new_missed(frm_i) + 1;
                end
            end
        end
        [~,hor_gt2res_num] = size(new_hor_id_gt_trans);
        for hor_i = 1:hor_gt2res_num
            %%%%% miss detection
            if ~isempty(hor_ids_gt)
                if hor_ids_gt.isKey(new_hor_id_gt_trans(hor_i)) ~= 1
                    new_missed(frm_i) = new_missed(frm_i) + 1;
                else
                    res_hor_label = hor_ids_gt(new_hor_id_gt_trans(hor_i));
                    %%%%% miss match            
                    if find(top_id_res == res_hor_label) == 0
                        new_missed(frm_i) = new_missed(frm_i) + 1;
                    end
                end  
            else
                new_missed(frm_i) = new_missed(frm_i) + size(new_hor_id_gt_trans,2);
            end
        end
        
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
        
        new_mme(frm_i) = gt_match_num - right;
        TP = right;

        recall(frm_i) = TP/length(top_id_gt_trans(top_id_gt_trans > 0));
        if length(top_id_res_trans(top_id_res_trans > 0)) == 0
            precision(frm_i) = 0;
        else
            precision(frm_i) = TP/length(top_id_res_trans(top_id_res_trans > 0));
        end
        
        new_missed(frm_i) = new_missed(frm_i);% / all_bbx;
        new_fp(frm_i) = new_fp(frm_i);% / all_bbx;
        new_mme(frm_i) = new_mme(frm_i);% / all_bbx;
        %accuracy = 2 * TP/(length(lab_top_GT) + length(lab_ego_GT));
        %eval_res.acc = accuracy;

%       top_id_res_match = top_id_res(top_det_match);
%       top_id_gt_match = top_id_gt(top_gt_match);
    end
    
        eval_res.pre{seq_i} = precision;
        eval_res.rec{seq_i} = recall;
        eval_res.new_missed{seq_i} = new_missed;
        eval_res.new_fp{seq_i} = new_fp;
        eval_res.new_mme{seq_i} = new_mme;
        eval_res.bbx_num{seq_i} = bbx_num;
        seq_i = seq_i + 1;
end

prec_avg = zeros(1,length(eval_res.pre));
reca_avg = zeros(1,length(eval_res.pre));
new_missed_avg = zeros(1,length(eval_res.pre));
new_fp_avg = zeros(1,length(eval_res.pre));
new_mme_avg = zeros(1,length(eval_res.pre));

CVIDF1_seq_all = zeros(600,1);
CVMA_seq_all = zeros(600,1);

for i = 1 : length(eval_res.pre)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    prec = eval_res.pre{i};
    reca = eval_res.rec{i};
    new_missed = eval_res.new_missed{i};
    new_fp = eval_res.new_fp{i};
    new_mme = eval_res.new_mme{i};
    bbx_num = eval_res.bbx_num{i};
    prec((isnan(prec)==1)) = 0;
    reca((isnan(reca)==1)) = 0;
    new_missed((isnan(new_missed)==1)) = 0;
    new_fp((isnan(new_fp)==1)) = 0;
    new_mme((isnan(new_mme)==1)) = 0;
    prec_avg(i) = mean(prec);
    reca_avg(i) = mean(reca);
    new_missed_avg(i) = sum(new_missed);
    new_fp_avg(i) = sum(new_fp);
    new_mme_avg(i) = sum(new_mme);
    
    CVIDF1_seq = (2 * prec .* reca)./(prec + reca);
    CVIDF1_seq(isnan(CVIDF1_seq(:))==1)=0;
    CVIDF1_seq_all = CVIDF1_seq(1:600) + CVIDF1_seq_all;
    CVMA_seq = 1 - (new_missed+new_fp+2*new_mme)./bbx_num;
    CVMA_seq(isnan(CVMA_seq(:))==1)=0;
    CVMA_seq_all = CVMA_seq_all + CVMA_seq(1:600);
    
end

plotDrawStyle={
struct('color',[0,0,0],'lineStyle','--'),...%dark red
struct('color',[0,1,0],'lineStyle','--'),...%orange
struct('color',[0,0,1],'lineStyle','-.'),...%Turquoise
struct('color',[1,0,0],'lineStyle',':'),...%purple    %%%%%%%%%%%%%%%%%%%%
};
point = {'^','d','s','p'};

CVIDF1_seq_avg = CVIDF1_seq_all/length(eval_res.pre);

figure(1)
xx = 1:1:600;
y_all = CVIDF1_seq_avg';
x = 1:50:600;
p=polyfit(xx,y_all,3); %����3Ϊ����ݵ�ĸ���-1��pΪ���ض���ʽ��ϵ��
yy=polyval(p,xx); %yyΪ��xx���и�ݶ���ʽ����õ���Ԥ��ֵ���У�
plot(xx,yy,'color',plotDrawStyle{4}.color,'lineStyle', plotDrawStyle{4}.lineStyle,'LineWidth',1.5); hold on;
plot(x,yy(x),'color',plotDrawStyle{4}.color,'Marker',point{4});
% plot(1:600,CVIDF1_seq_avg);
axis([0 600 0 1])
CVMA_seq_avg = CVMA_seq_all/length(eval_res.pre);
figure(2)
y_all = CVMA_seq_avg';
p=polyfit(xx,y_all,3); %����3Ϊ����ݵ�ĸ���-1��pΪ���ض���ʽ��ϵ��
yy=polyval(p,xx); %yyΪ��xx���и�ݶ���ʽ����õ���Ԥ��ֵ���У�
plot(xx,yy,'color',plotDrawStyle{4}.color,'lineStyle', plotDrawStyle{4}.lineStyle,'LineWidth',1.5); hold on;
plot(x,yy(x),'color',plotDrawStyle{4}.color,'Marker',point{4});
% plot(1:600,CVMA_seq_avg);
axis([0 600 0 1])

CVMOTA = 1 - (2*sum(new_mme_avg) + sum(new_fp_avg) + sum(new_missed_avg))/all_bbx;
CVIDF1 = (2*mean(prec_avg)*mean(reca_avg))/(mean(prec_avg) + mean(reca_avg));
mean(prec_avg)
mean(reca_avg)
disp(['CVMA--',num2str(CVMOTA)])
disp(['CVIDF1--',num2str(CVIDF1)])

