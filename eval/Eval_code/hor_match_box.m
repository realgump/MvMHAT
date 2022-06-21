function  hor_match = hor_match_box(hor_res_i,hor_gt_i,thr_IOU)

IOU_thr = thr_IOU;
objs_ego_GT = hor_gt_i(:,3:6);
objs_ego = hor_res_i(:,3:6);

ego_IOU = zeros(size(objs_ego_GT,1),size(objs_ego,1));
ego_match = zeros(size(objs_ego_GT,1),size(objs_ego,1));
    
    for m = 1 : size(objs_ego_GT,1)
        for n = 1 : size(objs_ego,1)            
         ego_IOU(m,n) = get_IOU(objs_ego_GT(m,:),objs_ego(n,:));   
        end
    end
    
    if size(ego_IOU,2)~=0  
    [max_row,col] = max(ego_IOU,[],2);    
    for i = 1: size(ego_IOU,1)
        try
        if max_row(i) > IOU_thr
            ego_match(i,col(i)) = 1;
        end
        catch err
        end
    end   
    [max_col,row] = max(ego_IOU,[],1);  
    for i = 1: size(ego_IOU,2)
        if max_col(i) > IOU_thr
            ego_match(row(i),i) = 1;
        end
    end
    
    end
%     det2gt = zeros(size(objs_ego_GT,1),1);
%     [gt,det] = find(ego_match == 1);
%     det2gt(gt) =  lab_ego_final(det);   
%     lab_ego_end = det2gt';   
    hor_match = ego_match;
end