function  top_match = top_match_box(top_res_i,top_gt_i,thr_dis)

dis_thr_top = thr_dis;
objs_top_GT = top_gt_i(:,3:6);
objs_top = top_res_i(:,3:6);

top_dis = zeros(size(objs_top_GT,1),size(objs_top,1));
top_match = zeros(size(objs_top_GT,1),size(objs_top,1));


    for i = 1 : size(objs_top_GT,1)
        for j = 1 : size(objs_top,1)  
            top_dis(i,j) = get_center_dis(objs_top_GT(i,:),objs_top(j,:)); 
        end
    end
    if size(top_dis,2)~=0  
        [min_row,col] = min(top_dis,[],2);    
        for i = 1: size(top_dis,1)
            if min_row(i) < dis_thr_top
                top_match(i,col(i)) = 1;
            end
        end   
        [min_col,row] = min(top_dis,[],1);   
        for i = 1: size(top_dis,2)
            if min_col(i) < dis_thr_top
                top_match(row(i),i) = 1;
            end
        end
        %----------------------------------------------------------------------
%         det2gt_top = zeros(size(objs_top_GT,1),1);
%         [gt,det] = find(top_match == 1);
%         det2gt_top(gt) =  lab_top_final(det);   
%         lab_top_end = det2gt_top';   
    else
%         lab_top_end = 0;      
    end
    
end