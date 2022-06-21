
filepath = ['F:\Jerome\Evaluation\Eval_CVMOT\Eval_Data\GT_mat\'];

addpath('F:\Jerome\CvMHT_baseline2.0\');
seq = configSeqs_benchmark;

for seq_i = 1 : length(seq)
    
change = 0;    
scene_name = seq{seq_i}.name;
num_hor_seq =  seq{seq_i}.num_hor; 

views = {'t','h1','h2', 'h3','h4'};

for view_i = 1 : num_hor_seq + 1
    change = 0;    
    if strcmp(scene_name,'V3_G1') && view_i == 1
        change = 1;
    end

    view = views{view_i};
    
    gt = struct2cell(load([filepath,scene_name,'_',view,'.mat']));
    gt = gt{1};
    
    fid = fopen(['../GT_txt/',scene_name,'_',view,'.txt'],'w');
  
    new_gt_all = cell(length(gt),1);
    for frame_i = 1 : length(gt)   
        if (seq_i == 26 && frame_i >= 1501) || (seq_i == 29 && frame_i >= 1301) || (seq_i == 30 && frame_i >= 1201)
            break
        end
        gt_frm = cell2mat(gt(frame_i,2)); % x,y,x2,y2,ID
        if isempty(gt_frm)
           continue; 
        end
        new_gt_frame = zeros(size(gt_frm,1),6);
        new_gt_frame(:,1) = frame_i; % frm
        new_gt_frame(:,2) = gt_frm(:,5); % ID
        if change ==1
            w = (gt_frm(:,3) - gt_frm(:,1));
            h = (gt_frm(:,4) - gt_frm(:,2));
            new_gt_frame(:,3:4) = floor(gt_frm(:,1:2) - 1/6 * [w,h]); %x,y       
            new_gt_frame(:,5) = floor(4/3 * w); %w
            new_gt_frame(:,6) = floor(4/3 * h); %h
            new_gt_all{frame_i} = new_gt_frame;
        else
            new_gt_frame(:,3:4) = gt_frm(:,1:2); %x,y     
            new_gt_frame(:,5) =  (gt_frm(:,3) - gt_frm(:,1)); %w
            new_gt_frame(:,6) =  (gt_frm(:,4) - gt_frm(:,2)); %h
            new_gt_all{frame_i} = new_gt_frame;
            
        end
    end   
    
    new_gt_mat = cell2mat(new_gt_all);
    new_gt_mat_sort = sortrows(new_gt_mat,2);

    for i = 1 : size(new_gt_mat_sort,1)
        
    new_gt_mat_sort_i = new_gt_mat_sort(i,:);

    fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',[new_gt_mat_sort_i,-1,-1,-1,-1]);

    end
        

    fclose(fid);
    
end

end
