function [inter_top_hor,top_id_trans,hor_id_trans] = transform_top_hor(top_id,hor_id)

top_id_trans = zeros(1,length(top_id));
hor_id_trans = zeros(1,length(hor_id));

[inter_top_hor,inter_top_idx,inter_hor_idx] = intersect(top_id,hor_id);

top_id_trans(inter_top_idx) = top_id(inter_top_idx);
hor_id_trans(inter_hor_idx) = hor_id(inter_hor_idx);

end