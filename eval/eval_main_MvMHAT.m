%dbstop if all error;
clear;
clc;
seqs = configSeqs_benchmark;
all_seq_len = 0;

%% cal allsequence num
for seq_i = 1 :  length(seqs)
    all_seq_len = all_seq_len +  seqs{seq_i}.num_hor + 1; 
end
top_seq_len = length(seqs);
hor_seq_len = all_seq_len - top_seq_len;

horSequences = cell(1,hor_seq_len);

views = {'1','2', '3','4'};
tmp = 1;
hor_tmp = 1;
for seq_i = 1 :  length(seqs)
    num_hor_seq =  seqs{seq_i}.num_hor; 
    for view_i = 1 : num_hor_seq 
            horSequences{hor_tmp}.name =  strcat(seqs{seq_i}.name,'_',views{view_i});
            horSequences{hor_tmp}.start =  strcat(seqs{seq_i}.startFrame);
            horSequences{hor_tmp}.end =  strcat(seqs{seq_i}.endFrame);
            hor_tmp = hor_tmp + 1;
        tmp = tmp + 1;
    end
end

Sequences = horSequences;


Res_path = './Eval_Data/Res/MvMHAT/';
GT_path = './Eval_Data/mot_gt/';
allMets = evaluation(Sequences,Res_path,GT_path,'MOT15');
MHA_eval = evaluation_MHA(seqs,Res_path,GT_path);

% MHA_eval = evaluation_MHA_plot(Sequences,Res_path,GT_path);
