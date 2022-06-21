% dbstop if all error;
clc;
clear;
% for Res_name = ["Res_app_11/"]
%     Res_mat_path = strcat('./Eval_Data/All_Res_txt/', Res_name);
%     imgPath = '/home/jerome/Documents/CVMOT/dataset/Images/';
%     addpath('/home/jerome/Documents/CVMOT/CVMOT_Tracker6.0');
%     seq = configSeqs;
% 
%     for seq_i = 1 : length(seq)
%         res_sort = cell(1,2);
%         scene_name = seq{seq_i}.name; % 'V1-S_square-G_3';
%         views = {'top','hor'};
%         dataset_directory = '/home/jerome/Documents/CVMOT/dataset';
%         data.cv_images = cell(2,1);
%         data.cv_im_directory = cell(2,1);
%         im_dir = strcat(dataset_directory,'/Images/',scene_name,'/frame_sel/');
%         DIRS = dir(strcat(im_dir, 'hor'));
%         [n, ~] = size(DIRS);
%         maxFrame = n - 2;
%         for frm_i = 1 : maxFrame
%             for view_i = 1:2
%                 resFilename = strcat(Res_mat_path, scene_name,'_', views{view_i},'.txt');
%                 res = dlmread(resFilename);
%                 res = res(:,1:6); 
% %                 res(:,5) = res(:,5) + res(:,3);
% %                 res(:,6) = res(:,6) + res(:,4);
%                 res_i = res(res(:,1) == frm_i,:);    
%                 img_path = strcat(im_dir,views{view_i},'/',sprintf('%04d.jpg',frm_i));
%                 img=imread(img_path);
%                 imshow(img); 
%                 text(round(size(img,2))-100,round(size(img,1))-50,sprintf('#%03d',frm_i),'FontSize',15,'FontWeight','bold','BackgroundColor','yellow');
%                 [numperson, ~] = size(res_i);
%                 for i = 1 : numperson
%                    person = res_i(i, :);
%                    ID = person(2);
%                    bboxWidth =  person(5);
%                    bboxHeight = person(6);
%                    bbox_size = bboxWidth * bboxHeight;
%                    if bbox_size > 0
% %                        line([person(3) person(5) person(5) person(3) person(3)],[person(4) person(4) person(6) person(6) person(4)],'color','g','LineWidth',1.5);
%                        text(round(person(3)) + 30,round(person(4)) - 20,num2str(ID),'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
%                        rectangle('Position',[person(3), person(4), bboxWidth,bboxHeight],'LineWidth',1,'EdgeColor','r');
%                        pause(1);
%                    end
%                    xlim=get(gca,'xlim');
%                    ylim=get(gca,'ylim');
%                    N=100;
%                    text(sum(xlim)/2,sum(ylim)/2+N,'what you want to input, just input here !','horiz','center')
%                end
% %                 im_directory = fullfile(dataset_directory,'Images/',scene_name,'frame_sel/',views{view_i});
% %                 images = dir(fullfile(im_directory,'*.jpg'));
% %                 data.cv_im_directory{view_i} = im_directory;
% %                 data.cv_images{view_i} = images;
%             end
%         end
%         outDir = cell(2,1);
%         for view_i = 1 : 2  
%             outDir{view_i} = sprintf('./output/%s/%s/%s/', Res_name, scene_name, views{view_i});
%         end
%         plotTracking(res_sort, data.cv_im_directory,data.cv_images,1,outDir);
%     end
% 
% end



% Res_mat_path = '.\Eval_Data\GT_txt\';
Res_mat_path = '.\Eval_Data\All_Res_txt\Res_v12\';
imgPath = 'F:\Jerome\CvMHTB_Dataset\Images\';
addpath('F:\Jerome\CvMHT_baseline2.0\');
seq = configSeqs_benchmark;

for seq_i = 1 : length(seq)
    res_sort = cell(1,2);
    scene_name = seq{seq_i}.name; % 'V1-S_square-G_3';
    views = {'t','h2'};
    dataset_directory = 'F:\Jerome\CvMHTB_Dataset\';
    data.cv_images = cell(2,1);
    data.cv_im_directory = cell(2,1);
    for view_i = 1:2
        resFilename = strcat(Res_mat_path, scene_name,'_', views{view_i},'.txt');
        res = dlmread(resFilename);
        res = res(:,1:6); 
        res(:,5) = res(:,5) + res(:,3);
        res(:,6) = res(:,6) + res(:,4);
        id = res(:,2);
        n = unique(id);
        res_sort{view_i} = sortrows(res,2);        
        im_directory = fullfile(dataset_directory,'Images\',scene_name,views{view_i});
        images = dir(fullfile(im_directory,'*.jpg'));
        data.cv_im_directory{view_i} = im_directory;
        data.cv_images{view_i} = images;
    end

    outDir = cell(2,1);
    for view_i = 1 : 2  
        outDir{view_i} = sprintf('./output/%s/%s/', scene_name, views{view_i});
    end
    plotTracking(res_sort, data.cv_im_directory,data.cv_images,1,outDir);
end



 



     
