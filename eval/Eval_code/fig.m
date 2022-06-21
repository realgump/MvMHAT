path = './compareMat/';
tracker = {'DMAN','GMMCP','MDP','Ours'};

plotDrawStyle={
struct('color',[0,0,0],'lineStyle','-'),...%dark red
struct('color',[0,1,0],'lineStyle','--'),...%orange
struct('color',[0,0,1],'lineStyle','-.'),...%Turquoise
struct('color',[1,0,0],'lineStyle',':'),...%purple    %%%%%%%%%%%%%%%%%%%%
};
point = {'^','d','s','p'};

figure(1)
for i = 1 : 4
MAT = cell2mat(struct2cell(load([path,tracker{i},'_hor.mat'])));
xx = 1:1:600;%xx = xx';
y_all = MAT;
x = 1:50:600;
y = MAT(1:50:600);
axis([1 600 0 1]);  
p=polyfit(xx,y_all,3); %����3Ϊ����ݵ�ĸ���-1��pΪ���ض���ʽ��ϵ��
yy=polyval(p,xx); %yyΪ��xx���и�ݶ���ʽ����õ���Ԥ��ֵ���У�
plot(xx,yy,'color',plotDrawStyle{i}.color,'lineStyle', plotDrawStyle{i}.lineStyle,'LineWidth',1.5); hold on;
plot(x,yy(x),'color',plotDrawStyle{i}.color,'Marker',point{i});

end
legend(tracker)

figure(2)
for i = 1 : 4
MAT = cell2mat(struct2cell(load([path,tracker{i},'_top.mat'])));
xx = 1:1:600;
y_all = MAT;
% xx = xx';
x = 1:50:600;
p=polyfit(xx,y_all,3); %����3Ϊ����ݵ�ĸ���-1��pΪ���ض���ʽ��ϵ��
yy=polyval(p,xx); %yyΪ��xx���и�ݶ���ʽ����õ���Ԥ��ֵ���У�
plot(xx,yy,'color',plotDrawStyle{i}.color,'lineStyle', plotDrawStyle{i}.lineStyle,'LineWidth',1.5); hold on;
plot(x,yy(x),'color',plotDrawStyle{i}.color,'Marker',point{i});
axis([0 600 0 1])
end

legend(tracker{1},tracker{2},tracker{3},tracker{4})