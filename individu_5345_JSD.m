clear;
clc;
close all;
fpi=fopen('control_pro_cancer.txt');           %对照样本
hline1 = textscan(fpi, '%s', 1, 'delimiter', '\n');  
field1=textscan(hline1{1}{1},'%s'); 
format='%s';      
% format=[format,' %s'];
for i=2:5           
 format=[format,' %f'];    
end
plines =textscan(fpi, format,1000000,'delimiter', '\t');  
pipi=plines{1};  
pprofile = [];   
for i = 2:5         
pprofile = [pprofile, plines{i}];  
end   
fclose(fpi);

fpi=fopen('case_pro_cancer.txt');                %病例样本
hline2 = textscan(fpi, '%s', 1, 'delimiter', '\n');
field2=textscan(hline2{1}{1},'%s');
clear format;
format='%s';   
% format=[format,' %s'];
for i=2:31
 format=[format,' %f'];   
end
mlines =textscan(fpi, format,1000000,'delimiter', '\t');
mipi=mlines{1};
mprofile = [];
for i = 2 :31
 mprofile = [mprofile, mlines{i}];
end
fclose(fpi);

%样本分组
psize=size(pprofile);
tempcontrol=pprofile;
reference_num=4;
stage_num=7;
patients_num=[4,6,4,4,4,4,4];
tempcase=zeros(psize(1),7,7);
tempcase(:,1,1:patients_num(1))=mprofile(:,1:4);    % Stage IA    
tempcase(:,2,1:patients_num(2))=mprofile(:,5:10);   % Stage IB
tempcase(:,3,1:patients_num(3))=mprofile(:,11:14);  % Stage IIA
tempcase(:,4,1:patients_num(4))=mprofile(:,15:18); % Stage IIB
tempcase(:,5,1:patients_num(5))=mprofile(:,19:22); % Stage IIIA
tempcase(:,6,1:patients_num(6))=mprofile(:,23:26); % Stage IIIB
tempcase(:,7,1:patients_num(7))=mprofile(:,27:30); % Stage IV
msize=size(tempcase);

%为参考样本拟合guass分布
jsd_t=zeros(1,stage_num);
jsd1_t=zeros(psize(1),stage_num);
jsd2_t=zeros(stage_num,7);

for t=1:stage_num
        for s=1:patients_num(t)  
            for i=1:psize(1)  
             mu=mean(tempcontrol(i,1:reference_num));   
             sigma=std(tempcontrol(i,1:reference_num)); 
             procancer_pdf(i,t,s)= normpdf(reshape(tempcase(i,t,s),1,1),mu,sigma); %case sample 
             normal_pdf(i,:)=normpdf(tempcontrol(i,1:reference_num),mu,sigma); %control sample 
             procancer_cdf(i,t,s)=normcdf(reshape(tempcase(i,t,s),1,1),mu,sigma);  
             if procancer_cdf(i,t,s)<(normcdf(mu-4*sigma,mu,sigma)) 
                procancer_cdf(i,t,s)=1-procancer_cdf(i,t,s);
             end
             normal_cdf(i,1:reference_num)=normcdf(tempcontrol(i,1:reference_num),mu,sigma);
            end
            [tmp_com_idx,index]=sort(procancer_cdf(:,t,s),1,'descend');
         for i=1:psize(1) 
             geneid_cdf_num(i,t,s)=pipi(index(i));
             procancer_cdf_num(i,t,s)=procancer_cdf(index(i),t,s);
             normal_cdf_num(i,t,s)=mean(normal_cdf(index(i),1:reference_num));
         end
        end
end
 
for t=1:stage_num
  %构造参考分布和扰动分布  
             normal_cumulative_area=normal_cdf_num(:,t,:);
             P=normal_cumulative_area;%P=P/sum(P);
             procancer_cumulative_area=procancer_cdf_num(:,t,:);
             Q=procancer_cumulative_area; %Q=Q/sum(Q);  
           
  %检测临界状态  
         for s=1:patients_num(t)
          for i=1:psize(1) 
           jsd(i,t,s)=0.5*(sum(P(i,1,s).*log(2.*P(i,1,s)./(P(i,1,s)+Q(i,1,s))))+sum(Q(i,1,s).*log(2.*Q(i,1,s)./(P(i,1,s)+Q(i,1,s))))); %JS散度   
           jsd1_t(i,t)=jsd1_t(i,t)+jsd(i,t,s);
           jsd2_t(t,s)=jsd2_t(t,s)+jsd(i,t,s)
           normal_cumulative_area1(i,s)=normal_cdf_num(i,t,s);
           P1(i)=mean(normal_cumulative_area1(i,:));
          end
         end
           jsd_t(t)=jsd_t(t)+sum(jsd1_t(:,t))/msize(1);
           [jsdgene_idx,indx]=sort(jsd1_t(:,t),'descend');
 %识别关键标志物基因
           for k=1:100
           jsd_genstage_id(k,t)=geneid_cdf_num(indx(k),t);
           jsd_top5gene_data(k,t)=jsd1_t(indx(k),t);
        end
end
%临界样本的关键标志物提取
 for k=1:1000
         for s=1:patients_num(t)
       jsd_gensamp_id(k,s)=geneid_cdf_num(indx(k),2,s);
         end
 end

figure(1);
t=[1:7];
plot(t,jsd_t,'r','LineWidth',3);
set(gca,'XTick',1:7);
B={ '0h' '6h' '9h' '12h' '18h' '24h' '48h'};
set(gca,'XTickLabel',B);
xlabel('Time');
ylabel('ICI score');
title('The average ICI score for Prostate Cancer ');

figure(2);
axes1 = axes('Parent',figure(2));
hold(axes1,'on');
mesh(jsd1_t,'LineWidth',3);
set(gca,'XTick',1:7);
B={ '0h' '6h' '9h' '12h' '18h' '24h' '48h'};
set(gca,'XTickLabel',B);
xlabel('Time');
ylabel('Gene');
zlabel('ICI score');
title('The global ICI score for Prostate Cancer ');

figure(3);
axes1 = axes('Parent',figure(3));
hold(axes1,'on');
mesh(jsd_top5gene_data,'LineWidth',3);
set(gca,'XTick',1:7);
B={ '0h' '6h' '9h' '12h' '18h' '24h' '48h'};
set(gca,'XTickLabel',B);
xlabel('Time');
ylabel('Gene');
zlabel('ICI score');
title('The ICI score of sJSD signal markers for Prostate Cancer ');
