
clear all 
clc

SearchAgents_no=50; % Number of search agents ��Ⱥ����

Function_name='F8'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper) �趨��Ӧ�Ⱥ���

Max_iteration=1000;
% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);  %�趨�߽��Լ��Ż�����

[Best_pos,Best_score,SSA_curve]=SSANew(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %��ʼ�Ż�

figure('Position',[269   240   660   290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
plot(SSA_curve,'Color','r')
hold on;
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('������ȸ�Ż��㷨SSA')
% 
display(['The best solution obtained by ������ȸ�Ż��㷨 is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by ������ȸ�Ż��㷨 is : ', num2str(Best_score)]);
   



%Copyright (c) 2020, JackXu


