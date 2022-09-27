
function [Best_pos,Best_score,curve]=SSANew(pop,Max_iter,lb,ub,dim,fobj)

ST = 0.8;
PD = 0.7;
SD = 0.2;

PDNumber = pop*PD; 
SDNumber = pop - pop*PD;
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%利用Tent映射策略种群初始化
X0=initializationNew(pop,dim,ub,lb);
X = X0;
%计算初始适应度值
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end
 [fitness, index]= sort(fitness);%排序
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%全局最优适应度值
for i = 1:pop
    X(i,:) = X0(index(i),:);
end
curve=zeros(1,Max_iter);
GBestX = X(1,:);%全局最优位置
X_new = X;
for i = 1: Max_iter
    
    BestF = fitness(1);
    WorstF = fitness(end);

    
    R2 = rand(1);
   for j = 1:PDNumber
      if(R2<ST)
          X_new(j,:) = X(j,:).*exp(-j/(rand(1)*Max_iter));
      else
          X_new(j,:) = X(j,:) + randn()*ones(1,dim);
      end     
   end
   for j = PDNumber+1:pop
%        if(j>(pop/2))
        if(j>(pop - PDNumber)/2 + PDNumber)
          X_new(j,:)= randn().*exp((X(end,:) - X(j,:))/j^2);
       else
          %产生-1，1的随机数
          A = ones(1,dim);
          for a = 1:dim
            if(rand()>0.5)
                A(a) = -1;
            end
          end 
          AA = A'*inv(A*A');     
          X_new(j,:)= X(1,:) + abs(X(j,:) - X(1,:)).*AA';
       end
   end
   Temp = randperm(pop);
   SDchooseIndex = Temp(1:SDNumber); 
   for j = 1:SDNumber
       if(fitness(SDchooseIndex(j))>BestF)
           X_new(SDchooseIndex(j),:) = X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - X(1,:));
       elseif(fitness(SDchooseIndex(j))== BestF)
           K = 2*rand() -1;
           X_new(SDchooseIndex(j),:) = X(SDchooseIndex(j),:) + K.*(abs( X(SDchooseIndex(j),:) - X(end,:))./(fitness(SDchooseIndex(j)) - fitness(end) + 10^-8));
       end
   end
   %边界控制
   for j = 1:pop
       for a = 1: dim
           if(X_new(j,a)>ub)
               X_new(j,a) =ub(a);
           end
           if(X_new(j,a)<lb)
               X_new(j,a) =lb(a);
           end
       end
   end 
   %更新位置
   for j=1:pop
    fitness_new(j) = fobj(X_new(j,:));
   end
   %% Step6 高斯变异
   %计算平均适应度值
   avgF = mean(fitness_new);
   for j = 1:pop
       if(fitness_new(j) < avgF)
           Temp = X(j,:).*(1 + randn());
           %边界控制
           Temp(Temp>ub) = ub(Temp>ub);
           Temp(Temp<lb) = lb(Temp<lb);
           ftemp = fobj(Temp);
           if(ftemp<fitness_new(j))
               fitness_new(j) = ftemp;
               X(j,:) = Temp;
           end
       else
           TentValue =  Tent(dim);%tent 扰动
           Temp = X(j,:).*(1 + TentValue);
           ftemp = fobj(Temp);
           %边界控制
           Temp(Temp>ub) = ub(Temp>ub);
           Temp(Temp<lb) = lb(Temp<lb);
           if(ftemp<fitness_new(j))
               fitness_new(j) = ftemp;
               X(j,:) = Temp;
           end
       end
   end
       
   X = X_new;
   fitness = fitness_new;
   [fitness, index]= sort(fitness);
   for j = 1:pop
      X(j,:) = X(index(j),:);
   end  
    if(fitness(1) < GBestF)
       GBestF = fitness(1);
        GBestX = X(1,:);   
    end 
   BestF = fitness(1);
   WorstF = fitness(end);
   
   
   curve(i) = GBestF;
end
Best_pos =GBestX;
Best_score = curve(end);
end

%Copyright (c) 2020, JackXu

