 function  [simpleSum]=findIRiniS(seqNum,spacer,spaceMulti,gap1,gap2)
 L=length(seqNum);
 segLen=6;spacer1=ceil(spacer*spaceMulti);intNum=floor(L/segLen);resNum=L-segLen*intNum; %  gap1=1;gap2=3;
% -----------初始查询，完全匹配长度  %-----------spacer:最远间隔距离,----spacer1:初始筛选的最远间隔----%---------总长除以segLen的商，即最多segLen块数 %---------余数
% --------------------segLen倍segLen块覆盖原序列。需分segLen组计算整数块次。为每组较快循环，先计算每组循环块数---------------gap:初始合并时，相同Code可以间隔的碱基数
cyclFir=resNum+1;Num=4^segLen; %NumSum=Num+1;      %--------------分成intNum块的前cyclFir组%-----------有Num个以编码数的结果为排序的组%-----------正反编码的和为4097

headNum1=cell(Num,1);
headNum2=cell(Num,1);%---------------posRem:每个code可能的个数设为最大可能值的1/time---------headNum1正向，headNum2反向
%------------posRem为最后一位，记录位置的个数 %---------记录最初得到的起始位置-----最远相距spacer，数组结构改变
%---储存正向计算结果，segLen：计算整条遍数，Num：结果范围(向量标号1~4096，结果范围0~4095)
%---1：储存正向计算结果，2，储存反向计算结果，cycFir：储存位置，以最多分组记,posRem:非空位置计数；
sumPerCode=zeros(Num,2);  %------------记录每个编码值的每个结果数,1为正向，2为反向
%----------------计算segLen块结果--------
for segLenN=1:segLen                   %------------大循环次数，起始位置逐次差1
%-------------循环次数的确定，前intNum组与后几组差1次    
    if segLenN<=cyclFir
        cyclTime=intNum;
    else
        cyclTime=intNum-1;
    end
    for cyclN=1:cyclTime  %--------------具体进行每块计算---
        startPos=segLenN+(cyclN-1)*segLen; seqPerCycl=seqNum(startPos:startPos+segLen-1);%----------每次计算起始位置,得到计算块-----
      
        codeValue1=scale6(seqPerCycl)+1;  %-------单位为6---正向计算编码数值,+1为与数组下标一致   
        sumPerCode(codeValue1,1)=sumPerCode(codeValue1,1)+1;%-----------------计数每个code值总的个数
        headNum1{codeValue1}(sumPerCode(codeValue1,1))=startPos; %----------记录每次正向编码数值所对应的起始位置
        
        codeValue2=Num-reverScale6(seqPerCycl);        %---反向互补计算---使储存位置与正向编码对应
        sumPerCode(codeValue2,2)=sumPerCode(codeValue2,2)+1;%-----------------计数每个code值总的个数
        headNum2{codeValue2}(sumPerCode(codeValue2,2))=startPos; %------------记录每次反向编码数值所对应的起始位置,与
    end
end    %-----------------headNum赋值完成

% -------------------含两对重复，去掉一对---------------
[numPair,~,~]=xlsread('F:\Matlab\ISR--Matlab\搜索用\正向与反向互补序号对应.xlsx');%----------------------存储正向与反向互补配对编码code
for k=1:Num                    %-----------------Num=4096
    if sumPerCode(k,1)>0 && headNum1{k}(1)>0 && k~=numPair(k)
        headNum1{numPair(k)}(1)=0;
    end
end
% -------------------------------------------------------------------------------------

%---------------------计算非空保留数对------------
nonNull=0;
for ini_N=1:Num
    if  sumPerCode(ini_N,1)>0 && sumPerCode(ini_N,2)>0 && headNum1{ini_N}(1)>0         %--------每个对应的code值对非空，且第一位不为0
        nonNull=nonNull+1;
    end
end

%---------正向与反向互补，两次重复去除---
posSum=cell(nonNull,2);%-----------元胞数组可动态使用,1列正向2列反向%-----------反向位置
nonNulN=0;  %---------------计数非空的code数
%-----------------结果汇总，非空且合并------------
for ini_N=1:Num
       iniSum1=sumPerCode(ini_N,1);iniSum2=sumPerCode(ini_N,2); %-----------记录每个code值里的数据个数       
       if  iniSum1>0 && iniSum2>0 && headNum1{ini_N}(1)>0          %--------每个对应的code值对非空且第一位不为0
           nonNulN=nonNulN+1;
           posSum{nonNulN,1}=zeros(iniSum1,2);posSum{nonNulN,2}=zeros(iniSum2,2);%---------------正反向数据开空间并赋值
           for  ini_Pos=1:iniSum1
                posSum{nonNulN,1}(ini_Pos,1)=headNum1{ini_N}(ini_Pos);
           end
           for  ini_Pos=1:iniSum2
                posSum{nonNulN,2}(ini_Pos,1)=headNum2{ini_N}(ini_Pos);
           end  
       end       
end
%---------正向与反向互补，两次重复去除---
  clear headNum1 headNum2   
  clear sumPerCode   %-------------释放大数组空间
% 
% % %----------------------------------------------------------------------重叠位置连接--------------------------------------------------------------------
    combinSum=cell(nonNull,2);%------------------储存合并后位置,1正向，2反向

 for combiN=1:nonNull
        %----------------------排序并末端位置赋值--------------------------------    
        posSum{combiN,1}(:,1)=sort(posSum{combiN,1}(:,1));
        posSum{combiN,1}(:,2)=posSum{combiN,1}(:,1)+segLen-1;%----------------对正向位置排序，并对各终点赋值
        posSum{combiN,2}(:,1)=sort(posSum{combiN,2}(:,1));
        posSum{combiN,2}(:,2)=posSum{combiN,2}(:,1)+segLen-1;%----------------对反向位置排序，并对各终点赋值
        
%-------------------------邻近位置连接------------------------------------------
    %---------------------------正向相同片段临近连接----------------------
    iH=1;iR=1;   %-----------定位起始位置 
    [rowH,~]=size(posSum{combiN,1});rowHcom=rowH+1;
    [rowR,~]=size(posSum{combiN,2});rowRcom=rowR+1;%-------------------取出列数为下一次作准备
    %---------------以下rowH、rowR纪录非0位置-----------------
    while iH>0 && iH<rowHcom
        jH=iH+1;   %------------先设定要比较的下一位点序号
        while jH<rowHcom && posSum{combiN,1}(jH,1)-posSum{combiN,1}(iH,2)<=gap1     %----------------间隔≤gap1
             posSum{combiN,1}(iH,2)=posSum{combiN,1}(jH,2);   %--------------终止位置变化          
             posSum{combiN,1}(jH,1)=0;            %--------------被连接处置0；
             rowH=rowH-1;    %-------------------置0一次，非0个数-1
             jH=jH+1;
        end
        iH=jH; %----------跳出上一循环时，jH不满足条件，下一循环在jH开始
    end  
    %------------------------反相相同片段临近连接-------------
    while iR>0 && iR<rowRcom
        jR=iR+1;   %------------先设定要比较的下一位点序号
        while jR<rowRcom && posSum{combiN,2}(jR,1)-posSum{combiN,2}(iR,2)<=gap1     %----------------间隔≤gap1
             posSum{combiN,2}(iR,2)=posSum{combiN,2}(jR,2);   %--------------终止位置变化          
             posSum{combiN,2}(jR,1)=0;       %------被连接处置0；
             rowR=rowR-1;%-------------------置0一次，非0个数-1
             jR=jR+1;
        end
        iR=jR; %----------跳出上一循环时，jH不满足条件，下一循环在jH开始
    end
%--------------------赋值给combinSum---------------------------------------
    combinSum{combiN,1}=zeros(rowH,2); combinSum{combiN,2}=zeros(rowR,2);
    ibH2=1;ibR2=1;
    
    for ibH=1:rowHcom-1
        if posSum{combiN,1}(ibH,1)>0
            combinSum{combiN,1}(ibH2,:)=posSum{combiN,1}(ibH,:);
            ibH2=ibH2+1;
        end
    end           %-------------------正向赋值完毕
    for ibR=1:rowRcom-1
        if posSum{combiN,2}(ibR,1)>0
            combinSum{combiN,2}(ibR2,:)=posSum{combiN,2}(ibR,:);
            ibR2=ibR2+1;
        end
    end      %-------------------反向赋值完毕
    
 end
 
  clear posSum
%------------------------------邻近位置连接-------------------------------


%-------------------将相隔远的互补序列第一次筛掉，距离>spacer1------------------   %----------------记录每次循环行列数----------------------
non02=0;       %--------------非0位置初始化为1
[nonNull2,~]=size(combinSum);    remSum=cell(nonNull2,1);
for rowHeadCode=1:nonNull2
    [rowPer1,~]=size(combinSum{rowHeadCode,1});rowPer1Rem=rowPer1;%-----------前向序列位置数
    [rowPer2,~]=size(combinSum{rowHeadCode,2});rowPer2Rem=rowPer2;%-----------反向序列位置数
    checM=zeros(rowPer1,rowPer2);              %-----------距离矩阵

%----------------------------得出距离矩阵,同时处理前向(Head)数据--------------------
    for rowHeadN=1:rowPer1 
        k=0;
        for rowRevN=1:rowPer2
            disL1=combinSum{rowHeadCode,1}(rowHeadN,1)-combinSum{rowHeadCode,2}(rowRevN,1);
            disL2=combinSum{rowHeadCode,1}(rowHeadN,2)-combinSum{rowHeadCode,2}(rowRevN,2);
            disLcomm=combinSum{rowHeadCode,2}(rowRevN,1)-combinSum{rowHeadCode,1}(rowHeadN,2);
            if disLcomm<1
                disLcomm=combinSum{rowHeadCode,1}(rowHeadN,1)-combinSum{rowHeadCode,2}(rowRevN,2);
            end
            
            if disLcomm>0 && disLcomm<=spacer1 || disL1==0 && disL2==0
                checM(rowHeadN,rowRevN)=1;  %----------位置符合的保留
                k=1;
            end
        end    
        if k==0
            combinSum{rowHeadCode,1}(rowHeadN,1)=0;
            rowPer1Rem=rowPer1Rem-1;
        end
    end
%------------------------比较距离，在前向数据存在不归0位点时，对某些反向数据不合格者首位归零------------------
    if  rowPer1Rem>0
        non02=non02+1;            %----------------------存在非0点，code值个数+1
        for rowRevN=1:rowPer2
            if sum(checM(:,rowRevN))==0
                combinSum{rowHeadCode,2}(rowRevN,1)=0;
                rowPer2Rem=rowPer2Rem-1;
            end
        end
    end
    remSum{rowHeadCode,1}=zeros(1,2);
    remSum{rowHeadCode,1}(1,1)=rowPer1Rem;
    remSum{rowHeadCode,1}(1,2)=rowPer2Rem;
end

% %--------------转移赋值---------------------------
disSum=cell(non02,2);disN=0;%--------------disN:disSum计行数
headRow=0;%-----------------为下一步记录前向总行数--------------
for rowHeadCode=1:nonNull2
    [rowPer1,~]=size(combinSum{rowHeadCode,1});%-----------前向序列行数
    [rowPer2,~]=size(combinSum{rowHeadCode,2});%-----------反向序列行数
    rowRem=remSum{rowHeadCode,1}(1,1);
    
    if rowRem>0
        disN=disN+1;                       %----------------------元胞行数变化
        disSum{disN,1}=zeros(rowRem,2);
        disSum{disN,2}=zeros(remSum{rowHeadCode,1}(1,2),2);%----------------对应位置初始化
        row1N=0;row2N=0;
        for rowheadN=1:rowPer1
            if combinSum{rowHeadCode,1}(rowheadN,1)>0
                row1N=row1N+1;%--------------行数变化
                disSum{disN,1}(row1N,:)=combinSum{rowHeadCode,1}(rowheadN,:);

            end 
        end 
        headRow=headRow+row1N;  %---------------记录前向序列总行数
        for rowRevN=1:rowPer2
            if combinSum{rowHeadCode,2}(rowRevN,1)>0
                row2N=row2N+1;%--------------行数变化
                disSum{disN,2}(row2N,:)=combinSum{rowHeadCode,2}(rowRevN,:);
            end
        end
    end
end
clear remSum  combinSum

%-------------正向位置单一化，反向位置可多重----------------
purSum=cell(headRow,2);%-----------下一步展开用
orderRem=zeros(headRow,3);%---------------用于排序-------
headRowN=0;  %-----------------记录行数
for rowHeadCode=1:disN
    [rowPer1,~]=size(disSum{rowHeadCode,1});
%--------------------前向位置赋值------------------------    
        if rowPer1==1
            headRowN=headRowN+1;
            purSum{headRowN,1}=disSum{rowHeadCode,1};%-------------------------正向赋值
            orderRem(headRowN,1)=purSum{headRowN,1}(1,1);%---------------第一列记录前向序列起始位置
            purSum{headRowN,2}=disSum{rowHeadCode,2};%-------------------------反向位置同时赋值  
        else
            for row1N=1:rowPer1
               headRowN=headRowN+1; 
               purSum{headRowN,1}=disSum{rowHeadCode,1}(row1N,:);%------------------正向赋值
               
               orderRem(headRowN,1)=purSum{headRowN,1}(1,1);%---------------第一列记录前向序列起始位置
               
               %-----------------反向赋值------------------------
               [rowPer2,~]=size(disSum{rowHeadCode,2}); 
               row2PurN=0;
               for row2N=1:rowPer2
                    disL1=disSum{rowHeadCode,1}(row1N,1)-disSum{rowHeadCode,2}(row2N,1);        %--------------包含两端不重复、全段重复序列
                    disL2=disSum{rowHeadCode,1}(row1N,2)-disSum{rowHeadCode,2}(row2N,2);
                    disLcomm=disSum{rowHeadCode,2}(row2N,1)-disSum{rowHeadCode,1}(row1N,2);
                    if disLcomm<1
                       disLcomm=disSum{rowHeadCode,1}(row1N,1)-disSum{rowHeadCode,2}(row2N,2);
                    end
                    if disLcomm>0 && disLcomm<=spacer1 || disL1==0 && disL2==0       %----------以正向位置为准，位置符合的保留
                       row2PurN=row2PurN+1;
                       purSum{headRowN,2}(row2PurN,:)=disSum{rowHeadCode,2}(row2N,:);
                   end
               end
            end
        end    
end
 clear disSum
 orderRem(:,2)=sort(orderRem(:,1));%------------------第二列为将前向位置排序----------------
for headRowN2=1:headRow          %-------------第一列
    for headRowN=1:headRow      %--------------第二列
        if orderRem(headRowN,1)==orderRem(headRowN2,2)
           orderRem(headRowN2,3)=headRowN;      %-----------------------第三列记录排序前位置
           break;
        end  
    end
end


%--------------------------------按顺序赋值------------------------
orderSum=cell(headRow,2);   %---------------储存排序后各个位置
for headRowN=1:headRow
    orderSum{headRowN,1}=purSum{orderRem(headRowN,3),1};
    orderSum{headRowN,2}=purSum{orderRem(headRowN,3),2};
end
clear purSum orderRem

%----------------------反向位置单一化，正向与反向一一对应--------------------
for headRowN=1:headRow
   [rowPer2,~]=size(orderSum{headRowN,2});
   if rowPer2>1
      for rowPer1N=2:rowPer2
          orderSum{headRowN,1}(rowPer1N,:)=orderSum{headRowN,1}(1,:);
      end
   end
end
%------------------------------------------------------------------

%---------------------------将Code值不同，但是位置邻近的片段连接：需满足正向、反向位置都相近-----------------------
for headRowN=1:headRow-1
    if sum(orderSum{headRowN,2}(:,1))>0           %------------两两对比，第一组数据需存在非零点
       [rowPer1,~]=size(orderSum{headRowN,2});%----------------取出第一组行数用于操作       
        for  headRowN2=headRowN+1:headRow
            if sum(orderSum{headRowN2,2}(:,1))>0  %--------------第二组数据存在非零点
               [rowPer2,~]=size(orderSum{headRowN2,2});%---------------取出第二组行数用于操作
                if orderSum{headRowN2,1}(1,1)<=orderSum{headRowN,1}(1,2)+gap2 &&  orderSum{headRowN2,1}(1,2)>=orderSum{headRowN,1}(1,2) %------------下一个位置需满足的条件
                    for rowPer1N=rowPer1:-1:1      %---------------从最后一个开始比较.------同一正向位置的，后来者与前一个比较--------
                        if orderSum{headRowN,2}(rowPer1N,1)>0
                            %-------------------同一正向位置相邻两个先比较-----------------------                            
                             if rowPer1N<rowPer1 %-------------从倒数第二个开始，首先与同Code值的下一行位置比较  
                                 if orderSum{headRowN,2}(rowPer1N+1,1)>0             
                                     if orderSum{headRowN,2}(rowPer1N+1,1)<=orderSum{headRowN,2}(rowPer1N,1)
                                         orderSum{headRowN,1}(rowPer1N,2)=orderSum{headRowN,1}(rowPer1N+1,2);%-----------同一Code值，前向位置的后位会有变化
                                         orderSum{headRowN,2}(rowPer1N,:)=orderSum{headRowN,2}(rowPer1N+1,:);%------------大包小的情况，反向位置变化-------------------
                                         orderSum{headRowN,2}(rowPer1N+1,1)=0;  %-------------被合并位置变为0-----------
                                     else %--------------相当于多加一个条件-------
                                         if  orderSum{headRowN,2}(rowPer1N+1,1)<=orderSum{headRowN,2}(rowPer1N,2)+gap2
                                              orderSum{headRowN,1}(rowPer1N,2)=orderSum{headRowN,1}(rowPer1N+1,2);%-----------同一Code值，前向位置的后位会有变化
                                              orderSum{headRowN,2}(rowPer1N,2)=orderSum{headRowN,2}(rowPer1N+1,2);%------------连续位置情况--反向位置变化-------------------
                                              orderSum{headRowN,2}(rowPer1N+1,1)=0;  %-------------被合并位置变为0-----------
                                         end
                                     end
                                 end
                             end
                             %-------------------------------------------------------------------------
                             
                             %--------------------不同正向位置再比较----------------
                             for rowPer2N=rowPer2:-1:1            %-----------------只考虑连续或交叉的情况------
                                 if orderSum{headRowN2,2}(rowPer2N,1)>0    
                                     if orderSum{headRowN2,2}(rowPer2N,1)<orderSum{headRowN,2}(rowPer1N,1) && orderSum{headRowN2,2}(rowPer2N,2)>=orderSum{headRowN,2}(rowPer1N,1)-gap2 && orderSum{headRowN2,2}(rowPer2N,2)<=orderSum{headRowN,2}(rowPer1N,2)
                                         orderSum{headRowN,1}(rowPer1N,2)=orderSum{headRowN2,1}(rowPer2N,2);%--------------正向位置变化
                                         orderSum{headRowN,2}(rowPer1N,1)=orderSum{headRowN2,2}(rowPer2N,1);%--------------反向位置变化
                                         orderSum{headRowN2,2}(rowPer2N,1)=0;%---------------被合并位置变为0
                                     end
                                 end
                             end
                        end
                    end                    
                else
                    break;             
                end
       
            end  
        end
    end
end

%-------------------对连接后的数组非0点计数---------------------
headRowL=0;%----------用于不同正向位置计数
for headRowN=1:headRow
    if sum(orderSum{headRowN,2}(:,1))>0
        headRowL=headRowL+1;
    end
end

simpleSum=cell(headRowL,2);%------------储存简并非0位置------------

headRowLN=0;
for headRowN=1:headRow
    if sum(orderSum{headRowN,2}(:,1))>0
        headRowLN=headRowLN+1;      %----------------记录正向位置数目的变化
        simpRow=0;        %--------------记录每一位置的行数变化
        [rowPer1,~]=size(orderSum{headRowN,2});
        for rowPer1N=1:rowPer1
            if orderSum{headRowN,2}(rowPer1N,1)>0
                simpRow=simpRow+1;    
                simpleSum{headRowLN,1}(simpRow,:)=orderSum{headRowN,1}(rowPer1N,:);
                simpleSum{headRowLN,2}(simpRow,:)=orderSum{headRowN,2}(rowPer1N,:);
            end
        end
    end
end

clear  orderSum
