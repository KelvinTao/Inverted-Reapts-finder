%-----------------------寻找反向重复序列，满足:至少有segLen个连续正确配对的碱基;动态规划后，正确配对的碱基比率>=prob------------------
%------panlinSum：回文位置--------simpleSum3：反向重复序列位置---------
function [simpleSum3,panlinSum]=IRFinderS(seqNum,segLen,spacer,spaceMulti,prob,gap1,gap2,jumpIni)

[simpleSum]=findIRiniS(seqNum,spacer,spaceMulti,gap1,gap2);%-------初步结果-------------
 [RowHead,~]=size(simpleSum);
 seqLen=length(seqNum);%------------循环时作界限-------
 
 %-----------------------先找出长度<2*segLen的回文序列-------------
 PanlinN=0;
 for rowHeadN=1:RowHead            %---------从前向后-----------
     [rowRev,~]=size(simpleSum{rowHeadN,2});
     for rowRevN=1:rowRev
         if  simpleSum{rowHeadN,1}(rowRevN,1)==simpleSum{rowHeadN,2}(rowRevN,1)&&simpleSum{rowHeadN,1}(rowRevN,2)==simpleSum{rowHeadN,2}(rowRevN,2) && simpleSum{rowHeadN,1}(rowRevN,2)-simpleSum{rowHeadN,1}(rowRevN,1)<11
             PanlinN=PanlinN+1;
         end
     end
 end
 panlinSum=cell(PanlinN,1);PanlinN=0;
  for rowHeadN=1:RowHead            %---------从前向后-----------
     [rowRev,~]=size(simpleSum{rowHeadN,2});
     for rowRevN=1:rowRev
         if  simpleSum{rowHeadN,1}(rowRevN,1)==simpleSum{rowHeadN,2}(rowRevN,1)&&simpleSum{rowHeadN,1}(rowRevN,2)==simpleSum{rowHeadN,2}(rowRevN,2) && simpleSum{rowHeadN,1}(rowRevN,2)-simpleSum{rowHeadN,1}(rowRevN,1)<11
              PanlinN=PanlinN+1;
              panlinSum{PanlinN,1}=simpleSum{rowHeadN,1}(rowRevN,:);
              simpleSum{rowHeadN,2}(rowRevN,1)=0;%-----------------回文序列不必做处理
         end
     end
 end
%----------------------------------------------------------------------------------- 
 
 
 for rowHeadN=1:RowHead%RowHead:-1:1
    headPos=simpleSum{rowHeadN,1}; revPos=simpleSum{rowHeadN,2};%--------短暂记录位置
    [rowRev,~]=size(revPos);
    for rowRevN=1:rowRev%rowRev:-1:1
         if revPos(rowRevN,1)>0
             baseJump=10;probMatc=1;%------------baseJump：每步延伸的碱基数,第一次为初始片段--
            %------------------------动态规划全过程---------------------------------
           %-------------第一步baseJump=0------------
             headSmaIf=1;   %---------------前小后大为1，前大后小为0---------------- 
            if revPos(rowRevN,2)>headPos(rowRevN,1)            %-------------大小关系确定，正向A序列在前-----反向B序列在后----------
                aS0=headPos(rowRevN,1);aE0=headPos(rowRevN,2);bS0=revPos(rowRevN,2);bE0=revPos(rowRevN,1);%------------aS0<aE0----------bS0>bE0----------
                
                meanPos=(aE0+bS0)/2;   %-------------计算正反序列中间位置
                AmidEnd=floor(meanPos); BmidEnd=AmidEnd+1;      %-------------计算初始正反序列最远可达的中间位置                
                
                aS1=aS0;aE1=aE0;bS1=bS0;bE1=bE0; %-----------baseJump=0,时用

                if aE1>AmidEnd
                    aE1=AmidEnd;
                end             %------------------正向片段末端位置检测-----------   
                if bE1<BmidEnd
                   bE1=BmidEnd;
                end           %-------------------反向片段前段位置检测-----------
                
            else
                headSmaIf=0;
                bS0=headPos(rowRevN,2);bE0=headPos(rowRevN,1);aS0=revPos(rowRevN,1);aE0=revPos(rowRevN,2);
                
                meanPos=(aE0+bS0)/2;   %-------------计算正反序列中间位置
                AmidEnd=floor(meanPos); BmidEnd=AmidEnd+1;      %-------------计算初始正反序列最远可达的中间位置      
                
                aS1=aS0;aE1=aE0;bS1=bS0;bE1=bE0;%-----------baseJump=0,时用
                
                if aE1>AmidEnd
                    aE1=AmidEnd;
                end             %------------------正向片段位置检测-----------   
                if bE1<BmidEnd
                   bE1=BmidEnd;
                end           %-------------------反向片段位置检测-----------                 
            end
            
            A=seqNum(aS1:aE1); B=seqNum(bS1:-1:bE1);
            [probMatc1,~,~]=IRCodlocDyMatc(A,B);%---------------第一次动态匹配
         if probMatc1>prob 
            k=0;%----------------probMatc记录匹配概率，其他记录最大间隔----%-----------------aS0,aE0;bS0,bE0为满足条件的位置---
            while  1               %----------------永真循环，下面有终止语句-----------
            %---------------确定起始比较位点-----------------------------------
                aS=aS0-baseJump;aE=aE0+baseJump; 
                if aS<1
                    aS=1;
                end
                if aE>AmidEnd
                    aE=AmidEnd;
                end             %------------------正向片段位置检测-----------
                bS=bS0+baseJump;bE=bE0-baseJump;  
                if bS>seqLen
                   bS=seqLen;
                end
                if bE<BmidEnd
                   bE=BmidEnd;
                end           %-------------------反向片段位置检测-----------
            %----------------------------------------------------------------------
                A=seqNum(aS:aE); B=seqNum(bS:-1:bE);        %------------截取正反向序列 
                [probMatc,Aend,Bend]=IRCodlocDyMatc(A,B);
                k=k+1;
                if probMatc>=prob 
                        aSN=aS-1+Aend(1);aEN=aS-1+Aend(2);bSN=bS+1-Bend(1);bEN=bS+1-Bend(2);%-----------------aS0,aE0;bS0,bE0为满足条件的位置
                     if  aSN==aS0&&aEN==aE0||bSN==bS0&&bEN==bE0 %----有一片段长度不增加,停止循环----
                            if ~(aSN==aS0&&aEN==aE0)
                                aS0=aSN;aE0=aEN;
                            end
                            if ~(bSN==bS0&&bEN==bE0)
                                bS0=bSN;bE0=bEN;
                            end
                            break;%-------------在匹配值较大的情况下，再进行判断                           
                     else          %-----------------长度都增加，继续循环-----------
                           aS0=aSN;aE0=aEN;bS0=bSN;bE0=bEN;
                           baseJump=baseJump+2;%------------步幅太小时，加大----
                           if k>=30
                               break;
                           end       %--------------------循环次数过多，停止循环-----------                           
                     end   
                else       %---------------匹配值较小
                        if baseJump>jumpIni
                            baseJump=jumpIni;
                        else
                            baseJump=baseJump-2;  %-----------步幅太大，缩小----最好最后能减为0-----
                        end
                end     
            end
            %-------------------------------------------------------------------------------------------------------------------------------
        else 
           probMatc=probMatc1;
           aS0=aS1;aE0=aE1;bS0=bS1;bE0=bE1;
        end  
       %-----------------------------储存结果------------------------------------------------
          if headSmaIf==1        %--------------前小后大正常赋值
             simpleSum{rowHeadN,1}(rowRevN,1)=aS0;simpleSum{rowHeadN,1}(rowRevN,2)=aE0;%-----------正向片段赋值----------
             simpleSum{rowHeadN,2}(rowRevN,1)=bE0;simpleSum{rowHeadN,2}(rowRevN,2)=bS0;%-----------反向片段赋值----------
          else            
             simpleSum{rowHeadN,1}(rowRevN,1)=bE0;simpleSum{rowHeadN,1}(rowRevN,2)=bS0;%-----------正向片段赋值----------
             simpleSum{rowHeadN,2}(rowRevN,1)=aS0;simpleSum{rowHeadN,2}(rowRevN,2)=aE0;%-----------反向片段赋值----------              
          end
             simpleSum{rowHeadN,3}(rowRevN)=probMatc;%-------------第三列元胞保存匹配概率；
       %-------------------------------------------------------------------------------------------- 
         end
     %----边做边排除------------------    
         
         
         
         
         
    end
end
%-------------------------最后结果进一步合并------------------------
for rowHeadN=1:RowHead            %---------从前向后-----------
    [rowRev,~]=size(simpleSum{rowHeadN,2});
    for rowRevN=1:rowRev
        if simpleSum{rowHeadN,2}(rowRevN,1)>0;
             rowCom2=rowHeadN+10;
             rowCom1=rowHeadN-10;
          if rowCom2>RowHead
             rowCom2=RowHead;
          end
          if rowCom1<1
             rowCom1=1;
          end
          %----------------------比较位置的选取----------------
          for i=rowCom1:rowCom2
              [rowRev2,~]=size(simpleSum{i,2});
              for j=1:rowRev2
                 if ~(i==rowHeadN && j==rowRevN)      %------------合并时除去本位置
                     if simpleSum{i,2}(j,1)>0
                         if simpleSum{i,1}(j,1)>=simpleSum{rowHeadN,1}(rowRevN,1) && simpleSum{i,1}(j,2)<=simpleSum{rowHeadN,1}(rowRevN,2) 
                             if simpleSum{i,2}(j,1)>=simpleSum{rowHeadN,2}(rowRevN,1) && simpleSum{i,2}(j,2)<=simpleSum{rowHeadN,2}(rowRevN,2)
                                simpleSum{i,2}(j,1)=0;
                             end
                         end
                     end
                 end
              end
           end
        end
    end
end
%-----------------------------------------------------

%--------------筛选出配对碱基数大于等于segLen的片段间隔≤spacer的片段---------
for rowHeadN=1:RowHead
    [rowRev,~]=size(simpleSum{rowHeadN,2});
    for rowRevN=1:rowRev
        if simpleSum{rowHeadN,2}(rowRevN,1)>0;
            if abs(simpleSum{rowHeadN,1}(rowRevN,1)-simpleSum{rowHeadN,1}(rowRevN,2))+1<segLen || abs(simpleSum{rowHeadN,2}(rowRevN,1)-simpleSum{rowHeadN,1}(rowRevN,2))+1<segLen
                simpleSum{rowHeadN,2}(rowRevN,1)=0;%-----------匹配长度
            elseif simpleSum{rowHeadN,2}(rowRevN,1)-simpleSum{rowHeadN,1}(rowRevN,2)>spacer || simpleSum{rowHeadN,1}(rowRevN,1)-simpleSum{rowHeadN,2}(rowRevN,2)>spacer
                simpleSum{rowHeadN,2}(rowRevN,1)=0;%-----------间隔
            end
        end
    end
end
%-------------------------------------------------------


rowHS=0;   %------------统计非0行数----------
for rowHeadN=1:RowHead
    [rowRev,~]=size(simpleSum{rowHeadN,2});
    for rowRevN=1:rowRev
        if simpleSum{rowHeadN,2}(rowRevN,1)>0;
            rowHS=rowHS+1;
        end
    end
end
simpleSum3=cell(rowHS,3);
rowH=0;%-------------非0位置个数
for rowHeadN=1:RowHead
    [rowRev,~]=size(simpleSum{rowHeadN,2});
    for rowRevN=1:rowRev
        if simpleSum{rowHeadN,2}(rowRevN,1)>0;
            rowH=rowH+1;
            simpleSum3{rowH,1}=simpleSum{rowHeadN,1}(rowRevN,:);
            simpleSum3{rowH,2}=simpleSum{rowHeadN,2}(rowRevN,:);
            simpleSum3{rowH,3}=simpleSum{rowHeadN,3}(rowRevN);
        end
    end
end
clear simpleSum
