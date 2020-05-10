function [MatcInDyPart,Aend,Bend]=IRCodlocDyMatc(A,B)           %----------------返回各个位置----------
misInter=-4;
La=length(A);%-----------用于计算总长度--
Lb=length(B);%-----------用于计算总长度--
col=La+1;
row=Lb+1;
%-----------------作出得分矩阵---------
F=zeros(row,col);Fmax=0;
iMax=row;jMax=col;
for j=2:col
    for i=2:row
        Ch1=F(i-1,j-1)+lettercmp(B(i-1),A(j-1)) ;
        Ch2=F(i-1,j)+misInter;
        Ch3=F(i,j-1)+misInter;
        F(i,j)=max([Ch1,Ch2,Ch3,0]);
%------------------记录最大得分处
        if  Fmax<F(i,j)
            Fmax=F(i,j);
            iMax=i;
            jMax=j;               
        end
%----------------------------------------        
    end
end
iMaxBend=iMax-1;jMaxAend=jMax-1;   %-------------------记录终点位置，最大值前一个为匹配的-----
%--------------反向读取，构建最佳匹配----
matcN=0; misInterN=0;%---------------------计数匹配/错配、插入的个数----错配和插入的位置---% misInterN=0;       %---------------------错配/插入的个数
while iMax>1 && jMax>1
      Sco=F(iMax,jMax);
      ScoDia=F(iMax-1,jMax-1);
      ScoLef=F(iMax,jMax-1);
      ScoUp=F(iMax-1,jMax);
      cmp=lettercmp(B(iMax-1),A(jMax-1)) ;
      iMax2=iMax;jMax2=jMax; %-----------------当两个标记都不改变时，则停止循环---用于极不相似的序列---------------
%----------------------回溯过程，至0值处终止
      if Sco==ScoUp+misInter  
          iMax=iMax-1;
          misInterN=misInterN+1;
      elseif Sco==ScoLef+misInter  
          jMax=jMax-1;  
          misInterN=misInterN+1;
      elseif Sco==ScoDia+cmp
         if cmp>0
             matcN=matcN+1;
         else
             misInterN=misInterN+1;
         end  
         iMax=iMax-1;
         jMax=jMax-1; 
      end 
%-----------------------------------------------------------         
%------------当不能到达0处时，跳出循环---------------
      if  iMax2==iMax && jMax2==jMax   %-------------已经不能回溯了---------
          break;
      end
end
%-----------------------------------------iMax,jMax最终值为起点位置----------------------
 %------------以碱基个数总数计算不同得分:各得分用小数-----
 MatcInDyPart=matcN/(matcN+misInterN); %----------------正确匹配的百分比----------------------
 Aend=[jMax ; jMaxAend];
 Bend=[iMax ; iMaxBend];
