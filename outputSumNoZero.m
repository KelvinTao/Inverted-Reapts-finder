function outSum=outputSumNoZero(sumI)
   row=0;
   [head,~]=size(sumI);
   for headN=1:head
       [col,~]=size(sumI{headN,1});
       for colN=1:col
           if sumI{headN,2}(colN,1)>0
              row=row+1;%----------找总行数
           end
       end
   end
   outSum=zeros(row,5);rowN=0;
   for headN=1:head
       [col,~]=size(sumI{headN,1});
       for colN=1:col
           if sumI{headN,2}(colN,1)>0
               rowN=rowN+1;
               outSum(rowN,1)=sumI{headN,1}(colN,1);
               outSum(rowN,2)=sumI{headN,1}(colN,2);
               outSum(rowN,3)=sumI{headN,2}(colN,1);
               outSum(rowN,4)=sumI{headN,2}(colN,2);
               outSum(rowN,5)=sumI{headN,3}(colN,1);
           end
       end
   end
end
