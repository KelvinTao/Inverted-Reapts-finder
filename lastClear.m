function outS=lastClear(s)
     [row,col]=size(s);row2=row;
     for rowN=1:row
        for rowNext=rowN+1:row
            if s(rowN,1)==s(rowNext,3)&&s(rowN,2)==s(rowNext,4)&&s(rowN,3)==s(rowNext,1)&&s(rowN,4)==s(rowNext,2)
                s(rowNext,1)=0;
                row2=row2-1;
                break;
            end
        end
     end
     outS=zeros(row2,col);row2N=0;
     for rowN=1:row
          if s(rowN,1)>0
                row2N=row2N+1;
                outS(row2N,:)=s(rowN,:);
          end
     end
end
