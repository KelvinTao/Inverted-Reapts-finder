function num=convert(alpha)
         switch alpha
             case 'C'
                 num=0;
             case 'T'
                 num=1;
             case 'A'
                 num=2;
             case 'G'
                 num=3;
             case 'N'
                 num=round(3*rand);%-------------对于不确定的碱基,随机取值
         end
end
