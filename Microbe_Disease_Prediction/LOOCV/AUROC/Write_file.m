function Write_file( F )
       [rows, cols] = size (F);    
       fid=fopen (  'PredHMDA.xls', 'w');
       for i = 1:rows
           for j = 1:cols
               fprintf(fid, '%s\t', F{i,j});
           end
           fprintf(fid, '\n');
       end
       fclose(fid);      


end

