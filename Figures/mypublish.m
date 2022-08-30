function mypublish(filename)
% A Custom publish function that includes the line
% number for each line in the function
% Usage: mypublish('myfun.m')

function_options.format='html'
function_options.evalCode=false;
function_options.showCode=true;

copyfile(filename, 'mytemp.m')
fid1 = fopen(filename);
fid2 = fopen('mytemp.m', 'w');
i=1;
while 1

    tline = fgetl(fid1);
    if ~ischar(tline)
        break;
    end
    if ~isempty(tline)
        tline = [num2str(i),' ', tline];
        fprintf(fid2,'%s',tline);
        fprintf(fid2,'%s\n','');
        i = i + 1;
    end
end
fclose(fid1);
fclose(fid2);
publish(filename, function_options); % publish original function without line numbers
publish('mytemp',function_options); % publish temp function with line numbers
delete('mytemp.m');  % delete temp M -file