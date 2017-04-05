%processes epr signals into something understandable
function [Gvalues,Avalues,Minutes,TextNames]=EPR

%get the list of txt files
a=dir('*.txt');

%first, extract the time in minutes from the .par file
for i=1:size(a,1)
    %get the .txt filename for i
    Name=a(i).name;
    
    %first,extract just the name, not the extension, then add par to the
    %end
    justname='[a-z_A-Z0-9~]*\.';
    filenamepart=regexp(Name,justname,'match');
    TimeFileName=strcat(filenamepart,'par');
    TextNames{i}=Name;
    
    %error catching to make sure the .par file exists
    b=dir(TimeFileName{1,1});
    if isempty(b)
        continue
    else
        filetext=fileread(TimeFileName{1,1});
    end
    
    %extract the time from the .par file, convert to number then minutes.
    expr='[0-9]*:[0-9]*';
    time=regexp(filetext,expr,'match');
    hours=str2num(time{1,1}(1,1:2));
    minutes=str2num(time{1,1}(1,4:5));
    Minutes(i)=60.*hours+minutes;
end

%next, extract the data.  It's in a different for loop just for easy of
%viewing.
for i=1:size(a,1)
    Name=a(i).name;
    
    fid=fopen(Name);
    frewind(fid);
    data=textscan(fid,'%f%f%s','HeaderLines',6);
    indices=data{1,1};
    values=data{1,2};
    cell=data{1,3};
    Gvalues(:,i)=values;
    for j=1:size(indices,1)
        Avalues(j,i)=str2num(cell{j,1});
    end
    
    fclose(fid);
end
end