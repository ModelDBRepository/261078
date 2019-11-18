% file_reader.m
% makes arrays available in matlab
% see name_array to see order that traces are stored in their respective arrays.
files=dir('*.dat');
data_array={};
data_index=1;
name_array={};
zero_valued={};
zero_index=1;
nonzero_indicies = {};
nonzero_index = 1;

ulb_array={};
upy_array={};
ulb_index=1;
upy_index=1;

for i=1:length(files)
    filename=files(i).name;
    if length(filename)>4
        if strcmp(filename(end-3:end),'.dat')
            cmd =['load ' filename ';'];
            eval(cmd)
            variable_name=filename(1:end-4);
            cmd=['data_array{' num2str(data_index) '}=' variable_name ';'];
            eval(cmd)
            name_array{data_index} = variable_name;
            cmd=['current_trace = ' variable_name ';'];
            eval(cmd)
            if sum(abs(current_trace))==0
                zero_valued{zero_index}=variable_name;
                zero_index = zero_index + 1;
            else
                nonzero_indicies{nonzero_index}=data_index;
                nonzero_index = nonzero_index + 1;
            end
            if strcmp(filename(1:3),'uLB')
                ulb_array{ulb_index}=current_trace;
                ulb_index=ulb_index + 1;
            end
            if strcmp(filename(1:3),'uPY')
                upy_array{upy_index}=current_trace;
                upy_index=upy_index + 1;
            end
            data_index = data_index+1;
        end
    end
end


figure
hold on
plot(ulb_array{7})
plot(ulb_array{9})
title('uLB3_01.dat, uLB4_01.dat','Interpreter','None')
