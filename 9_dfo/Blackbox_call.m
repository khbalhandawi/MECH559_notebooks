function [f,g,df,dg] = Blackbox_call(x,P)
    global index
    save_data = P{1};
    sur = P{2};
    use_gradients = P{3};
    p = P{4};
    sur_model = P{5};
    %% Run Blackbox
    index = index + 1;

    % Delete output files (from previous runs)
    try_remove('output.txt');

    command_str = ['-x ', vector2str(x), ' -p ', vector2str(p)];
    if use_gradients
        command_str = [command_str, ' -G'];
    end
    command = join(['bb ',command_str]);

    if ~(sur)
        %%%%%%%%%%%%%%%%%%%%%
        % Real model
        %%%%%%%%%%%%%%%%%%%%%
        status = system(command); % this runs the blackbox
        outs = read_output('output.txt');
        f = outs{1};
        g = outs{2};
        if use_gradients
            df = outs{3};
            dg = [outs{4}'; outs{5}'; outs{6}'];
        else
            df = [];
            dg = [];
        end

        %%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%
        % Surrogate model
        %%%%%%%%%%%%%%%%%%%%%
        [YX, ~] = predictor(x, sur_model);
        if ~isrow(YX) % make sure y is a row vector
            YX = YX';
        end

        f = YX(:,1);
        g = YX(:,2:end);
        df = [];
        dg = [];
        %%%%%%%%%%%%%%%%%%%%%
    end

    if save_data
        log_output('log.txt',x,g,f);
    end
end

function log_output(file, x, g, f)
    fid = fopen(file,'at');
    x_str = vector2str(x);
    g_str = vector2str(g);
    f_str = num2str(f);
    fprintf(fid,'%s %s %s\n',x_str,g_str,f_str);
    fclose(fid);
end

function a = read_output(file)
    fid = fopen(file,'r'); % read the output file
    
    a = {};
    n = 1;
    tline = fgetl(fid);
    while ischar(tline)
        a{n} = cell2mat(textscan(tline,'%f'));
        tline = fgetl(fid);
        n = n + 1;
    end

    fclose(fid);
end

function str = vector2str(array)
    fmt=[repmat('%f ',1,numel(array))];
    str = sprintf(fmt,array);
    str = str(1:end-1);
end

function try_remove(filename)
    % Delete files
    if exist(filename, "file") == 2
      delete(filename)
    end
end