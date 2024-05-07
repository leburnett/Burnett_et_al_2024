% %% read the file and extract cluster id and depth info
function ids_info = get_all_info(cluster_info_gen_file)
    f=fopen(cluster_info_gen_file);
    %first line is the header
    t=fgetl(f);
    ids_info=[];
    while ~feof(f)
        t=fgetl(f);
        strline=strsplit(t);
            id=str2num(strline{1});
            amp=str2num(strline{2});
%             kslabel=str2num(strline{4});
            chan=str2num(strline{6});
            depth=str2num(strline{7});
            numspikes=str2num(strline{10});
            ids_info=[ids_info; [id,amp, chan, depth, numspikes]];        
    end
    fclose(f);
end