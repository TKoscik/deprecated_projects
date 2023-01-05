function hdr = read_bdf_hdr(bdffile)

if exist(bdffile) ~= 2
    error('BDF file not found');
end

fid = fopen(bdffile, 'r');

hdr.id_ascii       = fread(fid, 1);
hdr.id_code        = fread(fid, 7, '*char')';
hdr.subject_id     = strtrim(fread(fid, 80, '*char')');
hdr.record_id      = strtrim(fread(fid, 80, '*char')');
hdr.start_date     = fread(fid, 8, '*char')';
hdr.start_time     = fread(fid, 8, '*char')';
hdr.hdr_bytes      = str2num(fread(fid, 8, '*char')');
hdr.data_format    = strtrim(fread(fid, 44, '*char')');
hdr.n_records      = str2num(fread(fid, 8, '*char')');
hdr.rec_dur_s      = str2num(fread(fid, 8, '*char')');
hdr.n_channels     = str2num(fread(fid, 4, '*char')');

for i = 1:hdr.n_channels
    hdr.channel_labels{i} = strtrim(fread(fid, 16, '*char')');
end
for i = 1:hdr.n_channels
    hdr.transducer{i} = strtrim(fread(fid, 80, '*char')');
end
for i = 1:hdr.n_channels
    hdr.channel_dim{i} = strtrim(fread(fid, 8, '*char')');
end
for i = 1:hdr.n_channels
    hdr.phys_min(i) = str2num(fread(fid, 8, '*char')');
end
for i = 1:hdr.n_channels
    hdr.phys_max(i) = str2num(fread(fid, 8, '*char')');
end
for i = 1:hdr.n_channels
    hdr.digi_min(i) = str2num(fread(fid, 8, '*char')');
end
for i = 1:hdr.n_channels
    hdr.digi_max(i) = str2num(fread(fid, 8, '*char')');
end
for i = 1:hdr.n_channels
    hdr.prefilter{i} = strtrim(fread(fid, 80, '*char')');
end
for i = 1:hdr.n_channels
    hdr.n_samples(i) = str2num(fread(fid, 8, '*char')');
end
for i = 1:hdr.n_channels
    hdr.reserved{i} = strtrim(fread(fid, 32, '*char')');
end

fclose(fid);