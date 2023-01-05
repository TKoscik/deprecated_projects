function [events, event_ts] = read_bdf_events(bdffile)

if exist(bdffile) ~= 2
    error('BDF file not found');
end

hdr = read_bdf_hdr(bdffile);

fid = fopen(bdffile, 'r');

event_ts = zeros(hdr.n_records * hdr.n_samples(end), 1);
for i = 1:hdr.n_records
    offset = hdr.hdr_bytes + sum(hdr.n_samples(1:end-1))*3 + sum(hdr.n_samples)*(i-1)*3;
    fseek(fid, offset, 'bof');
    temp = fread(fid, 3*hdr.n_samples(end));
    event_ts(((i-1)*hdr.n_samples(end)+1):(i*hdr.n_samples(end))) = temp(1:3:end) + temp(2:3:end)*256;    
end

event_ls = unique(event_ts);
event_ls = sort(event_ls(event_ls~=0));
for i = 1:numel(event_ls)
    events{i}.marker = event_ls(i);
    
    evtemp = find(event_ts == event_ls(i));
    evdiff = [evtemp(1); diff(evtemp)];
    evlocs = find(evdiff > 1);
    
    events{i}.sample = evtemp(evlocs);
    events{i}.time = (1/hdr.n_samples(end)) * (events{i}.sample-1);
    events{i}.duration = [diff(evlocs); numel(evtemp) - evlocs(numel(evlocs))];
end

fclose(fid);