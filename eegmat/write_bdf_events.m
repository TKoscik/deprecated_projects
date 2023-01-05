function write_bdf_events(infile, outfile, new_events)

if exist(infile) ~= 2
    error('BDF inpout file not found');
end
if exist(outfile) == 2
    error('BDF output file already exists');
end

hdr = read_bdf_hdr(infile);

fin = fopen(infile, 'r');
fout = fopen(outfile, 'a');

eventts = zeros(hdr.n_records*hdr.n_samples(end), 3);
for i = 1:numel(new_events)
  for j = 1:numel(new_events{i}.sample)
    eventts(new_events{i}.sample(j)) = new_events{i}.marker;
  end
end
% eventts = reshape(eventts', [], 1);

%%
hdr_vals = fread(fin, hdr.hdr_bytes);
fwrite(fout, hdr_vals);

for i = 1:hdr.n_records  
  record = fread(fin, sum(hdr.n_samples)*3);
  record((numel(record)-hdr.n_samples(end)*3+1):numel(record)) = reshape(eventts(((i-1)*hdr.n_samples(end)+1):(i*hdr.n_samples(end)),:)', [],1);  
  fwrite(fout, record);
end
fclose(fin);
fclose(fout);