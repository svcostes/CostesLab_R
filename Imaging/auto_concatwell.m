function auto_concatwell
cur_dir ='F:\Images';
imdir=dir([cur_dir '\' 'P*']);
imdir(~[imdir.isdir])=[];

A = datetime(2016,08,22,0,0,0,'Format','dd-MMM-yyyy HH:mm:ss');
for K = 1 : length(imdir)
  try   
      thisdir = imdir(K).name;
      average_well= dir(fullfile(cur_dir, thisdir, '*average_well*'));
      B=imdir(K).date;
      
      if isempty(average_well) & B>A 
        
        num= textscan(thisdir(2:end),'%d %s','delimiter','_');
        concatenate_well(double(num{1}));
        fprintf('The final file average_well has been created at F/%s \n',thisdir);
        
      end
      
  end
  
end