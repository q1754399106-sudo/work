function dims=ml2_dims(args)
if isempty(args),dims=[1 1];elseif numel(args)==1,dims=double(args{1});if isscalar(dims),dims=[dims dims];end
else,dims=cellfun(@double,args);end
end
