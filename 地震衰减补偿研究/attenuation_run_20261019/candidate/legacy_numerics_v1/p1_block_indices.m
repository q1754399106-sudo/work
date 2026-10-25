function ids=p1_block_indices(available,count,block,stream)
% Circular moving blocks; all traces share row indices.
starts=randi(stream,available,ceil(count/block),1);
ids=reshape(mod(starts-1+(0:block-1),available)'+1,[],1);
ids=ids(1:count);
end
