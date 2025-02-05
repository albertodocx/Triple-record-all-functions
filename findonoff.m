
function boutsidx = findonoff(idxvar)

  [~, locs] = findpeaks(double(idxvar));
    
    if idxvar(1) == 1
        locs = [1, locs];
    end
    
    if isempty(locs)
           disp("Warning: No events found")
           boutsidx = [];
    else
           % Get onset and offset of each event saved in locs:

       for kk = 1:length(locs)
           if locs(kk) == locs(end)
               onset(kk) = locs(kk);
               if idxvar(end) == 0
                   offset(kk) = locs(kk) + find(idxvar(locs(kk):end) == 0, 1, 'first') - 1 ;
               else
                   offset(kk) = length(idxvar)-1;
               end
           else
               onset(kk) = locs(kk);
               offset(kk) = locs(kk) + find(idxvar(locs(kk):end) == 0, 1, 'first') - 2 ;
           end
       end
    

           boutsidx = [onset.', offset.'];
    end


end