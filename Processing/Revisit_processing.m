%%
clear all; close all; clc;
load("optical_full5.mat")

% window = cell(1,n_debris); % Same as windows but just in case
% difference = cell(1,n_debris);
% idx = cell(1,n_debris);
% Trevisit = 7.2*3600;
% for u = 1:n_debris
%     times = timing{u};
%     idx{u} =1;
%     window{u} = 0;
%     % Define the maximum observation time per debris
%     %max_time = max(windows{u});
%     if isempty(timing{u})
%         continue;
%     end
%     tiempo = 0; 
%     t0=times(1);  
% 
%         for z = 2:length(times)
%             tiempo = abs(times(z-1)-t0);
%              if abs(times(z-1)-times(z))>Trevisit
%                  window{u}=[tiempo,window{u}];
%                  t0=times(z);
%                  difference{u} = [abs(times(z-1)-times(z)),difference{u}];
%                  idx{u} = [z,idx{u}];
%              elseif z==length(times)
%                  window{u}=[tiempo,window{u}];
%                  idx{u} = [z,idx{u}];
%              end
%         end 
% end

% Trevisit = 7.2*3600;
% 
% instant = cell(1,n_debris);
% for d = 1:n_debris
%     tim = timing{d};
%     debris_window = flip(window{d});
%     time_index = flip(idx{d});
% 
%     for t = 1:length(time_index)-1
%         for t2 = t+1:length(time_index)
%             if abs(tim(time_index(t2))-tim(time_index(t)))<=Trevisit
%                     if isempty(instant{d})
%                         instant{d} = [time_index(t),instant{d}];
%                         lastvisit = tim(time_index(t2));
%                     else
%                         if abs(lastvisit-tim(time_index(t)))<=Trevisit
%                             instant{d} = [time_index(t),instant{d}];
%                             lastvisit = tim(time_index(t2));
%                         end
%                     end
%             else
%                 break
%             end
%         end
%     end
% 
% end
% sum(~cellfun('isempty', instant))/n_debris
% %sum(cellfun(@(x) isvector(x) && length(x) > 100, instant))
% %%
% revisit = cell(1,n_debris);
% for deb = 1:n_debris
%     temps = timing{deb};
%     ins = instant{deb};
%     if ~isempty(ins)
%         t1 = temps(end);
%         t2 = temps(ins(end));
%         if abs(t1-t2)<=Trevisit
%             revisit{deb} = [ins,revisit{deb}];
%         end
%     end
% end
% sum(~cellfun('isempty', revisit))/n_debris

Trevisit = 7.2*3600;
revisit = cell(1,n_debris);
lastime = 0;

for deb = 1:n_debris
    temps = flip(timing{deb});
    for t0 = 1:length(temps)
        for t2 = t0:length(temps)
            if abs(temps(t0)-temps(t2))>1 && abs(temps(t0)-temps(t2))<=Trevisit
                if temps(t0)-lastime <=Trevisit
                    revisit{deb} = [1,revisit{deb}];
                    lastime = temps(t2);
                break;
                end
            end
        end
    end
end

check = zeros(1,n_debris);
for index = 1:n_debris
    debris = timing{index};
    rev = revisit{index};
    if length(rev)>= Tsim/Trevisit
        check(1,index) = 1;
    end
end
sum(check>0)

%%
percentage = [72 75 74 79 81];

mean_perc = mean(percentage);
std_perc = std(percentage);

FPC = sqrt((500000-1000)/(500000-1));
standard_perc = FPC*std_perc/sqrt(500000);

