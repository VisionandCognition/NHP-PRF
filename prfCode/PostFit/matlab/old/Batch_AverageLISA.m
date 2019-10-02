T = [0 1 2 4 5 10];

% for ti=1:length(T)
%     if ti==0%1
%         AverageResults_LISA('eddy',T(ti),1)
%     else
%         AverageResults_LISA('eddy',T(ti),0)
%     end
% end
% 
% T = [1 2 4 5 10];
% for ti=1:length(T)
%     if ti==0%1
%         AverageResults_LISA('danny',T(ti),1)
%     else
%         AverageResults_LISA('danny',T(ti),0)
%     end
% end

% all at once (saves time because less loading)
AverageResults_LISA('eddy',T,0)
AverageResults_LISA('danny',T,0)
