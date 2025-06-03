%  The code is based on the following paper:
%  Zhou, L., Jin, K. and Tang, Z.
%  Fault tolerant and load balanced deployment of access points in WLANs: 
%  a heuristic algorithm with multi-AP and single-AP search
%  International Journal of Web Information Systems
%  https://doi.org/10.1108/IJWIS-02-2025-0056

%% PSOMSAS for AP deployment problem 
function [gbest_Fitness,gbest,CNVG]=PSOMSAS_4APDP(aps_num,N,T,lb,ub,dim,fobj)
gbest=zeros(aps_num,dim);
gbest_Fitness=inf;

%Initialize the locations: X has size [N, aps_num, dim]
X=initialization(aps_num,N,dim,ub,lb);

CNVG=zeros(1,T);

t=0; % Loop counter
fitness = ones(N,1).*inf;
pbest_fitness = fitness;
pbest = zeros(N, aps_num, dim);
V = zeros(N, aps_num, dim);
Wmax = 0.9; Wmin = 0.4;
c1 = 2; c2 = 2;
Mt = 0.7;

while t<T
    w = Wmax - (Wmax - Wmin)*(t/T);
    for i=1:N
        % Check boundaries
        for ap=1:aps_num
            for d=1:dim
                if X(i,ap,d) > ub(d)
                    X(i,ap,d) = ub(d);
                elseif X(i,ap,d) < lb(d)
                    X(i,ap,d) = lb(d);
                end
            end
        end

        % Flatten AP information to 1D vector for fitness evaluation
        Xi = reshape(X(i,:,:), [1, aps_num*dim]);
        fitness(i) = fobj(Xi);

        % Update gbest
        if fitness(i) < gbest_Fitness
            gbest_Fitness = fitness(i);
            gbest = squeeze(X(i,:,:));
        end
        % Update pbest
        if fitness(i) < pbest_fitness(i)
            pbest_fitness(i) = fitness(i);
            pbest(i,:,:) = X(i,:,:);
        end
    end

    if rand < Mt
        % Multi-AP Search
        for i=1:N
            V(i,:,:) = w * V(i,:,:) ...
                       + c1 * rand() * (pbest(i,:,:) - X(i,:,:)) ...
                       + c2 * rand() * (gbest - X(i,:,:));
            X(i,:,:) = X(i,:,:) + V(i,:,:);
        end
    else
        % Single-AP Search
        for i = 1:N
            j = ceil(rand*N);
            while j==i
                j = ceil(rand*N);
            end
            ap = ceil(rand*aps_num);
            d = ceil(rand*dim);
            theta = (2*rand - 1)*pi;
            pre = X(i,ap,d);
            X(i,ap,d) = X(i,ap,d) + abs(X(j,ap,d)-X(i,ap,d)) * theta * cos(theta);
            if X(i,ap,d) > ub(d) || X(i,ap,d) < lb(d)
                X(i,ap,d) = pre;
            end
        end
    end
    t = t + 1;
    CNVG(t) = gbest_Fitness;
end


function Positions = initialization(aps_num, SearchAgents_no, dim, ub, lb)
Positions = zeros(SearchAgents_no, aps_num, dim);
for ap = 1:aps_num
    for d = 1:dim
        Positions(:,ap,d) = rand(SearchAgents_no,1)*(ub(d)-lb(d)) + lb(d);
    end
end