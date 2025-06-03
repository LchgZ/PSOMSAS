%  The code is based on the following paper:
%  Zhou, L., Jin, K. and Tang, Z.
%  Fault tolerant and load balanced deployment of access points in WLANs: 
%  a heuristic algorithm with multi-AP and single-AP search
%  International Journal of Web Information Systems
%  https://doi.org/10.1108/IJWIS-02-2025-0056

%% PSOMSAS for general problem
function [gbest_Fitness,gbest,CNVG]=PSOMSAS_4GP(N,T,lb,ub,dim,fobj)
gbest=zeros(1,dim);
gbest_Fitness=inf;

%Initialize the locations
X=initialization(N,dim,ub,lb);

CNVG=zeros(1,T);

t=0; % Loop counter
fitness = ones(N,1).*inf;
pbest_fitness = fitness;
pbest = zeros(N,dim);
V = zeros(N,dim);
Wmax = 0.9; Wmin = 0.4;
c1 = 2; c2 = 2;
Mt = 0.7;

while t<T
    w = Wmax - (Wmax - Wmin)*(t/T);
    for i=1:N
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness(i)=fobj(X(i,:));
        % Update the location
        if fitness(i)<gbest_Fitness
            gbest_Fitness=fitness(i);
            gbest=X(i,:);
        end
        if fitness(i)<pbest_fitness(i)
            pbest_fitness(i) = fitness(i);
            pbest(i,:) = X(i,:);
        end
    end
    
    % Update the location
    if rand < Mt
        % Multi-Dimension Search
        for i=1:size(X,1)
            V(i,:) = w*V(i,:) + c1*rand(1,dim).*(pbest(i,:)-X(i,:)) + c2*rand(1,dim).*(gbest-X(i,:));
            X(i,:) = X(i,:) + V(i,:);
        end
    else
        % Single-Dimension Search
        for i = 1:N
            j = ceil(rand*N);
            while j==i
                j= ceil(rand*N);
            end
            k = ceil(rand*dim);
            theta = (2*rand-1)*pi;
            pre = X(i,k);
            X(i,k) = X(i,k)+abs(X(j,k)-X(i,k))*theta*cos(theta);
            if X(i,k)>ub(k) || X(i,k)<lb(k)
                X(i,k) = pre;
            end
        end
        t=t+1;
        CNVG(t)=gbest_Fitness;
    end   
end

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % number of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end