%--------------------------------
% Metody Numeryczne - laboratorium 3. 
% Generator sieci polaczen miedzy stronami. 
% Grzegorz Fotyga, ETI, PG
%--------------------------------

function [Edges] = generate_network(N, density)

% randi(max_value, rows, columns)
Edges = randi(N, 2, N * density) ;

Edges(:, find((Edges(1,:) - Edges(2,:))==0))=[];     % usuniecie linkow - gdzie strona wskazuje na sama siebie.

Edges_temp = sortrows(Edges'); % sortowanie, zeby ulatwic usuniecie zdublowanych polaczen.
% sum(A,dim) returns the sum along dimension dim. For example, if A is a matrix, then sum(A,2) returns a column vector containing the sum of each row.
% abs(x) = |x|
d_connections = find(sum(abs(Edges_temp(2:end,:) - Edges_temp(1:end-1,:)),2) == 0);   % wyszukanie zdublowanych polaczen.
Edges_temp(d_connections,:) = [];    % usuniecie zdublowanych polaczen.
Edges = Edges_temp';

Sites_uniq = unique(Edges(1,:));
nr_of_end_sites = 0;

for k = 1:N
    if ~ismember(k, Sites_uniq)
        nr_of_end_sites = nr_of_end_sites + 1;
        Edges = [Edges [[ones(1,N-1) * k]; [1:k-1 k+1:N]]];
    end
end

% disp(['* Liczba stron koncowych: ', num2str( nr_of_end_sites)])

end
