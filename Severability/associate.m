%{
keyfile=fopen('WordAssoc/key2.txt');
cues=textscan(keyfile, '%d %s');
fclose(keyfile);
A=adjacency_matrix('WordAssoc/words.adj');
P=diag(sum(A,2).^-1)*A;
%}
%%{
time=2;
start=4774;
neighbors=find(A(start,:));
for i=neighbors
    [commrun commsize sev]=sev_node([start i ], A, time, -1, 75);
    [sev retention mixing ] = sev0(P(commrun(1:commsize),commrun(1:commsize))^time);
    disp([sev retention mixing]);
    
    for i=[commrun(1:commsize)]
        disp(cues{2}{i});
    end
end

%}
