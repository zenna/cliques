hold off
%plotyy(1:24,tsev/2,1:24,tcommsize)
plot(1:24,tcommsize/960+0.5,'bo')
hold on
plot(1:24,tsev/2,'c-')

P=diag(sum(adj).^-1)*adj;


comm=[20 40 160 480];

Q=P(1:comm(1),1:comm(1));
for i=1:24
[sev(i)] = sev0(Q^i);
end
plot(1:24,sev/2,'r')

Q=P(1:comm(2),1:comm(2));
for i=1:24
[sev(i)] = sev0(Q^i);
end
plot(1:24,sev/2,'g')

Q=P(1:comm(3),1:comm(3));
for i=1:24
[sev(i)] = sev0(Q^i);
end
plot(1:24,sev/2,'m')

Q=P(1:comm(4),1:comm(4));
for i=1:24
[sev(i)] = sev0(Q^i);
end
plot(1:24,sev/2,'k')

%%
legend('Optimal Community','Optimal Severability','20-severability', ...
    '40-severability','160-severability','480-severability')
ylabel('Severability (or rescaled community size)')
xlabel('Markov Time')
title('Severability curves of 20-40-160-480 multilevel random graph, with errors');
