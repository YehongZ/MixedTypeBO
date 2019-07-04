% Plot the averaged IR for synthetic functions
clear all, close all

cost = 100;
N = 2500;
T = 10;

load syn_true_max

mixed50 = [];

for seq = 1:10
    
    fname = ['syn', num2str(seq)];
    groundtruth = syn_max(seq);
    
    load([fname, '_cr50_result'])
    tmp = groundtruth - [zeros(T, 1), result.umax_f];
    mixed50 = [mixed50; tmp(:, (1:cost:N+1))];

end

xl='Cost'; % 3. xlabel
yu = 'log10 Average IR'; % 4. ylabel

mark={'-gv'}; % 6. mark
leg={'MT-PES [50,1]'}; % 5. leg

af = [0:cost:N];
bu = log10([mean(mixed50, 1)]);
err = [std(mixed50)]./sqrt(100);
    
c = mFig(af,bu,xl,yu,mark,leg);
c.lpos = 'northeast'; %position of the legent
c.title = '(a) Synthetic - [Target, aux1]';
c.range = [0, N, -2.2, 0.5];
c.err = err;
% mPlot('plot',c,'');
mPlot('errorbar',c,'');




