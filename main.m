function main()
clear; clc; close all; rng(42);
total_channels=160; reuse_dist=4; grid_size=27;
fixed_caps=[0,5,6,7,8,9,10]; load_points=4:2:18;
call_duration_avg=98; warmup_time=120; sim_time=420; total_time=warmup_time+sim_time;

R=zeros(5,length(fixed_caps),length(load_points));
fprintf('Building interference matrix for %dx%d grid...\n',grid_size,grid_size);
[interf_cells, ring_channels] = build_interference_list(grid_size, reuse_dist);

for f_idx=1:length(fixed_caps)
    Nf=fixed_caps(f_idx); Nd=total_channels-Nf*(reuse_dist^2);
    fprintf('\n=== Nf=%d (Nd=%d) ===\n',Nf,Nd);
    for L_idx=1:length(load_points)
        traffic_load=load_points(L_idx);
        lambda_global=(traffic_load/call_duration_avg)*(grid_size^2);
        [R(1,f_idx,L_idx),R(2,f_idx,L_idx),R(3,f_idx,L_idx),R(4,f_idx,L_idx),R(5,f_idx,L_idx)]=...
            run_sim(lambda_global,grid_size,interf_cells,ring_channels,Nf,Nd,call_duration_avg,total_time,warmup_time);
        fprintf('  Load %.0f -> Block:%.1f%%, Carried:%.2f, DynAcc:%.1f%%, Switch:%.1f%%, Occ:%.1f%%\n',...
            traffic_load,R(1,f_idx,L_idx),R(2,f_idx,L_idx),R(3,f_idx,L_idx),R(4,f_idx,L_idx),R(5,f_idx,L_idx));
    end
end

% Compute Table I data for plotting
tableI_data = [];
idx_fixed=find(fixed_caps==10,1); idx_hybrid=find(fixed_caps==8,1);
if ~isempty(idx_fixed) && ~isempty(idx_hybrid)
    target_blocking=[1,2,5];
    fprintf('\n===== TABLE I =====\nBlock%% | Channels Required (FIXSYS) | (Hybrid: F=8+2 dynamic)\n');
    for tb=1:3
        blk_hybrid = squeeze(R(1,idx_hybrid,:));
        carried_hybrid = squeeze(R(2,idx_hybrid,:));
        [blk_unique,ia_unique]=unique(blk_hybrid,'stable');
        carried_unique=carried_hybrid(ia_unique);
        carried_at_target=interp1(blk_unique,carried_unique,target_blocking(tb),'linear','extrap');
        A=carried_at_target;
        Bt=target_blocking(tb)/100;
        Cr=1;
        while true
            num=A^Cr/factorial(Cr); den=0;
            for i=0:Cr, den=den+A^i/factorial(i); end
            if num/den<=Bt || Cr>30, break; end
            Cr=Cr+1;
        end
        tableI_data = [tableI_data; target_blocking(tb), Cr];
        fprintf('  %2d%% -> %2d\n',target_blocking(tb),Cr);
    end
    fprintf('注: 混合系统8固定+2动态, 纯固定系统需更多信道达相同承载\n');
end

plot_all(load_points,fixed_caps,R,tableI_data);

fprintf('\n\n===== TABLE =====\nNf\tLoad\tBlock%%\tCarried\tDynAcc%%\tSwitch%%\tOcc%%\n');
for f_idx=1:length(fixed_caps)
    for L_idx=1:length(load_points)
        fprintf('%d\t%.0f\t%.1f\t%.2f\t%.1f\t\t%.1f\t\t%.1f\n',...
            fixed_caps(f_idx),load_points(L_idx),R(1,f_idx,L_idx),R(2,f_idx,L_idx),R(3,f_idx,L_idx),R(4,f_idx,L_idx),R(5,f_idx,L_idx));
    end
end
end

function [blocking_pct,carried_per_chan,dyn_attempts_pct,switched_pct,occupancy_pct]=...
    run_sim(lambda_global,grid_size,interf_cells,ring_channels,Nf,Nd,call_duration_avg,total_time,warmup_time)
num_cells=grid_size^2;
central_start = floor((grid_size-15)/2) + 1;
central_end = central_start + 14;
central_cells = false(num_cells, 1);
for i = 1:grid_size
    for j = 1:grid_size
        if i >= central_start && i <= central_end && j >= central_start && j <= central_end
            idx = (i-1)*grid_size + j;
            central_cells(idx) = true;
        end
    end
end
num_central = 225;

fixed_in_use=zeros(num_cells,1);
cell_dyn=cell(num_cells,1); for i=1:num_cells, cell_dyn{i}=false(Nd,1); end

stat_arr=0; stat_blk=0; stat_dyn=0; stat_sw=0; stat_car=0; ga_central=0;

max_ev=1000000; eq=zeros(max_ev,5); eqc=0;
n_arr=min(ceil(total_time*lambda_global*1.5),500000);
arr_t=cumsum(exprnd(1/lambda_global,n_arr,1));
arr_c=randi(num_cells,n_arr,1); ai=1;

eqc=1; eq(1,:)=[arr_t(1),1,arr_c(1),0,0]; ai=2;
ct=0; ga=0;

while eqc>0 && ct<total_time
    evt=eq(1,:); eq(1,:)=eq(eqc,:); eqc=eqc-1;
    idx=1;
    while true
        l=2*idx; r=2*idx+1; s=idx;
        if l<=eqc && eq(l,1)<eq(s,1), s=l; end
        if r<=eqc && eq(r,1)<eq(s,1), s=r; end
        if s==idx, break; end
        tmp=eq(idx,:); eq(idx,:)=eq(s,:); eq(s,:)=tmp; idx=s;
    end
    t=evt(1); if t>total_time, break; end
    dt=t-ct; ct=t; is_stat=(t>=warmup_time);
    if is_stat, stat_car=stat_car+ga_central*dt; end
    
    if evt(2)==1
        cid=evt(3);
        if central_cells(cid)
            stat_arr=stat_arr+1;
        end
        served=false;
        if fixed_in_use(cid)<Nf
            fixed_in_use(cid)=fixed_in_use(cid)+1; served=true; ch_t=1; ch_i=0;
        elseif Nd>0
            ch_assigned = -1;
            for ch = 1:Nd
                available = true;
                if cell_dyn{cid}(ch)
                    available = false;
                end
                if available
                    for k=1:length(interf_cells{cid})
                        if cell_dyn{interf_cells{cid}(k)}(ch)
                            available = false;
                            break;
                        end
                    end
                end
                if available
                    ch_assigned = ch;
                    break;
                end
            end
            if ch_assigned >= 0
                cell_dyn{cid}(ch_assigned)=true; served=true; ch_t=2; ch_i=ch_assigned;
                if is_stat && central_cells(cid), stat_dyn=stat_dyn+1; end
            end
        end
        if served
            ga=ga+1;
            if central_cells(cid)
                ga_central=ga_central+1;
            end
            dur=max(1,round(exprnd(call_duration_avg)));
            eqc=eqc+1; eq(eqc,:)=[t+dur,2,cid,ch_t,ch_i];
            idx=eqc;
            while idx>1
                p=floor(idx/2);
                if eq(p,1)>eq(idx,1), tmp=eq(p,:); eq(p,:)=eq(idx,:); eq(idx,:)=tmp; idx=p;
                else break; end
            end
        elseif is_stat && central_cells(cid)
            stat_blk=stat_blk+1;
        end
        if ai<=n_arr && arr_t(ai)<=total_time
            eqc=eqc+1; eq(eqc,:)=[arr_t(ai),1,arr_c(ai),0,0];
            idx=eqc;
            while idx>1
                p=floor(idx/2);
                if eq(p,1)>eq(idx,1), tmp=eq(p,:); eq(p,:)=eq(idx,:); eq(idx,:)=tmp; idx=p;
                else break; end
            end
            ai=ai+1;
        end
    else
        cid=evt(3); ch_t=evt(4); ch_i=evt(5);
        ga=ga-1;
        if central_cells(cid)
            ga_central=ga_central-1;
        end
        if ch_t==1
            fixed_in_use(cid)=fixed_in_use(cid)-1;
            ad=find(cell_dyn{cid},1);
            if ~isempty(ad)
                cell_dyn{cid}(ad)=false; fixed_in_use(cid)=fixed_in_use(cid)+1;
                for q=1:eqc
                    if eq(q,2)==2 && eq(q,3)==cid && eq(q,4)==2 && eq(q,5)==ad
                        eq(q,4)=1; eq(q,5)=0; break;
                    end
                end
                if is_stat && central_cells(cid), stat_sw=stat_sw+1; end
            end
        elseif ch_i>0
            cell_dyn{cid}(ch_i)=false;
        end
    end
end
sd=max(ct-warmup_time,1);
blocking_pct=(stat_blk/max(stat_arr,1))*100;
carried_per_chan=(stat_car/sd)/num_central;
dyn_attempts_pct=(stat_dyn/max(stat_arr,1))*100;
switched_pct=(stat_sw/max(stat_arr,1))*100;
occupancy_pct=(carried_per_chan/10)*100;
end

function [interf_cells, ring_channels] = build_interference_list(grid_size, reuse_dist)
num_cells=grid_size^2; interf_cells=cell(num_cells,1);
[X,Y]=meshgrid(1:grid_size,1:grid_size); coords=[X(:),Y(:)]; min_d2=reuse_dist^2;
for i=1:num_cells
    nb=[];
    for j=1:num_cells
        if i==j, continue; end
        dx=coords(i,1)-coords(j,1); dy=coords(i,2)-coords(j,2);
        if dx^2+dy^2<min_d2, nb(end+1)=j; end
    end
    interf_cells{i}=nb;
end
Nd_max = 160;
ring_channels = cell(1, 1);
ring_channels{1} = 1:Nd_max;
end

function plot_all(load_points,fixed_caps,R,tableI_data)
set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10,'DefaultLineLineWidth',1.5);
n_plots=length(fixed_caps);
co=lines(n_plots); mk={'o','s','d','^','v','p','h'};
pl_str=arrayfun(@(x) num2str(x), fixed_caps, 'UniformOutput', false);

figure('Name','Paper Results','Position',[50,50,1400,900]);

subplot(3,3,1); hold on; grid on;
for k=1:n_plots
    plot(load_points,squeeze(R(2,k,:)),['-' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'MarkerFaceColor',co(k,:),'DisplayName',sprintf('Nf=%s',pl_str{k}));
end
xlabel('Offered (Erl/cell)'); ylabel('Carried (Erlangs/cell)'); title('Fig3a: Carried vs Offered'); legend('Location','northwest'); box on;

subplot(3,3,2); hold on; grid on;
for k=1:n_plots
    semilogy(load_points,squeeze(R(1,k,:)),['--' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'DisplayName',sprintf('Nf=%s',pl_str{k}));
end
xlabel('Offered (Erl/cell)'); ylabel('Blocking (%)'); title('Fig3b: Blocking vs Offered');
set(gca,'YScale','log'); ylim([0.1,100]); legend('Location','northwest'); box on;

subplot(3,3,3); hold on; grid on;
for k=1:n_plots
    plot(squeeze(R(1,k,:)),squeeze(R(2,k,:)),['-' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'MarkerFaceColor',co(k,:),'DisplayName',sprintf('Nf=%s',pl_str{k}));
end
xlabel('Blocking (%)'); ylabel('Carried (Erlangs/cell)'); title('Fig4: Carried vs Blocking'); legend('Location','southeast'); box on;

subplot(3,3,4); hold on; grid on;
for k=1:n_plots
    if fixed_caps(k)>0
        fc=squeeze(R(2,k,:)).*(1-squeeze(R(3,k,:))/100);
        plot(squeeze(R(1,k,:)),fc,['-' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'MarkerFaceColor',co(k,:),'DisplayName',sprintf('Nf=%s',pl_str{k}));
    end
end
xlabel('Blocking (%)'); ylabel('Fixed Carried (Erlangs/cell)'); title('Fig5: Fixed Ch Carried'); legend('Location','best'); box on;

subplot(3,3,5); hold on; grid on;
for k=1:n_plots
    if fixed_caps(k)>0 && fixed_caps(k)<10
        plot(squeeze(R(1,k,:)),squeeze(R(3,k,:)),['-' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'MarkerFaceColor',co(k,:),'DisplayName',sprintf('Nf=%s',pl_str{k}));
    end
end
xlabel('Blocking (%)'); ylabel('Dyn Access (%)'); title('Fig6: Dyn Ch Access'); legend('Location','best'); box on;

subplot(3,3,6); hold on; grid on;
for k=1:n_plots
    if fixed_caps(k)>0 && fixed_caps(k)<10
        plot(squeeze(R(1,k,:)),squeeze(R(4,k,:)),['-' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'MarkerFaceColor',co(k,:),'DisplayName',sprintf('Nf=%s',pl_str{k}));
    end
end
xlabel('Blocking (%)'); ylabel('Switched (%)'); title('Fig7: Switching Activity'); legend('Location','best'); box on;

subplot(3,3,7); hold on; grid on;
for k=1:n_plots
    plot(squeeze(R(1,k,:)),squeeze(R(5,k,:)),['-' mk{mod(k-1,length(mk))+1}],'Color',co(k,:),'MarkerSize',6,'MarkerFaceColor',co(k,:),'DisplayName',sprintf('Nf=%s',pl_str{k}));
end
xlabel('Blocking (%)'); ylabel('Occupancy (%)'); title('Fig8: Ch Occupancy'); legend('Location','best'); box on;

% Subplot 8: Table I
subplot(3,3,8); axis off;
if ~isempty(tableI_data)
    t_str = {'TABLE I: Channels Required'; ''; 'Block%   FIXSYS   Hybrid(8+2)'};
    for r = 1:size(tableI_data,1)
        t_str{end+1} = sprintf('  %2d%%      %2d', tableI_data(r,1), tableI_data(r,2));
    end
    t_str{end+1} = '';
    t_str{end+1} = 'Hybrid: 8 fixed + 2 dynamic';
    t_str{end+1} = 'FIXSYS: 10 fixed + 0 dynamic';
    text(0.1, 0.9, t_str, 'FontSize', 10, 'VerticalAlignment', 'top', 'FontName', 'Courier');
end
end
