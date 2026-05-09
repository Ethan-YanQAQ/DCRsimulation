% main.m
% 复现论文 "Increasing channel occupancy in large-scale mobile radio systems Dynamic channel REassignment" 
clear; clc; close all;

%% ==================== 参数初始化 ====================
total_channels = 160;       
reuse_dist = 4;             
grid_size = 28;             
% 观测中间 16x16 的区域 (从 7 到 22)
obs_range = 7:22;           
fixed_caps = [0, 1, 2, 3, 4, 5, 8, 10]; 
load_points = 4:2:12;       % 话务量测试点 (Erlangs/cell)
call_duration = 98;         
move_time_threshold = 60;   
sim_time = 600;            

%% ==================== 结果记录容器 ====================
res_blocking_prob = zeros(length(fixed_caps), length(load_points));
res_traffic_carried = zeros(length(fixed_caps), length(load_points));
res_dyn_attempts_pct = zeros(length(fixed_caps), length(load_points));
res_switched_pct = zeros(length(fixed_caps), length(load_points));
res_fixed_carried = zeros(length(fixed_caps), length(load_points));

%% ==================== 主仿真循环 ====================
for f_idx = 1:length(fixed_caps)
    Nf = fixed_caps(f_idx);
    Nd = total_channels - Nf * (reuse_dist^2); 
    
    fprintf('test.m开始仿真 Nf = %d... \n', Nf);
    
    for L_idx = 1:length(load_points)
        traffic_load = load_points(L_idx);
        lambda_cell = traffic_load / call_duration; 
        lambda_global = lambda_cell * (grid_size^2);
        
        fixed_pool = zeros(7, 7, 4, 4, Nf);
        dynamic_pool = zeros(7, 7, Nd);
        active_calls = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
        global_call_id = uint32(0);
        
        stat_total_arrivals = 0;
        stat_blocked = 0;
        stat_dyn_attempts = 0;
        stat_switched = 0;
        stat_carried_samples = 0;
        stat_carried_fixed_samples = 0;
        sample_counts = 0;
        
        for t = 1:sim_time
            % 1. 状态更新
            call_ids = keys(active_calls);
            for i = 1:length(call_ids)
                cid = call_ids{i};
                c = active_calls(cid);
                c.dur = c.dur - 1;
                c.move_timer = c.move_timer - 1;
                
                if c.dur <= 0 || c.move_timer <= 0
                    [fixed_pool, dynamic_pool, active_calls, switched_flag] = ...
                        release_and_reassign(c, fixed_pool, dynamic_pool, active_calls, Nf);
                    
                    if switched_flag && is_in_obs(c.x, c.y, obs_range)
                        stat_switched = stat_switched + 1;
                    end
                    
                    if c.move_timer <= 0 && c.dur > 0
                        c.x = max(1, min(grid_size, c.x + randi([-1, 1])));
                        c.y = max(1, min(grid_size, c.y + randi([-1, 1])));
                        c.move_timer = move_time_threshold;
                        [fixed_pool, dynamic_pool, c, assigned, ~] = ...
                            assign_channel(c, fixed_pool, dynamic_pool, Nf, Nd);
                        if assigned, active_calls(c.id) = c; end
                    end
                else
                    active_calls(cid) = c;
                end
            end
            
            % 2. 呼叫生成 (全网均匀分布)
            num_new_calls = poissrnd(lambda_global);
            for k = 1:num_new_calls
                x = randi([1, grid_size]);
                y = randi([1, grid_size]);
                global_call_id = global_call_id + 1;
                c = struct('id', global_call_id, 'x', x, 'y', y, ...
                    'dur', max(1, round(exprnd(call_duration))), ...
                    'move_timer', move_time_threshold, 'pool_type', 0, 'ch_idx', 0);
                
                if is_in_obs(x, y, obs_range), stat_total_arrivals = stat_total_arrivals + 1; end
                
                [fixed_pool, dynamic_pool, c, assigned, pool_type] = ...
                    assign_channel(c, fixed_pool, dynamic_pool, Nf, Nd);
                
                if assigned
                    active_calls(c.id) = c;
                    if pool_type == 2 && is_in_obs(x, y, obs_range), stat_dyn_attempts = stat_dyn_attempts + 1; end
                elseif is_in_obs(x, y, obs_range)
                    stat_blocked = stat_blocked + 1;
                end
            end
            
            % 3. 采样统计 (16x16 区域)
            if mod(t, 20) == 0 
                [on_calls, fixed_on] = count_calls_in_obs(active_calls, obs_range);
                stat_carried_samples = stat_carried_samples + on_calls;
                stat_carried_fixed_samples = stat_carried_fixed_samples + fixed_on;
                sample_counts = sample_counts + 1;
            end
        end
        
        % 指标计算
        num_obs_cells = length(obs_range)^2;
        if stat_total_arrivals > 0
            res_blocking_prob(f_idx, L_idx) = (stat_blocked / stat_total_arrivals) * 100;
            res_dyn_attempts_pct(f_idx, L_idx) = (stat_dyn_attempts / stat_total_arrivals) * 100;
            res_switched_pct(f_idx, L_idx) = (stat_switched / stat_total_arrivals) * 100;
        end
        if sample_counts > 0
            res_traffic_carried(f_idx, L_idx) = (stat_carried_samples / sample_counts) / num_obs_cells / 10;
            if Nf > 0
                res_fixed_carried(f_idx, L_idx) = (stat_carried_fixed_samples / sample_counts) / num_obs_cells / Nf;
            end
        end
    end
    fprintf('test.m仿真完成 Nf = %d... \n', Nf);
end

%% ==================== 绘图 ====================
figure('Position', [100, 100, 1000, 600]);
subplot(2,2,1); hold on; grid on;
for i=1:length(fixed_caps)
    plot(load_points, res_blocking_prob(i,:), '-o', 'DisplayName', sprintf('Nf=%d', fixed_caps(i)));
end
xlabel('Traffic Offered (Erlangs/cell)'); ylabel('Blocking (%)'); title('Fig 3: Blocking'); legend('Location','best');

subplot(2,2,2); hold on; grid on;
for i=1:length(fixed_caps)
    if fixed_caps(i) < 10
        plot(res_blocking_prob(i,:), res_switched_pct(i,:), '-d', 'DisplayName', sprintf('Nf=%d', fixed_caps(i)));
    end
end
xlabel('Blocking (%)'); ylabel('Switched (%)'); title('Fig 7: Switching Activity');

subplot(2,2,3); hold on; grid on;
for i=1:length(fixed_caps)
    if fixed_caps(i) < 10
        plot(res_blocking_prob(i,:), res_dyn_attempts_pct(i,:), '-s', 'DisplayName', sprintf('Nf=%d', fixed_caps(i)));
    end
end
xlabel('Blocking (%)'); ylabel('Dynamic Attempts (%)'); title('Fig 6: Dynamic Access');

subplot(2,2,4); hold on; grid on;
for i=1:length(fixed_caps)
    plot(res_blocking_prob(i,:), res_traffic_carried(i,:), '-x', 'DisplayName', sprintf('Nf=%d', fixed_caps(i)));
end
xlabel('Blocking (%)'); ylabel('Traffic Carried (E/chan)'); title('Fig 4: Carried Traffic');

%% ==================== 函数 ====================
function [bx, by, cx, cy] = get_cluster_coords(x, y)
    bx = ceil(x / 4); by = ceil(y / 4);
    cx = mod(x - 1, 4) + 1; cy = mod(y - 1, 4) + 1;
end

function in_obs = is_in_obs(x, y, obs_range)
    in_obs = any(obs_range == x) && any(obs_range == y);
end

function [fp, dp, c, ok, pt] = assign_channel(c, fp, dp, Nf, Nd)
    ok = false; pt = 0; [bx, by, cx, cy] = get_cluster_coords(c.x, c.y);
    if Nf > 0
        for i = 1:Nf
            if fp(bx, by, cx, cy, i) == 0
                fp(bx, by, cx, cy, i) = c.id; c.pool_type = 1; c.ch_idx = i;
                ok = true; pt = 1; return;
            end
        end
    end
    if Nd > 0
        for i = 1:Nd
            if dp(bx, by, i) == 0
                dp(bx, by, i) = c.id; c.pool_type = 2; c.ch_idx = i;
                ok = true; pt = 2; return;
            end
        end
    end
end

function [fp, dp, ac, sw] = release_and_reassign(c, fp, dp, ac, Nf)
    sw = false; [bx, by, cx, cy] = get_cluster_coords(c.x, c.y);
    if c.pool_type == 2
        dp(bx, by, c.ch_idx) = 0; remove(ac, c.id); return;
    end
    rel_idx = c.ch_idx; fp(bx, by, cx, cy, rel_idx) = 0; remove(ac, c.id);
    if Nf == 0, return; end
    
    % REassignment: 尝试从动态池拉回
    Nd = size(dp, 3);
    for i = 1:Nd
        d_id = dp(bx, by, i);
        if d_id ~= 0 && isKey(ac, d_id)
            dc = ac(d_id);
            if dc.x == c.x && dc.y == c.y
                fp(bx, by, cx, cy, rel_idx) = d_id; dp(bx, by, i) = 0;
                dc.pool_type = 1; dc.ch_idx = rel_idx; ac(d_id) = dc;
                sw = true; return;
            end
        end
    end
    % 否则向下压缩固定信道序列
    for i = Nf:-1:rel_idx+1
        mid = fp(bx, by, cx, cy, i);
        if mid ~= 0
            fp(bx, by, cx, cy, rel_idx) = mid; fp(bx, by, cx, cy, i) = 0;
            mc = ac(mid); mc.ch_idx = rel_idx; ac(mid) = mc; break;
        end
    end
end

function [total, fixed] = count_calls_in_obs(ac, obs)
    total = 0; fixed = 0; ids = keys(ac);
    for i = 1:length(ids)
        c = ac(ids{i});
        if is_in_obs(c.x, c.y, obs)
            total = total + 1;
            if c.pool_type == 1, fixed = fixed + 1; end
        end
    end
end