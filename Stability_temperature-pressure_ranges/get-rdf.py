
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
PT = '0/traj-3000000-150'
filename = F"./su7ultra-relax/{PT}/npt-run-0-20000.lammpstrj"  # 修改为您的实际路径
skip_frames = 180  # 跳过的帧数（前N帧）
dr = 0.005          # 距离分辨率(Å)

def read_lammps_dump(filename, ref_type, target_type):
    """读取LAMMPS轨迹文件，处理三斜盒子"""
    frames = []
    with open(filename, 'r') as f:
        frame_count = 0
        while True:
            # 查找TIMESTEP行
            while True:
                line = f.readline()
                if not line:
                    return frames  # 文件结束
                if "ITEM: TIMESTEP" in line.strip():
                    break
                    
            # 读取时间步
            try:
                timestep = int(f.readline().strip())
            except:
                continue
                
            # 读取原子数
            while True:
                line = f.readline()
                if not line:
                    return frames
                if "ITEM: NUMBER OF ATOMS" in line.strip():
                    break
                    
            try:
                num_atoms = int(f.readline().strip())
            except:
                continue
                
            # 读取盒子边界 - 处理三斜盒子
            while True:
                line = f.readline()
                if not line:
                    return frames
                if "ITEM: BOX BOUNDS" in line.strip():
                    box_style = line.strip().split()[3:]  # 获取盒子类型
                    break
            
            # 读取三行盒子尺寸
            try:
                box_data = []
                for i in range(3):
                    box_line = f.readline().split()
                    if len(box_line) < 2:
                        break
                    box_data.append([float(x) for x in box_line])
                
                if len(box_data) != 3:
                    continue
                    
                # 提取盒子参数
                xlo, xhi = box_data[0][0], box_data[0][1]
                ylo, yhi = box_data[1][0], box_data[1][1]
                zlo, zhi = box_data[2][0], box_data[2][1]
                xy = box_data[0][2] if len(box_data[0]) > 2 else 0.0
                xz = box_data[1][2] if len(box_data[1]) > 2 else 0.0
                yz = box_data[2][2] if len(box_data[2]) > 2 else 0.0
                
                # 计算盒子尺寸
                Lx = xhi - xlo
                Ly = yhi - ylo
                Lz = zhi - zlo
                
                # 存储盒子参数
                box_params = {
                    'xlo': xlo, 'xhi': xhi,
                    'ylo': ylo, 'yhi': yhi,
                    'zlo': zlo, 'zhi': zhi,
                    'xy': xy, 'xz': xz, 'yz': yz,
                    'Lx': Lx, 'Ly': Ly, 'Lz': Lz,
                    'inv_Lx': 1.0/Lx, 'inv_Ly': 1.0/Ly, 'inv_Lz': 1.0/Lz
                }
                
            except Exception as e:
                print(f"读取盒子尺寸错误: {e}")
                continue
                
            # 读取原子数据头
            while True:
                line = f.readline()
                if not line:
                    return frames
                if "ITEM: ATOMS" in line.strip():
                    break
                    
            headers = line.strip().split()[2:]
            
            # 获取列索引
            id_idx = type_idx = xs_idx = ys_idx = zs_idx = None
            
            # 尝试不同可能的列名
            for i, h in enumerate(headers):
                if h in ['id']:
                    id_idx = i
                elif h in ['type']:
                    type_idx = i
                elif h in ['xs', 'x', 'xsu']:
                    xs_idx = i
                elif h in ['ys', 'y', 'ysu']:
                    ys_idx = i
                elif h in ['zs', 'z', 'zsu']:
                    zs_idx = i
            
            # 检查必要列
            if None in (id_idx, type_idx, xs_idx, ys_idx, zs_idx):
                # 打印列名以帮助调试
                print(f"警告: 缺少必要的原子数据列: 找到的列 {headers}")
                print(f"当前列: id_idx={id_idx}, type_idx={type_idx}, xs_idx={xs_idx}, ys_idx={ys_idx}, zs_idx={zs_idx}")
                continue
            
            # 读取原子数据
            atoms = []
            for _ in range(num_atoms):
                data = f.readline().split()
                if len(data) < max(id_idx, type_idx, xs_idx, ys_idx, zs_idx) + 1:
                    continue
                try:
                    atom_id = int(data[id_idx])
                    atom_type = int(data[type_idx])
                    xs = float(data[xs_idx])
                    ys = float(data[ys_idx])
                    zs = float(data[zs_idx])
                    
                    # 将缩放坐标转换为绝对坐标
                    x = xlo + xs * Lx + ys * xy + zs * xz
                    y = ylo + ys * Ly + zs * yz
                    z = zlo + zs * Lz
                    
                    atoms.append((atom_id, atom_type, x, y, z))
                except Exception as e:
                    print(f"原子数据转换错误: {e}")
                    continue
            
            # 只保留指定类型的原子
            selected_atoms = [a for a in atoms if a[1] in (ref_type, target_type)]
            if not selected_atoms:
                continue
                
            # 按原子ID排序
            selected_atoms.sort(key=lambda x: x[0])
            coords = np.array([(x, y, z) for _, _, x, y, z in selected_atoms])
            types = np.array([atom_type for _, atom_type, _, _, _ in selected_atoms])
            
            frames.append({
                'timestep': timestep,
                'box_params': box_params,
                'coords': coords,
                'types': types,
                'ref_mask': (types == ref_type),
                'target_mask': (types == target_type)
            })
            
            frame_count += 1
            if frame_count % 100 == 0:
                print(f"已读取 {frame_count} 帧...")
    
    return frames

def vectorized_minimum_image_distance(ref_coords, target_coords, box_params):
    """向量化计算三斜盒子中的最小镜像距离"""
    # 计算位移向量
    dx = target_coords[:, 0] - ref_coords[:, 0, None]
    dy = target_coords[:, 1] - ref_coords[:, 1, None]
    dz = target_coords[:, 2] - ref_coords[:, 2, None]
    
    # 获取盒子参数
    Lx, Ly, Lz = box_params['Lx'], box_params['Ly'], box_params['Lz']
    xy, xz, yz = box_params['xy'], box_params['xz'], box_params['yz']
    inv_Lx, inv_Ly, inv_Lz = box_params['inv_Lx'], box_params['inv_Ly'], box_params['inv_Lz']
    
    # 应用倾斜因子修正
    if xy != 0:
        dy -= np.round(dx * xy * inv_Lx) * Lx / xy
    if xz != 0:
        dz -= np.round(dx * xz * inv_Lx) * Lx / xz
    if yz != 0:
        dz -= np.round(dy * yz * inv_Ly) * Ly / yz
    
    # 应用周期性边界条件
    dx -= np.round(dx * inv_Lx) * Lx
    dy -= np.round(dy * inv_Ly) * Ly
    dz -= np.round(dz * inv_Lz) * Lz
    
    # 计算距离
    return np.sqrt(dx**2 + dy**2 + dz**2)

def calculate_rdf(frames, r_max, dr, ref_type, target_type, skip_frames=0):
    """计算径向分布函数 - 优化版本"""
    if not frames:
        raise ValueError("没有有效的帧可用于计算RDF")
    
    # 初始化直方图
    r_min = 0.0
    nbins = int((r_max - r_min) / dr)
    hist = np.zeros(nbins, dtype=float)
    norm = 0.0  # 归一化因子
    
    total_frames = len(frames)
    
    # 跳过前N帧
    if skip_frames > len(frames):
        skip_frames = len(frames) - 1
        print(f"警告: 跳过帧数超过总帧数，调整为跳过 {skip_frames} 帧")
    
    # 只使用跳过后的帧
    valid_frames = frames[skip_frames:]
    num_valid_frames = len(valid_frames)
    
    if num_valid_frames == 0:
        raise ValueError(f"跳过 {skip_frames} 帧后没有剩余帧可用于计算！")
    
    print(f"开始处理 {num_valid_frames} 帧数据 (跳过前 {skip_frames} 帧)...")
    
    # 计算抽样率 - 基于系统大小动态调整
    sample_size = 5000  # 最多处理5000帧
    sample_rate = max(1, num_valid_frames // sample_size)
    print(f"抽样率: 每 {sample_rate} 帧处理1帧")
    
    # 处理每一帧
    processed_frames = 0
    n_ref_total = 0
    n_target_total = 0
    
    # 使用进度条
    progress = tqdm(total=num_valid_frames, desc="计算RDF")
    
    for frame_idx, frame in enumerate(valid_frames):
        # 更新进度条
        progress.update(1)
        
        # 抽样处理
        if frame_idx % sample_rate != 0:
            continue
            
        coords = frame['coords']
        box_params = frame['box_params']
        ref_mask = frame['ref_mask']
        target_mask = frame['target_mask']
        
        # 检查是否有参考原子和目标原子
        n_ref = np.sum(ref_mask)
        n_target = np.sum(target_mask)
        if n_ref == 0 or n_target == 0:
            continue
        
        # 获取参考原子和目标原子的坐标
        ref_coords = coords[ref_mask]
        target_coords = coords[target_mask]
        
        # 向量化计算距离
        dists = vectorized_minimum_image_distance(ref_coords, target_coords, box_params)
        
        # 过滤掉自身距离 (设为NaN)
        np.fill_diagonal(dists, np.nan)
        
        # 展平距离数组并过滤有效距离
        flat_dists = dists.flatten()
        valid_dists = flat_dists[(flat_dists > 1e-5) & (flat_dists < r_max) & ~np.isnan(flat_dists)]
        
        # 计算直方图
        if valid_dists.size > 0:
            bin_indices = (valid_dists / dr).astype(int)
            valid_bins = bin_indices[bin_indices < nbins]
            np.add.at(hist, valid_bins, 1)
        
        # 更新统计信息
        n_ref_total += n_ref
        n_target_total += n_target
        processed_frames += 1
        
    progress.close()
    
    if processed_frames == 0:
        raise ValueError("没有处理任何帧！")
    
    print(f"实际处理了 {processed_frames} 帧数据")
    print(f"平均参考原子数: {n_ref_total/processed_frames:.1f}")
    print(f"平均目标原子数: {n_target_total/processed_frames:.1f}")
    
    # 计算平均归一化因子
    avg_vol = np.mean([f['box_params']['Lx'] * f['box_params']['Ly'] * f['box_params']['Lz'] 
                      for f in frames])
    avg_norm = (n_ref_total * n_target_total) / (processed_frames * avg_vol)
    
    if avg_norm <= 0:
        raise ValueError("归一化因子无效，没有找到有效的原子对")
    
    # 计算径向分布函数
    r = np.linspace(dr/2, r_max - dr/2, nbins)  # 中心点位置
    shell_vol = 4 * np.pi * (r**2) * dr  # 球壳体积
    
    # 防止零除错误
    shell_vol[shell_vol == 0] = np.inf
    
    # 计算RDF
    rdf = hist / (shell_vol * avg_norm * sample_rate)
    
    return r, rdf

# ==============================
# 用户可修改参数
# ==============================
ref_type = 1      # 参考原子类型 (H)
target_type = 1   # 目标原子类型 (H)
r_max = 10.0      # 最大距离(Å)
# ==============================

# 检查文件是否存在
if not os.path.exists(filename):
    print(f"错误: 文件 '{filename}' 不存在！")
    exit()

print(f"开始处理文件: {filename}")

# 读取轨迹数据
print("正在读取轨迹文件...")
try:
    frames = read_lammps_dump(filename, ref_type, target_type)
    print(f"成功读取 {len(frames)} 帧数据")
    if len(frames) == 0:
        print("错误：没有找到有效的帧！")
        print("可能的原因:")
        print(f"1. 原子类型不匹配 (当前设置 ref_type={ref_type}, target_type={target_type})")
        print("2. 文件格式不正确")
        print("3. 文件为空或损坏")
        exit()
    
    # 检查第一帧的数据
    first_frame = frames[0]
    print(f"第一帧信息:")
    print(f"  时间步: {first_frame['timestep']}")
    box = first_frame['box_params']
    print(f"  盒子尺寸: Lx={box['Lx']:.3f}Å, Ly={box['Ly']:.3f}Å, Lz={box['Lz']:.3f}Å")
    print(f"  倾斜因子: xy={box['xy']:.3f}, xz={box['xz']:.3f}, yz={box['yz']:.3f}")
    print(f"  原子总数: {len(first_frame['coords'])}")
    print(f"  参考原子数: {np.sum(first_frame['ref_mask'])}")
    print(f"  目标原子数: {np.sum(first_frame['target_mask'])}")
    
    print(f"参考原子类型: {ref_type}, 目标原子类型: {target_type}")
    print(f"计算参数: r_max={r_max}Å, dr={dr}Å, skip_frames={skip_frames}")
    
except Exception as e:
    print(f"读取文件时出错: {e}")
    exit()

# 计算RDF
print("正在计算径向分布函数...")
try:
    r, rdf = calculate_rdf(frames, r_max, dr, ref_type, target_type, skip_frames)
    print("RDF计算完成！")
    
    # 保存结果
    np.savetxt(f"{PT}_rdf_result.dat", np.column_stack((r, rdf)), 
               header="Distance(Å) g(r)", fmt='%.6f')
    
    # 绘制图形
    plt.figure(figsize=(10, 6))
    plt.plot(r, rdf, 'b-', linewidth=2)
    plt.xlabel('Distance (Å)', fontsize=14)
    plt.ylabel('g(r)', fontsize=14)
    plt.title(f'Radial Distribution Function (Type {ref_type}-{target_type})', fontsize=16)
    plt.grid(alpha=0.3)
    plt.xlim(0, r_max)
    plt.savefig(f'{PT}_rdf_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"计算完成！结果已保存到 {PT}_rdf_result.dat 和 {PT}_rdf_plot.png")

except Exception as e:
    print(f"计算RDF时出错: {e}")
    print("可能的问题：")
    print("1. 轨迹文件中没有找到指定类型的原子")
    print("2. 盒子尺寸为零或无效")
    print("3. 距离参数设置不当")
    print("4. 内存不足")