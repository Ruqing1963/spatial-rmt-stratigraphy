# 真实 Blue Lias 磁化率+岩性 (Weedon et al. 2018/2019, PANGAEA 896875, CC-BY-4.0)
# Southam Quarry, Warwickshire, 早侏罗世Hettangian-Sinemurian, Bucklandi带
# 格式: height_m, MAGS, lithology, bed
raw = """33.00 0.0345 Light_marl 44a
32.97 0.0395 Light_marl 44a
32.94 0.0225 Light_marl/Limestone 44a/43
32.91 0.0220 Limestone 43
32.88 0.0190 Limestone 43
32.85 0.0195 Limestone 43
32.82 0.0200 Limestone 43
32.79 0.0195 Limestone 43
32.76 0.0255 Limestone 43
32.73 0.0285 Light_marl 42
32.70 0.0365 Light_marl 42
32.67 0.0425 Light_marl 42
32.64 0.0315 Dark_marl 42
32.61 0.0355 Dark_marl 42
32.58 0.0330 Dark_marl 42
32.55 0.0335 Dark_marl 42
32.52 0.0370 Dark_marl 42
32.49 0.0310 Dark_marl 42
32.46 0.0305 Dark_marl 42
32.43 0.0335 Dark_marl 42
32.40 0.0365 Dark_marl 42
32.37 0.0340 Dark_marl 42
32.34 0.0325 Dark_marl 42
32.31 0.0415 Dark_marl 42
32.28 0.0340 Dark_marl 42
32.25 0.0350 Dark_marl 42
32.22 0.0400 Dark_marl 42
32.19 0.0395 Dark_marl 42
32.16 0.0405 Dark_marl 42
32.13 0.0340 Dark_marl 42
32.10 0.0395 Dark_marl 42
32.07 0.0435 Light_marl 42
32.04 0.0365 Light_marl 42
32.01 0.0410 Light_marl 42
31.98 0.0305 Light_marl 42
31.95 0.0380 Light_marl 42
31.92 0.0375 Light_marl 42
31.89 0.0340 Light_marl 42
31.86 0.0370 Light_marl 42
31.83 0.0360 Light_marl 42
31.80 0.0245 Limestone 41e
31.77 0.0215 Limestone 41e
31.74 0.0185 Limestone 41e
31.71 0.0180 Limestone 41e
31.68 0.0180 Limestone 41e
31.65 0.0190 Limestone 41e
31.62 0.0255 Limestone 41e
31.59 0.0280 Light_marl 41d
31.56 0.0320 Light_marl 41d
31.53 0.0300 Light_marl 41d
31.50 0.0280 Light_marl 41d
31.47 0.0315 Light_marl/Laminated_shale 41d
31.44 0.0410 Laminated_shale 41d
31.41 0.0340 Light_marl 41d
31.38 0.0190 Limestone 41c
31.35 0.0195 Limestone 41c
31.32 0.0205 Limestone 41c
31.29 0.0180 Limestone 41c
31.26 0.0185 Limestone 41c
31.23 0.0315 Dark_marl 41b
31.20 0.0320 Light_marl 41b
31.17 0.0210 Light_marl 41a
31.14 0.0170 Light_marl 41a
31.11 0.0190 Light_marl 41a
31.08 0.0235 Light_marl 41a
31.05 0.0205 Light_marl 41a
31.02 0.0330 Light_marl 40i
30.99 0.0345 Light_marl 40i
30.96 0.0360 Light_marl 40i
30.93 0.0320 Laminated_shale 40h
30.90 0.0370 Laminated_shale 40h
30.87 0.0500 Laminated_shale 40h
30.84 0.0410 Laminated_shale 40h
30.81 0.0420 Laminated_shale 40h
30.78 0.0360 Light_marl 40g
30.75 0.0335 Light_marl 40g
30.72 0.0415 Light_marl 40g
30.69 0.0315 Light_marl 40g
30.66 0.0320 Light_marl 40g
30.63 0.0420 Light_marl 40g
30.60 0.0490 Laminated_shale 40f
30.57 0.0480 Laminated_shale 40f
30.54 0.0400 Laminated_shale 40f
30.51 0.0455 Laminated_shale 40f
30.48 0.0325 Light_marl 40e
30.45 0.0335 Light_marl 40e
30.42 0.0225 Light_marl 40e
30.39 0.0210 Light_marl 40e
30.36 0.0225 Light_marl 40e
30.33 0.0285 Light_marl 40e
30.30 0.0295 Light_marl 40e
30.27 0.0290 Light_marl 40e
30.24 0.0270 Light_marl 40e"""

import numpy as np
h=[];chi=[];lith=[];bed=[]
for ln in raw.strip().split("\n"):
    p=ln.split()
    h.append(float(p[0])); chi.append(float(p[1])); lith.append(p[2]); bed.append(p[3])
h=np.array(h); chi=np.array(chi)
# 按高度升序(地层从下往上)
idx=np.argsort(h); h=h[idx]; chi=chi[idx]
lith=[lith[i] for i in idx]; bed=[bed[i] for i in idx]
np.savez('/home/claude/bluelias_real.npz', height=h, chi=chi,
         lith=np.array(lith), bed=np.array(bed))
print(f"Blue Lias 真实段: {len(h)}点, 高度 {h.min():.2f}-{h.max():.2f}m, 跨度{h.max()-h.min():.2f}m")
print(f"采样间隔: {np.median(np.diff(h)):.3f}m")
print(f"岩性类别: {set(lith)}")
# 统计灰岩层(旋回标志)
ls=[l for l in lith if 'Limestone' in l]
print(f"灰岩层数据点: {len(ls)}")
# 唯一bed数(每个bed是一个沉积单元)
print(f"地层单元(bed)数: {len(set(bed))}")
