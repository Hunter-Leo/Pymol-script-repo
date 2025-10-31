# PyMOL Scripts

PyMOL 可视化脚本集合，用于分子结构分析和可视化。

## 依赖

- PyMOL
- ImageMagick (用于生成 GIF)

**安装方法:**
```bash
conda activate pymol
conda install -c conda-forge imagemagick
```

## 脚本列表

### 01 - 分子可视化与 Free Cysteine 分析

**脚本:** `01/visualize_molecule.py`

**功能:**
- 多链着色显示（Mirabo 品牌配色）
- 标记所有 Cysteine 残基
- 识别并高亮 Free Cysteine（未形成二硫键）
- 多种可视化样式（cartoon, surface, sticks, spheres, lines）
- 显示侧链（可选）
- 生成多角度高清截图
- 生成平滑旋转动画并聚焦到 Free Cysteine

**快速开始:**
```bash
# 查看完整帮助文档
conda run -n pymol python pymol-scripts/01/visualize_molecule.py --help

# 基本用法（默认 cartoon 视图，低质量）
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file>

# 使用表面视图
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view surface

# 高质量渲染
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view surface --quality high
```

**使用示例:**
```bash
# 指定输出目录
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> ./output/

# 显示侧链
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> ./output/ --show-sidechain

# 不同视图样式
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view cartoon   # 卡通视图（默认）
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view surface   # 表面视图
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view sticks    # 棍状模型
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view spheres   # 球形模型
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --view lines     # 线条模型

# 指定渲染质量（仅影响 GIF 动画）
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --quality low     # 快速
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --quality medium  # 平衡
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> --quality high    # 高质量

# 组合使用
conda run -n pymol python pymol-scripts/01/visualize_molecule.py <pdb_file> ./output/ --view surface --quality high --show-sidechain
```

**参数说明:**
- `<pdb_file>` - 输入的 PDB/CIF 文件路径（必需）
- `[output_dir]` - 输出目录，默认为当前目录（可选）
- `--view STYLE` - 可视化样式（可选，默认 cartoon）
  - `cartoon`: 二级结构表示，显示 α-螺旋和 β-折叠，适合理解蛋白质整体折叠
  - `surface`: 分子表面表示，显示溶剂可及表面，适合分析结合位点和蛋白质界面
  - `sticks`: 棍状模型，显示所有原子和键，适合详细检查原子位置
  - `spheres`: 空间填充球体（CPK 模型），显示范德华半径，适合理解分子体积
  - `lines`: 简单线条表示，连接成键原子，适合快速概览
- `--quality LEVEL` - 渲染质量，仅影响 GIF 动画（可选，默认 low）
  - `low`: 600x450 px, 100 DPI, 无光线追踪，速度快，文件小
  - `medium`: 800x600 px, 150 DPI, 有光线追踪，平衡质量和速度
  - `high`: 1200x900 px, 200 DPI, 有光线追踪，发表级质量
  - 注意：静态 PNG 图片始终为 1920x1080 @ 300 DPI
- `--show-sidechain` - 显示氨基酸侧链（线条表示）（可选）

**输出文件:**
- `<molecule>_front.png` - 正面视图（1920x1080 @ 300 DPI, ray traced）
- `<molecule>_top.png` - 顶部视图（90° X 轴旋转）
- `<molecule>_side.png` - 侧面视图（90° Y 轴旋转）
- `<molecule>_overview.png` - 全景视图（缩放显示完整分子）
- `<molecule>_free_cys_N.png` - Free Cysteine 特写（如存在，带标签）
- `<molecule>_rotation.gif` - 旋转动画（360° 旋转 + 平滑 zoom in 到 Free Cysteine）
- `<molecule>_session.pse` - PyMOL session 文件（保留所有渲染设置）

**视图样式对比:**

| 样式 | 适用场景 | 优点 | 缺点 |
|------|---------|------|------|
| **cartoon** | 蛋白质整体结构分析、演示 | 清晰显示二级结构，易于理解折叠 | 不显示原子细节 |
| **surface** | 结合位点、蛋白质界面分析 | 显示表面特征和空腔 | 隐藏内部结构 |
| **sticks** | 原子级详细分析、小分子 | 显示所有原子和键 | 大分子时视觉混乱 |
| **spheres** | 分子体积、堆积分析 | 真实反映原子大小 | 遮挡严重，难以看清内部 |
| **lines** | 快速概览、大型复合物 | 视觉简洁，性能好 | 缺少细节信息 |

**推荐组合:**
- 结构展示：`--view cartoon`（默认）
- 表面分析：`--view surface --quality high`
- 详细检查：`--view sticks --show-sidechain`
- 发表图片：`--view cartoon --quality high` 或 `--view surface --quality high`


**交互式查看:**
```bash
# 运行脚本后，使用 session 文件打开 PyMOL（保留所有颜色和标记）
pymol ./output/<molecule>_session.pse
```

**Free Cysteine 检测:**
- 检测距离阈值：2.5 Å（S-S 原子间距离）
- 输出包含每个 Free Cysteine 到最近 SG 原子的距离
- 典型二硫键距离：2.0-2.1 Å

**颜色标记:**
- 不同链：Mirabo 品牌配色（深蓝、亮绿、橙、紫、青、金、洋红、浅蓝）
- 普通 Cysteine：黄色棍状
- Free Cysteine：红色棍状（加粗）
- 侧链：线状显示（可选）

**示例:**
```bash
# 快速预览（cartoon 视图）
conda run -n pymol python pymol-scripts/01/visualize_molecule.py fold_model.cif

# 表面视图分析结合位点
conda run -n pymol python pymol-scripts/01/visualize_molecule.py fold_model.cif ./output/ --view surface

# 高质量发表图片
conda run -n pymol python pymol-scripts/01/visualize_molecule.py fold_model.cif ./output/ --view cartoon --quality high

# 详细原子级分析
conda run -n pymol python pymol-scripts/01/visualize_molecule.py fold_model.cif ./output/ --view sticks --show-sidechain

# 完整分析流程
conda run -n pymol python pymol-scripts/01/visualize_molecule.py fold_model.cif ./output/ --view surface --quality high --show-sidechain
```

**工作流程建议:**
1. 首次分析：使用默认 cartoon 视图快速了解结构
2. 表面分析：使用 surface 视图识别结合位点和 Free Cysteine 位置
3. 详细检查：使用 sticks 视图检查感兴趣区域的原子细节
4. 最终输出：根据需求选择合适视图生成高质量图片
