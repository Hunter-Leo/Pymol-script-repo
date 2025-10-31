请帮我写一个 pymol 的脚本，实现对分子的可视化，要求如下：
1. 需要标记分子的不同链，不同链用不同的颜色表示；
2. 需要标记 cystein 的位置，并特别标记出 free cystein；
3. 需要标记 side chain;
4. 每一次运行，给我几个标准的不同视角的截图，展示不同的角度和距离。如果存在 free cystein，需要额外增加一些图片 zoom in 到 free cystein 的位置；
5. 最后帮我生成一个动态旋转的分子视图，最后 zoom in 到 free cystein 的位置， 如果存在 free cystein， 并保存成高清的 gif 格式。
6. 脚本存放到 pymol-scripts/01/ 目录下， 一次运行只处理一个分子， 分子的文件名作为参数传入，并在该目录下生成对应的截图和 gif 文件。