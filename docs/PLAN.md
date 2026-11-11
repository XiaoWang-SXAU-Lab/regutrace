# ReguTrace：本地 TF–Target 证据工具（Phase 1）

## 目标

在不依赖 PlantTFDB 等在线服务的前提下，用本地的 genome FASTA、GFF3、TF–motif
映射和 PWM 库，构建可复现的 **TF–Target 候选表**。工具输出的是“候选 + 证据级别”，
不是体内 binding 的最终结论。

## 非目标（Phase 1 不做）

- 不从蛋白序列做 de novo motif 预测；蛋白序列只用于后续 TF 家族分类。
- 不做 GRN / 多层网络整合；这些在 TF–Target 证据层稳定后再做。
- 不把 prediction 和 ChIP/DAP/点对点验证混成一个分数。

## 证据级别（固定词表）

| level | 含义 | 例子 |
|---|---|---|
| A | 点对点实际验证 | EMSA、reporter、ChIP-qPCR、突变体 |
| B | 体内广谱 | ChIP-seq、ChIP-chip |
| C | 体外广谱 | DAP-seq、PBM、SELEX |
| D | 计算预测 | promoter PWM、motif、同源 |
| E | 功能/关联证据 | 表达相关、模块共表达、扰动效应 |

注意：DAP-seq 是 in vitro；ChIP-seq 是 in vivo。两者不能混标。

## 输入

- `--genome`：genome FASTA（可 `.gz`）
- `--gff`：GFF3（可 `.gz`）
- `--motifs`：JASPAR PFM 文本（多记录）
- `--tf-map`：可选 TSV，列 `tf_id, tf_symbol, family, motif_id`；省略时用
  motif header name 作为 `tf_id`/`tf_symbol`
- `--promoter-len`：上游长度，默认 1000 bp
- `--downstream`：TSS 下游长度，默认 0；PlantRegMap 使用 500 上游 + 100 下游
- `--pvalue`：PWM hit 的 p 值阈值，默认 1e-4
- `--background`：A,C,G,T 背景概率，默认 0.25,0.25,0.25,0.25
- `--annotation`：可选，TSV `gene_id, symbol`
- `--user-list`：可选，TSV `tf, target, evidence_type, evidence_level, method, source, pmid, condition, score, notes`
- `--out`：输出 TSV

## 输出（每条证据一行）

```text
tf_id  tf_symbol  target_id  target_symbol  evidence_type  evidence_level
method  source  pmid  condition  promoter_start  promoter_end  motif_id
hit_start  hit_end  strand  score  p_value  notes
```

- 预测记录：`evidence_type=prediction_motif`, `evidence_level=D`,
  `method=jaspar_pwm`, `source=promoter_scan`。
- 用户 list：`evidence_type=user_curated`（或用户指定）；`evidence_level` 必须显式
  给出且只能是 A/B/C/D/E，缺失或非法值会拒绝导入；预测字段留空。

## 处理流程

1. 读 genome FASTA、GFF3；按 transcript 选代表 TSS（每个 gene 取最长 transcript，
   缺失时退回 gene 坐标）。
2. 提取 promoter：默认取 TSS 上游 `--promoter-len`；可加 `--downstream` 纳入 TSS
   下游。`+` 链取上游；`-` 链取下游，并反向互补，使 promoter 始终按基因转录方向
   5'→3'。
3. 读 JASPAR PFM，转 log-odds；背景概率可配置。
4. 对每个 motif 构建 score 的 null 分布（DP，按 0.1 精度离散化），得到每个 hit 的
   p 值；扫描 promoter 双链。
5. 按 TF–motif 映射，输出 TF–Target 候选；可选合并用户 TF–Target list。
6. 输出 TSV；控制台打印 gene/motif/TF/hit 统计。

## 数据包（Phase 2）

- Arabidopsis TAIR10/Araport11 genome + GFF3（Ensembl Plants/TAIR）。
- JASPAR plant PFMs；后续可加 CIS-BP。
- 公共 DAP-seq/ChIP-seq 先整理成 TF–Target 边表（processed，不在安装包内塞原始
  bigWig/peak）。
- 用 DAP/ChIP 做 PWM 预测的 benchmark，按 TF family 校准阈值。
- 所有数据包记录版本、来源、许可、日期；不依赖在线服务。

## 阶段计划

- **Phase 1（本仓库）**：promoter + PWM + TF–Target 候选 + 用户 list 合并；
  synthetic tests；最小 CLI。
- **Phase 2**：Arabidopsis data pack；DAP/ChIP processed 边表；benchmark/阈值校准。
- **Phase 3**：接入 rsWGCNA worker/dashboard；按模块/条件展示证据；GRN 后续。

## 坐标与 p 值

- `promoter_start` / `promoter_end`：1-based inclusive 基因组坐标。
- `hit_start` / `hit_end`：1-based inclusive，相对于按基因方向 5'→3' 的 promoter。
- `p_value`：PWM 打分的离散化 null 分布（DP）得到，近似值；Phase 1 未做 FDR。
- `evidence_level=D` 表示纯预测；用户 list 保留用户给出的 `evidence_level`。

## Phase 1 状态（2026-11-11）

- 已实现：`regutrace` CLI（genome/GFF/JASPAR PFM/promoter 提取/PWM 扫描/
  用户 TF–target list 合并）。
- 测试：9 个 core 单元测试 + 1 个 synthetic CLI E2E 测试通过。
- 真实数据冒烟：Ensembl Plants release 60（27,655 genes）+ JASPAR 2024 plant
  PFMs；HSFA1E motif `MA1758.2` 扫描前 200 个 1 kb promoters，得到 67 个
  prediction hits，耗时约 3.4 s（含 genome/GFF 解析）。
- 当前缺口：`tf-map` 可省略并按 motif name 自动填充，但 name→gene ID/家族映射
  仍未做；未做 FDR、DAP/ChIP 边表和全基因组多 motif 性能优化。

## Benchmark（Biostrings 对照，2026-11-11）

脚本：`scripts/benchmark.R`。为了公平比较，两个工具都在**负样本 FPR
校准**后比较（target FPR = 1%），而不是直接比较“各自 max score 的 80%”。

合成数据：1,200 个 1 kb promoters；600 个正样本（300 个 consensus site + 300 个
弱 site，弱 site 分数为 max 的 50–85%），600 个负样本；motif = HSFA1E
`MA1758.2`；双链扫描。

| tool | selected threshold | TP | FP | FN | precision | recall | F1 | AP | time |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| regutrace | 21.08 (raw log2 odds) | 300 | 4 | 300 | 0.987 | 0.5 | 0.664 | 0.970 | **0.05 s** |
| Biostrings::matchPWM | 0.975 (unit scale) | 300 | 2 | 300 | 0.993 | 0.5 | 0.665 | 0.970 | 4.19 s |

- 两边都检出了全部 300 个 consensus site；300 个弱 site 在 1% FPR 阈值下都被
  漏掉，这是阈值灵敏度问题，不是工具差异。
- 在 matched FPR 下准确率几乎一致（AP 0.970 vs 0.970）；Rust CLI 约快 **80×**
  （0.05 s vs 4.19 s，1,200 promoters × 1 motif，含进程启动/文件解析）。
- 真实 1,000 个 Arabidopsis promoters：regutrace 4 hits，Biostrings 2 hits，
  Jaccard 0.5；耗时 0.05 s vs 2.06 s。两者在 1% FPR 下都非常严格。
- 注意：Biostrings 的 timing 包含 R 循环开销；下一阶段应加 FIMO/HOMER
  容器对照和 DAP-seq/ChIP-seq biological benchmark。

复现：

```bash
Rscript scripts/benchmark.R target/release/regutrace.exe \
  _rsWGCNA-local/regutrace-data/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt \
  MA1758.2 _rsWGCNA-local/regutrace-benchmark 1200 1000 \
  _rsWGCNA-local/regutrace-benchmark/real_promoters.fa
```

## 限制

- PWM 预测假阳性高；很多植物 TF 家族 motif 退火/共享，预测是 family-level。
- promoter 定义（长度、TSS、5' UTR、多转录本）会显著影响结果；输出必须带坐标。
- 没有 ATAC/ChIP 时，预测不能区分“可结合”和“实际结合”。
- 蛋白序列只做 TF 分类，不做 DNA 特异性预测。
